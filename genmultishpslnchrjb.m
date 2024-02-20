function [outFileNameList] = genmultishpslnchrjb(path2jsonlab,jobParamsFile,varargin)
%GENMULTISHPSLNCHRJB(J,P)
%Generates a .slurm file containing a LAUNCHER job for SHAPES estimated
%matched filtering. J is the path to the jsonlab package. If set to
%'', the jsonlab package is assumed to be in the Matlab search path. Note
%that this search path must be accessible from every compute node, which
%cannot always be guaranteed. Hence it is safer to specify the path
%explicitly.
%
%P is a JSON file containing the following job parameters. Text in <>
%should be replaced by an appropriate value. Examples are shown for some of
%the parameters.
%   {
% "jobName":"<a name for this job which will become the prefix for file
%            names>",
% "path2project":"<path to the project directory>",
% "path2drase":"<path to the DRASE directory>",
% "path2pso":"<path to the PSO directory>",
% "path2shapes":"<path to the SHAPES directory>",
% "inFilePSD":"<path to the file containing pwelch PSD training data>",
% "inFileshpsPSD":"<path to the file containing SHAPES estimated PSD
% training data>",
% "outFilePSD":"<path to the file containing pwelch PSD training data to be
% run by rungwpso>",
% "outFileshpsPSD":"<path to the file containing SHAPES estimated PSD training
% data to be run by rungwpso>",
% "inFileData":"<path to the file containing detector strain data>",
% "outDir":"<path to directory where all output files will be stored, a
% subdirectory for the date is created under this directory. A subdirectory
% under outDir/<date> with a UID>",
% "scrtchDir": <path to directory where all temporary files, e.g., SLURM
%               file, will be stored>",
% "genDataParamsjson":"<path to JSON file containing parameters for
%                generating line data>",
% "draseParamsFile":"<path to JSON file containing parameters for each
%                SHAPES run using drase>",
% "jbTime": <Time per Matlab job in hours>,
% "nNodes": <Number of compute nodes for this job>,
% "qType":"<job queue: leave empty for default/ls6,skx-normal for stampede2>",
% "email":"<email address>"
%   }
% }
% Optional Input Argument
%GENMULTISHPSLNCHRJB(J,P,U)
% U is an unique ID number that is used when creating a folder under the
% outfile directory. If specified, it uses that value with the format %05d,
% i.e. U = 1, out directory is <outDir>/00001. Can be used to overwrite
% preexisiting data(provided the date is the same) if an already used UID is
% specified.

%Modified from genlineshpslnchrjb, Aug 2023 for ANA use

% path2jsonlab = 'C:\Users\tcruz\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\JSONLab_ a toolbox to encode_decode JSON files\jsonlab-2.0';
% jobParamsFile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\JSON\multi_mtchdfltr_PC_Job_params.json';
addpath(path2jsonlab)

% userUID = input("UID:");
jobParams = loadjson(jobParamsFile);

%File Naming Convention
[paramsFile,outdataFilePrfx] = ana_basics(jobParams,varargin{1});
paramsFileshps = [paramsFile,'shps'];
[interFilePath,interFileName,~] = fileparts(jobParams.inFilePSD);
[outFilePath,outFileName,~] = fileparts(jobParams.inFileshpsPSD);
inFileNameList = cell(jobParams.inFileDataRange(2),1);
interFileNameList = cell(jobParams.inFileDataRange(2),1);
outFileNameList = cell(jobParams.inFileDataRange(2),1);
paramsFileList = cell(jobParams.inFileDataRange(2),1);
paramsFileshpsList = cell(jobParams.inFileDataRange(2),1);
for fileCount = 1:jobParams.inFileDataRange(2)
    %Time series initial file
    inFileNameList{fileCount} = [jobParams.inFileDataPrFx,...
        num2str(fileCount),'.mat'];
    %Pwelch data file (inFilePSD-goes into SHAPES)
    interFileNameList{fileCount} = [interFilePath,filesep,...
        interFileName,'_n',num2str(fileCount),'.mat'];
    %SHAPES estimated PSD data file (inFileshpsPSD)
    outFileNameList{fileCount} = [outFilePath,filesep,...
        outFileName,'_n',num2str(fileCount),'.mat'];
    %Pwelch params files
    paramsFileList{fileCount} = [paramsFile,'_n',...
        num2str(fileCount),'.mat'];
    %SHAPES params files
    paramsFileshpsList{fileCount} = [paramsFileshps,'_n',...
        num2str(fileCount),'.mat'];
    %Copying root file to others, calculated params do not change
    copyfile([paramsFile,'.mat'],paramsFileList{fileCount});
    copyfile([paramsFileshps,'.mat'],paramsFileshpsList{fileCount});
end

%% Job Text File
fidJbFile = fopen([jobParams.scrtchDir,filesep,jobParams.jobName,'_jbfile.txt'],'w');
disp(['Job File created: ',jobParams.scrtchDir,filesep,jobParams.jobName,'_jbfile.txt'])
%Store list of output files in .txt file for post-processing codes
fidOutFileList = fopen([outdataFilePrfx,'_inFilesList.txt'],'w');
disp(['Input File list created: ',outdataFilePrfx,'_inFilesList.txt'])
fidparamsFileList = fopen([outdataFilePrfx,'_paramsFilesList.txt'],'w');
disp(['Parameter File list created: ',outdataFilePrfx,'_paramsFilesList.txt'])

nJobs = 1;
for nCount = 1:fileCount
    fprintf(fidJbFile,'matlab -batch ');
    %path to jsonlab,
    fprintf(fidJbFile,' "addpath ''%s''; ', path2jsonlab);
    %path to DRASE
    fprintf(fidJbFile,' addpath ''%s''; ', jobParams.path2drase);
    %path to SHAPES, PSO, and project
    fprintf(fidJbFile,' setpath(''%s''); ', jobParamsFile);
    %Data Load
    fprintf(fidJbFile,' load_mtchdfltrdata(''%s'',''%s''); ',...
        inFileNameList{nCount},interFileNameList{nCount});
    %SHAPES call 
    fprintf(fidJbFile,' run_drase4lines(''%s'',''%s'',''%s'',''%s''); ',...
        jobParamsFile,outdataFilePrfx,interFileNameList{nCount},outFileNameList{nCount});
    %PSD Interpolation
    fprintf(fidJbFile,' createPSD(''%s''); ',... welch
        interFileNameList{nCount});
    fprintf(fidJbFile,' createPSD(''%s''); ',... shapes
        outFileNameList{nCount});
    %Matched filtering conditioning
    fprintf(fidJbFile,' cond_mtchdfltrdata(''%s'',''%s'',''%s'');', ...
        interFileNameList{nCount},paramsFileList{nCount},jobParams.injSig);
    %SHAPES conditioning
    fprintf(fidJbFile,' cond_mtchdfltrdata(''%s'',''%s'',''%s'');" \n', ...
        outFileNameList{nCount},paramsFileshpsList{nCount},jobParams.injSig);
    fprintf(fidOutFileList,'%s  %s  %s\n',...
        inFileNameList{nCount},interFileNameList{nCount},...
        outFileNameList{nCount});
    fprintf(fidparamsFileList,'%s  %s\n',...
        paramsFileList{nCount},paramsFileshpsList{nCount});
    nJobs = nJobs +1;
end
fclose(fidJbFile);
fclose(fidOutFileList);
%% Slurm File Generation
genslurm(jobParams,nJobs)