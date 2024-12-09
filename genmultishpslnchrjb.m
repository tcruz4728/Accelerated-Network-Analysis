function [outFileNameList] = genmultishpslnchrjb(path2jsonlab,jobParamsFile,varargin)
%GENMULTISHPSLNCHRJB(J,P)
%Generates a .slurm file containing a LAUNCHER job for matched filtering
%data loading, conditioning, and estimating pwelch PSDs using SHAPES.
%J is the path to the jsonlab package. If set to '', the jsonlab package is
%assumed to be in the Matlab search path. Note that this search path must
%be accessible from every compute node, which cannot always be guaranteed.
%Hence it is safer to specify the path explicitly.
%
%P is a JSON file containing the following job parameters. Text in <>
%should be replaced by an appropriate value. Examples are shown for some of
%the parameters.
%   {
% "jobName":"<a name for this job which will become the prefix for file
%            names>",
% "path2project":"<path to the project directory>",
% "path2drase":"<path to the DRASE directory>",
% "path2aline":"path to the ALINE directory>",
% "path2pso":"<path to the PSO directory>",
% "path2shapes":"<path to the SHAPES directory>",
% "inFileDataPrfx":"<file prefix and location of data realizaitons>",
% "inFileDataRange":"<1x2 array giving the range of values from which the
%       full file names of data realizations, i.e. [1, 5] will have 5 data
%       realizations loaded in and estimated.>",
% "inFilePSD":"<path to the file containing the pwelch estimated PSD
%       training data>",
% "inFileshpsPSD":"<path to the file containing SHAPES estimated PSD
%       training data>",
% "outDir":"<path to directory where all output files will be stored, a
%       subdirectory for the current date is created under this directory. 
%       Another subdirectory under outDir/<date> with a UID is also set.>",
% "scrtchDir": <path to directory where all temporary files, e.g., SLURM
%       files, job text files, will be stored>",
% "psoParamsjson":"<path to JSON file containing parameters for rungwpso's
%       PSO run.>",
% "signalParamsjson":"<path to JSON file containing parameters for a
%       realized signal or an injected signal>",
% "genDataParamsjson":"<path to JSON file containing parameters for
%       generating line data>",
% "draseParamsFile":"<path to JSON file containing parameters for each
%       SHAPES run using drase>",
% "injSig":"<Injection signal control parameter, if empty no signal
%       injection is performed, else the injected signal will have 
%       parameters set by signalParams.json>",
% "jbTime":"<Time per Matlab job in hours>",
% "nNodes":"<Number of compute nodes for this job>",
% "qType":"<job queue: leave empty for default/ls6,skx-normal for 
%       stampede2>",
% "email":"<email address>"
%   }
% }
% Optional Input Argument
%GENMULTISHPSLNCHRJB(J,P,U)
% U is an unique ID number that is used when creating a folder under the
% output file directory. If specified, it uses that value with the format
% %05d, i.e. U = 1, out directory is <outDir>/00001. Can be used to
% overwrite preexisiting data(provided the date is the same) if an already
% used UID is specified.

%Modified from genlineshpslnchrjb, Aug 2023 for ANA use

addpath(path2jsonlab)

jobParams = loadjson(jobParamsFile);

%% File Naming Convention
[paramsFile,outdataFilePrfx,filepaths] = ana_basics(jobParams,varargin{1},[],2);
paramsFileshps = [paramsFile,'shps'];
[interFilePath,interFileName,~] = fileparts(jobParams.inFilePSD);
[outFilePath,outFileName,~] = fileparts(jobParams.inFileshpsPSD);
inFileNameList = cell(jobParams.inFileDataRange(2),1);
interFileNameList = cell(jobParams.inFileDataRange(2),1);
outFileNameList = cell(jobParams.inFileDataRange(2),1);
paramsFileList = cell(jobParams.inFileDataRange(2),1);
paramsFileshpsList = cell(jobParams.inFileDataRange(2),1);
for fileCount = jobParams.inFileDataRange(1):jobParams.inFileDataRange(2)
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
    fprintf(fidJbFile,' load_mfdata(''%s'',''%s'',''%s''); ',...
        inFileNameList{nCount},interFileNameList{nCount},jobParams.injSig);
    %SHAPES call 
    fprintf(fidJbFile,' drase4lines(''%s'',''%s'',''%s'',''%s'',''%s''); ',...
        jobParamsFile,outdataFilePrfx,filepaths.end,...
        interFileNameList{nCount},outFileNameList{nCount});
    %PSD Interpolation
    fprintf(fidJbFile,' createPSD(''%s''); ',... welch
        interFileNameList{nCount});
    fprintf(fidJbFile,' createPSD(''%s''); ',... shapes
        outFileNameList{nCount});
    %Matched filtering conditioning
    fprintf(fidJbFile,' cond_mfdata(''%s'',''%s'');', ...
        interFileNameList{nCount},paramsFileList{nCount});
    %SHAPES conditioning
    fprintf(fidJbFile,' cond_mfdata(''%s'',''%s'');" \n', ...
        outFileNameList{nCount},paramsFileshpsList{nCount});
    fprintf(fidOutFileList,'%s  %s  %s\n',...
        inFileNameList{nCount},interFileNameList{nCount},...
        outFileNameList{nCount});
    fprintf(fidparamsFileList,'%s  %s\n',...
        paramsFileList{nCount},paramsFileshpsList{nCount});
    nJobs = nJobs +1;
end
fclose(fidJbFile);
fclose(fidparamsFileList);
fclose(fidOutFileList);
%% Slurm File Generation
genslurm(jobParams,nJobs)