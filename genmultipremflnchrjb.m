function [outFileNameList] = genmultipremflnchrjb(path2jsonlab,jobParamsFile,varargin)
%GENMULTIPREMFLNCHRJB(J,P)
%Generates a .slurm file containing a LAUNCHER job for matched filtering
%data loading, conditioning, and (optional) estimating pwelch PSDs using
%SHAPES. Necessary to use before genmultimflnchrjb.m.
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
% "inFile":"<path to the file containing the pwelch estimated PSD
%       training data>",
% "outFile":"<path to the file containing SHAPES estimated PSD
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
%GENMULTIPREMFLNCHRJB(J,P,U)
% U is an unique ID number that is used when creating a folder under the
% output file directory. If specified, it uses that value with the format
% %05d, i.e. U = 1, out directory is <outDir>/00001. Can be used to
% overwrite preexisiting data(provided the date is the same) if an already
% used UID is specified.
% Optional Input Arguments
% S is the SHAPES control parameter which determines whether matched
% filtering Jobs will be created for a SHAPES estimated PSD in addition to
% the norm. S = 1 is the default which creates the SHAPES jobs, set to 0
% to disable SHAPES jobs.

% Modified from genlineshpslnchrjb, Aug 2023 for ANA use
% Modified name from genmultipremflnchrjb.m to better match function
% purpose, Jul 2025

addpath(path2jsonlab)

jobParams = loadjson(jobParamsFile);

userUID = 1;
shpsCtrl = 1;

%Override the file name if optional input given
nreqArgs = 2;
for lpargs = 1:(nargin-nreqArgs)
    if ~isempty(varargin{lpargs})
        switch lpargs
            case 1
                userUID = varargin{lpargs};
            case 2
                shpsCtrl = varargin{lpargs};
        end
    end
end

%% File Naming Convention
[paramsFile,outdataFilePrfx,filepaths] = ana_basics(jobParams,userUID,[],2);

[interFilePath,interFileName,~] = fileparts(jobParams.inFile);
[outFilePath,outFileName,~] = fileparts(jobParams.outFile);
inFileNameList = cell(jobParams.inFileDataRange(2),1);
interFileNameList = cell(jobParams.inFileDataRange(2),1);
outFileNameList = cell(jobParams.inFileDataRange(2),1);
paramsFileList = cell(jobParams.inFileDataRange(2),1);
if shpsCtrl == 1
    paramsFileshps = [paramsFile,'shps'];
    paramsFileshpsList = cell(jobParams.inFileDataRange(2),1);
end

for fileCount = jobParams.inFileDataRange(1):jobParams.inFileDataRange(2)
    %Time series initial file
    inFileNameList{fileCount} = [jobParams.inFileDataPrFx,...
        num2str(fileCount),'.mat'];
    %Pwelch data file (inFile-goes into SHAPES)
    interFileNameList{fileCount} = [interFilePath,filesep,...
        interFileName,'_n',num2str(fileCount),'.mat'];
    %Pwelch params files
    paramsFileList{fileCount} = [paramsFile,'_n',...
        num2str(fileCount),'.mat'];
    %Copying root file to others, calculated params do not change
    copyfile([paramsFile,'.mat'],paramsFileList{fileCount});
    %Create File lists for SHAPES files
    if shpsCtrl == 1
        %SHAPES estimated PSD data file (outFile)
        outFileNameList{fileCount} = [outFilePath,filesep,...
            outFileName,'_n',num2str(fileCount),'.mat'];
        %SHAPES params files
        paramsFileshpsList{fileCount} = [paramsFileshps,'_n',...
            num2str(fileCount),'.mat'];
        copyfile([paramsFileshps,'.mat'],paramsFileshpsList{fileCount});
    end
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
for nCount = jobParams.inFileDataRange(1):jobParams.inFileDataRange(2)
    for runType = 1:(1+shpsCtrl)
    fprintf(fidJbFile,'matlab -batch ');
    %path to jsonlab,
    fprintf(fidJbFile,' "addpath ''%s''; ', path2jsonlab);
    %path to DRASE
    fprintf(fidJbFile,' addpath ''%s''; ', jobParams.path2drase);
    %path to SHAPES, PSO, and project
    fprintf(fidJbFile,' setpath(''%s''); ', jobParamsFile);
    %Data Load
    fprintf(fidJbFile,' load_mfdata(''%s'',''%s'',''%s''); ',...
        inFileNameList{nCount},interFileNameList{nCount},num2str(jobParams.injSig));
    %PSD Interpolation
    fprintf(fidJbFile,' createPSD(''%s''); ',... welch
        interFileNameList{nCount});
    %Matched filtering conditioning
    fprintf(fidJbFile,' cond_mfdata(''%s'',''%s'');', ...
        interFileNameList{nCount},paramsFileList{nCount});
    fprintf(fidOutFileList,'%s  %s',...
        inFileNameList{nCount},interFileNameList{nCount});
            fprintf(fidparamsFileList,'%s',...
        paramsFileList{nCount});
    %Repeat process for SHAPES on PSD
    if shpsCtrl == 1
        %SHAPES call
        fprintf(fidJbFile,' drase4lines(''%s'',''%s'',''%s'',''%s'',''%s''); ',...
            jobParamsFile,[],filepaths.end,...
            interFileNameList{nCount},outFileNameList{nCount});
        fprintf(fidJbFile,' createPSD(''%s''); ',... shapes
            outFileNameList{nCount});
        %SHAPES conditioning
        fprintf(fidJbFile,' cond_mfdata(''%s'',''%s'');" \n', ...
            outFileNameList{nCount},paramsFileshpsList{nCount});
        fprintf(fidOutFileList,'  %s\n',...
            outFileNameList{nCount});
        fprintf(fidparamsFileList,'  %s\n',...
            paramsFileshpsList{nCount});
    else
        %Necessary next lines for file lists
        fprintf(fidJbFile,'" \n');
        fprintf(fidOutFileList,'\n');
        fprintf(fidparamsFileList,'\n');
    end
    nJobs = nJobs +1;
    end
end
fclose(fidJbFile);
fclose(fidparamsFileList);
fclose(fidOutFileList);
%% Slurm File Generation
genslurm(jobParams,nJobs)