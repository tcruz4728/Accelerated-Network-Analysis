function [] = genmultipsolnchrjb(path2jsonlab,jobParamsFile,varargin)
%GENMULTIPSOLNCHRJB(J,P)
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
%GENMULTIPSOLNCHRJB(J,P,U)
% U is an unique ID number that is used when creating a folder under the
% outfile directory. If specified, it uses that value with the format %05d,
% i.e. U = 1, out directory is <outDir>/00001. Can be used to overwrite
% preexisiting data(provided the date is the same) if an already used UID is
% specified.

%Modified from genlineshpslnchrjb, Aug 2023 for ANA use

%Add path to jsonlab
addpath(path2jsonlab);

%Load job parameters
jobParams = loadjson(jobParamsFile);

%Optionals for ana_basics setup
userUID = 1;
datad = [];
%Override the file name if optional input given
nreqArgs = 2;
for lpargs = 1:(nargin-nreqArgs)
    if ~isempty(varargin{lpargs})
        switch lpargs
            case 1
                userUID = varargin{lpargs};
            case 2
                datad = varargin{lpargs};
        end
    end
end

% progCtrl = 1; % Always set to on for debugging

%% ANA prep package
[paramsFile,outdataFilePrfx] = ana_basics(jobParams,userUID,datad,[],1);
paramsFileshps = [paramsFile,'shps'];
paramsFileList = cell(jobParams.inFileDataRange(2),1);
paramsFileshpsList = cell(jobParams.inFileDataRange(2),1);
dataFileList = cell(jobParams.inFileDataRange(2),1);
shpsDataFileList = cell(jobParams.inFileDataRange(2),1);
for fileCount = 1:jobParams.inFileDataRange(2)
    paramsFileList{fileCount} = [paramsFile,'_n',...
        num2str(fileCount),'.mat'];
    paramsFileshpsList{fileCount} = [paramsFileshps,'_n',...
        num2str(fileCount),'.mat'];
    dataFileList{fileCount} = [outdataFilePrfx,'_n',...
        num2str(fileCount),'C'];
    shpsDataFileList{fileCount} = [outdataFilePrfx,'_n',...
        num2str(fileCount),'shps_C'];
end
% fidShpsOutFileList = fopen([outdataFilePrfx,'shpsoutFilesList.txt'],'r');
%% Constuct job file for Launcher
fidJbFile = fopen([jobParams.scrtchDir,filesep,jobParams.jobName,'_jbfile.txt'],'w');
disp(['Job File: ',jobParams.jobName,'_jbfile.txt',' file created in ',jobParams.scrtchDir,filesep])
%Store list of output files in .txt file for post-processing codes
fidOutFileList = fopen([outdataFilePrfx,'_outFilesList.txt'],'w');
disp(['Output File list created: ',outdataFilePrfx,'_outFilesList.txt'])

nJobs = 1;
for nCount = 1:fileCount
    for runType = 1:2
        %  PSO and drase command on input file.
        fprintf(fidJbFile,'matlab -batch ');
        %path to jsonlab,
        fprintf(fidJbFile,' "addpath ''%s''; ', path2jsonlab);
        %path to DRASE
        fprintf(fidJbFile,' addpath ''%s''; ', jobParams.path2drase);
        %path to SHAPES, PSO, and project
        fprintf(fidJbFile,' setpath(''%s''); ', jobParamsFile);
        %Call to matched filtering code
        switch runType
            case 1 %pwelch run
                fprintf(fidJbFile, ' rungwpso(''%s'',''%s'');" \n', ...
                    paramsFileList{nCount},dataFileList{nCount});
            case 2 %shapes run
                fprintf(fidJbFile, ' rungwpso(''%s'',''%s'');" \n', ...
                    paramsFileshpsList{nCount},shpsDataFileList{nCount});           
        end
        % Count number of jobs
        nJobs = nJobs + 1;
    end
    fprintf(fidOutFileList,'%s  %s\n',...
        dataFileList{nCount},...
        shpsDataFileList{nCount});
        
end
% fclose(fidOutFileList);
fclose(fidJbFile);
%% Slurm file generation
% genslurm(jobParams,nJobs,anabasicsstr,1,0.5)
genslurm(jobParams,nJobs)
end