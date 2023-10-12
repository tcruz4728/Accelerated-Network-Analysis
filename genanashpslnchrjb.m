function [] = genanashpslnchrjb(path2jsonlab,jobParamsFile,varargin)
%GENLINESHPSLNCHRJB(J,P)
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
%GENANASHPSLNCHRJB(J,P,U)
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
if ~isempty(varargin{1})
    userUID = varargin{1};
end
% progCtrl = 1; % Always set to on for debugging

%% ANA prep package
[paramsFile,outdataFilePrfx] = ana_basics(jobParams,userUID,[],1);
dataFile = [outdataFilePrfx,'C'];
shpsDataFile = [outdataFilePrfx,'shps_C'];
paramsFileshps = [paramsFile,'shps'];

%% Constuct job file for Launcher
fidJbFile = fopen([jobParams.scrtchDir,filesep,jobParams.jobName,'_jbfile.txt'],'w');
anabasicsstr = [jobParams.jobName,'anasetup'];
fidbasicsJbFile = fopen([jobParams.scrtchDir,filesep,anabasicsstr,'_jbfile.txt'],'w');
%Store list of output files in .txt file for post-processing codes
fidOutFileList = fopen([outdataFilePrfx,'outFilesList.txt'],'w');


nJobs = 1;
for nCount = 1:3
    fid = fidJbFile;
    if nCount == 1
        fid = fidbasicsJbFile;
    end
    %  PSO and drase command on input file.
    fprintf(fid,'matlab -batch ');
    %path to jsonlab,
    fprintf(fid,' "addpath ''%s''; ', path2jsonlab);
    %path to DRASE
    fprintf(fid,' addpath ''%s''; ', jobParams.path2drase);
    %path to SHAPES, PSO, and project

    %Call drasesetup which sets up parameters for drase run
    %     fprintf(fidJbFile,' rungwpsosetup(''%s'', ''%s'');" \n', ...
    %         paramsFile,outdataFilePrfx);
    switch nCount
        case 1
            fprintf(fid, ' [~,~]=ana_basics(''%s'',''%05d'',[],[],1);" \n', ...
                jobParamsFile,userUID);
            fprintf(fidOutFileList,'%s\n',paramsFile);
        case 2
            fprintf(fid,' setpath(''%s''); ', jobParamsFile);
            fprintf(fid, ' rungwpso(''%s'',''%s'');" \n', ...
                [paramsFile,'.mat'],dataFile); %pwelch ')
            fprintf(fidOutFileList,'%s\n',dataFile);
        case 3
            fprintf(fid,' setpath(''%s''); ', jobParamsFile);
            fprintf(fid, ' rungwpso(''%s'',''%s'');" \n', ...
                [paramsFileshps,'.mat'],shpsDataFile); %shps on pwelch ')
            fprintf(fidOutFileList,'%s\n',shpsDataFile);
    end
    %     Count number of jobs
    nJobs = nJobs + 1;
end
% fclose(fidOutFileList);
fclose(fidJbFile);
%% Slurm file generation
genslurm(jobParams,nJobs,anabasicsstr,1,0.5)
genslurm(jobParams,nJobs)
end