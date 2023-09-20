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
if ~isempty(varargin{1})
    userUID = varargin{1};
end

%Add path to jsonlab
addpath(path2jsonlab);

%Load job parameters
jobParams = loadjson(jobParamsFile);
%ALINE and DRASE function and JSON load
addpath(genpath(jobParams.path2project));
addpath(genpath(jobParams.path2drase));
addpath(genpath(jobParams.path2shapes));

%Read job name: it will be used for output file names
jobName = jobParams.jobName;

%Defining Path and Folder Creation
[outDir,~] = dpfc(jobParams,userUID);
paramsFile = [outDir,'params'];
paramsFileshps = [paramsFile,'shps'];

% genDataParams = loadjson(jobParams.genDataParamsjson);
% draseParams = loadjson(jobParams.draseParamsjson);
psoParams = loadjson(jobParams.psoParamsjson);
signalParams = loadjson(jobParams.signalParamsjson);

%ReadMe Log File
% fidReadMe = fopen([outDir,'README.txt'],'a');
% fprintf(fidReadMe,'%s\n\n',datestr(datetime('now')));
% fclose(fidReadMe);

%Generate Prefix for output files
filetagstr = filetagana(psoParams,signalParams);
% lineDataFileName = [datascrtchDir,jobParams.jobName,'_',...
%     genDataParams.name,genDataParams.indCtrl];
% tmpdataFilePrfx = [datascrtchDir,jobName,'_',filetagstr];
outdataFilePrfx = [outDir,jobName,'_',filetagstr];

%% Load  data and split it, store files in output directory
% load(jobParams.inFile,'dataY');
% matData = line2mat(dataY,lineDataFileName,genDataParams,draseParams);
params = gwpsoparams(jobParamsFile,paramsFile);
copyfile([paramsFile,'.mat'],[paramsFileshps,'.mat']);

outData = load_mtchdfltrdata(jobParamsFile);
% run_test_drase4lines(jobParams,userUID);
createPSD(outData.PSD,outData.freqVec,outData.tlen,outData.sampFreq,jobParams.outFilePSD);
load(jobParams.inFileshpsPSD,"PSD")
createPSD(PSD,outData.freqVec,outData.tlen,outData.sampFreq,jobParams.outFileshpsPSD);

cond_mtchdfltrdata(jobParamsFile,params,paramsFile,...
    outData.sampFreq*outData.tlen,1);
cond_mtchdfltrdata(jobParamsFile,params, ...
    paramsFileshps,outData.sampFreq*outData.tlen,2);

%Constuct job file for Launcher
scrtchDir = jobParams.scrtchDir;
lnchrJbFName = [scrtchDir,filesep,jobName,'_jbfile.txt'];
fidJbFile = fopen(lnchrJbFName,'w');

%Store list of output files in .txt file for post-processing codes
fidOutFileList = fopen([outdataFilePrfx,'outFilesList.txt'],'w');
%Count number of jobs in the job file: can be different from
%nDataSplitFiles if this is a repeat run

nJobs = 2;
% [~,dataInFileName,~] = fileparts(datasplitFilesList);
for nCount = 1:2
    %     shpsOutFilesList{nCount} = [outDir,dataInFileName{nCount},'_C','.mat'];
    %     fprintf(fidOutFileList,'%s\n',shpsOutFilesList{nCount});
    %     %Check if the output file exists...
    %     if ~isempty(dir(shpsOutFilesList{nCount}))
    %         %File exists, so do not repeat the run
    %         continue;
    %     end
    %  PSO and drase command on input file.
    fprintf(fidJbFile,'matlab -batch ');
    %path to jsonlab,
    fprintf(fidJbFile,' "addpath ''%s''; ', path2jsonlab);
    %path to DRASE
    fprintf(fidJbFile,' addpath ''%s''; ', jobParams.path2drase);
    %path to SHAPES, PSO, and project
    %     fprintf(fidJbFile,' setpath(''%s''); ', jobParamsFile);
    %Call drasesetup which sets up parameters for drase run
    %     fprintf(fidJbFile,' rungwpsosetup(''%s'', ''%s'');" \n', ...
    %         paramsFile,outdataFilePrfx);
    switch nCount
        case 1
            fprintf(fidJbFile, ' rungwpso(''%s'',''%s'');" \n', ...
                [paramsFile,'.mat'],outdataFilePrfx); %pwelch ')
            fprintf(fidOutFileList,'%s\n',outdataFilePrfx);
        case 2
            fprintf(fidJbFile, ' rungwpso(''%s'',''%s'');" \n', ...
                [paramsFileshps,'.mat'],[outdataFilePrfx,'shps_']); %pwelch ')
            fprintf(fidOutFileList,'%s\n',[outdataFilePrfx,'shps_']);
    end
%     Count number of jobs
nJobs = nJobs + 1;
end
% fclose(fidOutFileList);
fclose(fidJbFile);

%% Generate SLURM file
nNodes = jobParams.nNodes;
qType = jobParams.qType;
if isempty(qType)
    qType = 'normal';
end
emailAdd = jobParams.email;

jbTime = jobParams.jbTime;
%Convert to hours, minutes, and seconds
tJobh = floor(jbTime);
jbTime = (jbTime-tJobh)*60;
tJobm = floor(jbTime);
jbTime = (jbTime - tJobm)*60;
tJobs = floor(jbTime);

slurmFileName = [scrtchDir,filesep,jobName,'.slurm'];
slurmFileId = fopen(slurmFileName,'w');
fprintf(slurmFileId,'#!/bin/bash\n');
fprintf(slurmFileId,'#-------------------------------------------------------\n');
fprintf(slurmFileId,'#\n');
fprintf(slurmFileId,'#SBATCH -J %s\n',jobName);
fprintf(slurmFileId,'#SBATCH -N %d\n',nNodes);
fprintf(slurmFileId,'#SBATCH -n %d\n',nJobs);
fprintf(slurmFileId,'#SBATCH -p %s\n',qType); 
fprintf(slurmFileId,'#SBATCH -o %s.o%%j\n',['TEMP',filesep,jobName]);
fprintf(slurmFileId,'#SBATCH -e %s.e%%j\n',['TEMP',filesep,jobName]);
fprintf(slurmFileId,'#SBATCH -t %d:%d:%d\n',tJobh,tJobm,tJobs);
fprintf(slurmFileId,'#SBATCH --mail-user=%s\n',emailAdd);
fprintf(slurmFileId,'#SBATCH --mail-type=all\n');
fprintf(slurmFileId,'#          <------ Account String ----->\n');
fprintf(slurmFileId,'# <--- (Use this ONLY if you have MULTIPLE accounts) --->\n');
fprintf(slurmFileId,'#SBATCH -A PHY21049\n');
% fprintf(slurmFileId,'#SBATCH -A GWEMpopStudy\n');
fprintf(slurmFileId,'#------------------------------------------------------\n');
fprintf(slurmFileId,'ml matlab\n');
fprintf(slurmFileId,'ml launcher\n');
fprintf(slurmFileId,'export LAUNCHER_WORKDIR=%s\n',scrtchDir);
fprintf(slurmFileId,'export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins\n');
fprintf(slurmFileId,'export LAUNCHER_RMI=SLURM\n');
fprintf(slurmFileId,'export LAUNCHER_JOB_FILE=%s\n\n',lnchrJbFName);
fprintf(slurmFileId,'$LAUNCHER_DIR/paramrun\n');
fclose(slurmFileId);

end