%ANA define path and folder creation
%define and set paths for gendatashpslnchrjb.m and combsplitdatashps.m for
%PC use.
% clearvars
% close all
%Folder Pathing
% filepaths.root = 'C:\Users\Thomas Cruz\Documents\GitHub';
% jobParamsFile = [filepaths.root,'\Accelerated-Network-Analysis\JSON\PCmtchdfltrtest_Job_params.json'];
addpath(['JSON',filesep])
jobParams = loadjson(jobParamsFile);
% setpath(jobParamsFile);

filepaths.figs = [jobParams.path2project,filesep,'Figures',filesep,date,filesep];
filepaths.psd = [jobParams.path2project,filesep,'TMPPSDDATA',filesep];
filepaths.end = [jobParams.path2project,filesep,'TMPMTCHDDATA',filesep,date,filesep];
outDir = [jobParams.outDir,filesep,date,filesep];

%[Optional] UID Folder Creation and output file directory update
if isempty(userUID)
    userUID = input("Type UID value: ");
end

% Data Folder creation
cd(jobParams.path2project)
[folderstatus.linedata,~,~] = mkdir('TMPMTCHDDATA',date);
[folderstatus.fig,~,~] = mkdir('Figures',date);

outDir = mkdirUID(userUID,outDir);
filepaths.end = mkdirUID(userUID,filepaths.end);
paramsFile = [outDir,'params'];
paramsFileshps = [paramsFile,'shps'];