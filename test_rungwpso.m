% rungwpso test script
%Set UID for particular run or to overwrite a previous run
userUID = input("Type UID value: ");
datad = input("Date of Data (leave empty for current): ");
%Writes a progress file 
progCtrl = input("For a progress file, type 1: ");
path2jsonlab = 'C:\Users\tcruz\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\JSONLab_ a toolbox to encode_decode JSON files\jsonlab-2.0';
addpath(path2jsonlab)
addpath("References\")
%Job File 
jobParamsFile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\JSON\mf_PC_Job_params.json';
%Load job parameters
jobParams = loadjson(jobParamsFile);
%% Runs necessary prep-functions for rungwpso
[paramsFile,outdataFilePrfx,filepaths,progressFile] = ana_basics(jobParams,userUID,datad,[],progCtrl,12);
fidprog = fopen(progressFile,'a');
proglines = struct('nd','done.',...
    'pr','Pwelch GW PSO Run...',...
    'sr','SHAPES GW PSO Run...',...
    'p','Post-processing...');
%% Run PSO on GW data 
dataFile = [outdataFilePrfx,'C'];
shpsDataFile = [outdataFilePrfx,'shps_C'];
paramsFileshps = [paramsFile,'shps'];
    progstatus(proglines.pr,fidprog,progCtrl)
rungwpso(paramsFile,dataFile) %pwelch 
    progstatus(proglines.nd,fidprog,progCtrl)
    progstatus(proglines.sr,fidprog,progCtrl)
rungwpso(paramsFileshps,shpsDataFile) %shapes estimate
    progstatus(proglines.nd,fidprog,progCtrl)
%% Post Processing
    progstatus(proglines.p,fidprog,progCtrl)
combFileName = comb_anashpsjb(path2jsonlab,jobParamsFile,userUID);
close all
postprocessing(combFileName,filepaths,jobParams.injSig);
    progstatus(proglines.nd,fidprog,progCtrl)
fclose(fidprog);