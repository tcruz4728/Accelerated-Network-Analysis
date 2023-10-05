% rungwpso test script
%Set UID for particular run or to overwrite a previous run
userUID = input("Type UID value: ");
datad = input("Date of Data (leave empty for current): ");
%Writes a progress file 
progCtrl = input("For a progress file, type 1: ");

%Job File 
jobParamsFile = 'C:\Users\Thomas Cruz\Documents\GitHub\Accelerated-Network-Analysis\JSON\PCmtchdfltrtest_Job_params.json';
%Load job parameters
jobParams = loadjson(jobParamsFile);
%% Runs necessary prep-functions for rungwpso
[paramsFile,outdataFilePrfx,filepaths,progressFile] = ana_basics(jobParams,userUID,datad,[],progCtrl,12);
fidprog = fopen(progressFile,'a');
proglines = struct('nd','done.',...
    'pr','\nPwelch GW PSO Run...',...
    'sr','\nSHAPES GW PSO Run...',...
    'p','\Post-processing...');
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
postprocessing(jobParams,dataFile,shpsDataFile,filepaths);
    progstatus(proglines.nd,fidprog,progCtrl)
fclose(fidprog);