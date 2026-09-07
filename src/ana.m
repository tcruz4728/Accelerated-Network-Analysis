function results = ana(params_jb,filepaths,opts)

arguments
    params_jb                   (1,1) struct = struct()
    filepaths                   (1,1) struct = struct()
    opts.ParameterOverride      (1,1) struct = struct()
    opts.UseDelayedSegments    (1,1) logical = false
    opts.nRuns                  (1,1) double = 6
    opts.RunID                 (1,1) string  = ""
    opts.ProgressMonitoring    (1,1) logical = true
end
ProgressMonitoring = true;

%% Runs necessary prep-functions for rungwpso

[paramsFile,outdataFilePrfx,~,progressFile] = ana_basics(params_jb,...
    "ProgressMonitoring",true,...
    "UseLegacyFolders",{false,filepaths},...
    "EarlyReturn",[false,true]);
fidprog = fopen(progressFile,'a');
proglines = struct('nd','done.',...
    'pr','Pwelch GW PSO Run...',...
    'sr','SHAPES GW PSO Run...',...
    'p','Post-processing...');
%% Run PSO on GW data 
dataFile = [outdataFilePrfx,'C'];
shpsDataFile = [outdataFilePrfx,'shps_C'];
paramsFileshps = [paramsFile,'shps'];
    progstatus(proglines.pr,fidprog,ProgressMonitoring)
rungwpso(paramsFile,dataFile) %pwelch 
    progstatus(proglines.nd,fidprog,ProgressMonitoring)
    progstatus(proglines.sr,fidprog,ProgressMonitoring)
rungwpso(paramsFileshps,shpsDataFile) %shapes estimate
    progstatus(proglines.nd,fidprog,ProgressMonitoring)
%% Post Processing
results.params = paramsFile;
results.shps_params = paramsFileshps;
results.progress = progressFile;
fclose(fidprog);