% test_rungwpso_dataDL
% rungwpso test script for data run on ls6 to check mf portion

path2jsonlab = 'C:\Users\tcruz\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\JSONLab_ a toolbox to encode_decode JSON files\jsonlab-2.0';
addpath(path2jsonlab)
addpath("References\")
%Job File 
jobParamsFile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\JSON\multi_shps_PC_Job_params.json';
%Load job parameters
jobParams = loadjson(jobParamsFile);
setpath(jobParams)
%% Update parameters in paramsFiles
paramsFile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\SCRATCH\23-Mar-2024\00001\params_n1.mat';
paramsFileshps = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\SCRATCH\23-Mar-2024\00001\paramsshps_n1';
filepaths.end = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\SCRATCH\23-Mar-2024\00001';
%% Run PSO on GW data 
dataFile = ['test','C'];
shpsDataFile = ['test','shps_C'];

rungwpso(paramsFile,dataFile) %pwelch 
rungwpso(paramsFileshps,shpsDataFile) %shapes estimate
%% Post Processing
combFileName = comb_anashpsjb(path2jsonlab,jobParams,1,[],...
    ['SCRATCH\',dataFile],['SCRATCH\',shpsDataFile]);
close all
postprocessing(combFileName,filepaths,jobParams.injSig);
