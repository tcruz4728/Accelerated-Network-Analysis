% test_postprocessing
% This script tests the function postprocessing using 3 files, a job
% parameter json file, a data file containing the output data from rungwpso
% on pwelch data, and a data file containing the output data from rungwpso
% on shapes estimated pwelch data.

jobParamsFile = 'C:\Users\Thomas Cruz\Documents\GitHub\Accelerated-Network-Analysis\JSON\mtchdfltr_ls6test_Job_params.json';
dataFile = 'C:\Users\Thomas Cruz\Documents\GitHub\Accelerated-Network-Analysis\SCRATCH\mtchdfltr_ls6test_Job_tau_fs4096N500tsl256Tin138SNR10_C.mat';
shpsDataFile = 'C:\Users\Thomas Cruz\Documents\GitHub\Accelerated-Network-Analysis\SCRATCH\mtchdfltr_ls6test_Job_tau_fs4096N500tsl256Tin138SNR10_shps_C.mat';

jobParams = loadjson(jobParamsFile);
postprocessing(jobParams,dataFile,shpsDataFile)