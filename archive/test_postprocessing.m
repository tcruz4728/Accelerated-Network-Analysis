% test_postprocessing
% This script tests the function postprocessing using 3 files, a job
% parameter json file, a data file containing the output data from rungwpso
% on pwelch data, and a data file containing the output data from rungwpso
% on shapes estimated pwelch data.

jobParamsFile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\JSON\multi_mf_PC_Job_params.json';
% dataFile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\SCRATCH\multi_mf_ls6_Job_tau_fs4096stp500tsL512ta138snr200_n1C.mat';
% dataFile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\SCRATCH\18-Nov-2024\00001\multi_ls6_Job_tau_fs4096stp500tsL512ta138snr30_F.mat
dataListFile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\SCRATCH\15-Jan-2026\00002\Dcomb_ls6_Job_tau_fs4096stp500tsL512ta138snr30_outFilesList.txt'; 
% shpsDataFile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\SCRATCH\multi_mf_ls6_Job_tau_fs4096stp500tsL512ta138snr200_n1shps_C.mat';
if ~isempty(dataListFile)
    fid = fopen(dataListFile,'r');
    C = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);
    dataList = C{1};
    N = length(dataList);
    combData = cell(N,5);
    for i = 1:N
        S = load(dataList{i}, 'inputData','psdData','estpsdData','outData','estoutData');
        combData{i,1} = S.inputData;
        combData{i,2} = S.psdData;
        combData{i,3} = S.estpsdData;
        combData{i,4} = S.outData;
        combData{i,5} = S.estoutData;
    end
else
    load(dataFile)
end

jobParams = loadjson(jobParamsFile);

%%
% postprocessing(jobParams,dataFile,shpsDataFile)
[~,filepaths]=dpfc(jobParams,1);
close all
for realization = 1:length(combData)
    postprocessing(combData(realization,:),filepaths.end,jobParams.injSig,realization)
end

%% Histogram plots
est_amp          = zeros(1,N);
est_shps_amp     = zeros(1,N);
est_phase        = zeros(1,N);
est_shps_phase   = zeros(1,N);
est_time         = zeros(1,N);
est_shps_time    = zeros(1,N);
est_fitness      = zeros(1,N);
est_shps_fitness = zeros(1,N);

for j = 1:N
    est_amp(j)          = combData{j,4}.outStruct.bestAmp;
    est_shps_amp(j)     = combData{j,5}.outStruct.bestAmp;
    est_phase(j)        = combData{j,4}.outStruct.bestPhase;
    est_shps_phase(j)   = combData{j,5}.outStruct.bestPhase;
    est_time(j)         = combData{j,4}.outStruct.bestTime;
    est_shps_time(j)    = combData{j,5}.outStruct.bestTime;
    est_fitness(j)      = combData{j,4}.outStruct.bestFitness;
    est_shps_fitness(j) = combData{j,5}.outStruct.bestFitness;
end

figure(10)
tiledlayout(2,2)   % 2x2 grid [web:1]

% 1) Amplitude
nexttile
histogram(est_amp,      'FaceColor','b','FaceAlpha',0.5)
hold on
histogram(est_shps_amp, 'FaceColor','r','FaceAlpha',0.5)
hold off
title('Amplitude SNR')
xlabel('Amplitude')
ylabel('Count')
legend('Estimate','Shapes','Location','best')

% 2) Phase
nexttile
histogram(est_phase,      'FaceColor','b','FaceAlpha',0.5)
hold on
histogram(est_shps_phase, 'FaceColor','r','FaceAlpha',0.5)
hold off
title('Phase')
xlabel('Phase')
ylabel('Count')
legend('Estimate','Shapes','Location','best')

% 3) Time
nexttile
histogram(est_time,      'FaceColor','b','FaceAlpha',0.5)
hold on
histogram(est_shps_time, 'FaceColor','r','FaceAlpha',0.5)
hold off
title('Time')
xlabel('Time')
ylabel('Count')
legend('Estimate','Shapes','Location','best')

% 4) Fitness
nexttile
histogram(est_fitness,      'FaceColor','b','FaceAlpha',0.5)
hold on
histogram(est_shps_fitness, 'FaceColor','r','FaceAlpha',0.5)
hold off
title('Fitness')
xlabel('Fitness')
ylabel('Count')
legend('Estimate','Shapes','Location','best')


