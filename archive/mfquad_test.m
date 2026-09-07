sigparams = loadjson(['JSON',filesep,'realizations.json']);
nRealizations = sigparams.nRealizations;
addpath(sigparams.path2sensitivitydata)
inFileData = struct(...
    'nSamples',sigparams.sampFreq*sigparams.datalen,...
    'dataY',zeros(nRealizations,sigparams.sampFreq*sigparams.datalen),...
    'dsstPSD',[],...
    'dsstfreqVec',[],...
    'tIntrvl',1/sigparams.sampFreq);

for lpruns = 1:nRealizations
    [inFileData.dataY(lpruns,:),inFileData.dsstPSD,inFileData.dsstfreqVec] = LIGOnoise(inFileData.nSamples,sigparams.sampFreq,1,'sample');
end

close all
%Design Sensitivity PSD plot
figure;
semilogy(inFileData.dsstfreqVec,inFileData.dsstPSD)
title('Design Sensitivity PSD')
axis tight
xlabel('Frequency (Hz)')
ylabel('Amplitude Spectrum (1/sqrt(Hz))')
%% Pre-Allocation and Data structure creation
psoParams = loadjson('JSON\pso.json');
signalParams = loadjson('JSON\signal.json');
fileName =  'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\SCRATCH\mfgwparam_test';
params = gwpsoparams(psoParams,signalParams,fileName);
inData = struct(...
    'design',cell(nRealizations,1),...
    'cond',cell(nRealizations,1),...
    'params',cell(nRealizations,1));
[inData(:).design] = deal(inFileData);
% [inData(:).params] = deal(params);
% mfData = zeros(nRealizations,2);
tempinData = inData;

negFStrt = 1-mod(inFileData.nSamples,2);
kNyq = floor(inFileData.nSamples/2)+1;
% Compute two-sided PSD from design sensitivity PSD for signal injection
dsstPSDtotal = [inFileData.dsstPSD, inFileData.dsstPSD((kNyq-negFStrt):-1:2)];
params.signal.data = sigInj(params,dsstPSDtotal);
mf1 = zeros(nRealizations,size(dsstPSDtotal,2));
mf2 = zeros(nRealizations,size(dsstPSDtotal,2));
maxMF = zeros(nRealizations,1);
%% Realization Loop
for lpruns = 1:nRealizations
    tempinData(lpruns).design.dataY = inFileData.dataY(lpruns,:);
    %Injected signal in strain
    [tempinData(lpruns).cond] = load_mfdata(tempinData(lpruns).design,[],params.signal); 
    %No injected signal
    % tempinData(lpruns).cond = load_mfdata(tempinData(lpruns).design,[]);
    
    tempinData(lpruns).cond.interpPSD = createPSD(tempinData(lpruns).cond.PSD,...
        tempinData(lpruns).cond.freqVec,...
        tempinData(lpruns).cond.tlen,...
        tempinData(lpruns).cond.sampFreq);

    %Injected signal in before whitening data
    tempinData(lpruns).params = cond_mfdata(tempinData(lpruns).cond,fileName,1);
    %No injected signal
       % tempinData(lpruns).params = cond_mfdata(tempinData(lpruns).cond,'mfgwparam_test');
    fftq0 = gen2PNwaveform(tempinData(lpruns).params,0,0,tempinData(lpruns).params.gwCoefs,1);
    fftq1 = fftq0.*tempinData(lpruns).params.phaseDiff;
    
    %Compute fitness value after maximizing by matched filtering
    mf1(lpruns,:) = matchedfiltering(tempinData(lpruns).params.fftdataYbyPSD, fftq0);
    mf2(lpruns,:) = matchedfiltering(tempinData(lpruns).params.fftdataYbyPSD, fftq1);
    maxMF(lpruns,:) = max(sqrt(mf1(lpruns,:).^2+mf2(lpruns,:).^2));
    disp(['Peak:', num2str(maxMF(lpruns))])
end

%% Plot analysis of Results
figure;
plot(params.dataX,params.signal.data)
title('Injected Signal')
axis tight
xlabel('Time(s)')

figure; 
plot(params.dataX,sqrt(mf1.^2+mf2.^2))
title('Matched Filtering Quadrature')
axis tight
xlabel('Time(s)')

figure;
histogram(maxMF,nRealizations)
title('Estimated Signal SNR Histogram')
axis tight
xlabel('Standard Deviation')
ylabel('SNR')
% medMF = median(sqrt(mf1.^2+mf2.^2));
% meanMF = mean(sqrt(mf1.^2+mf2.^2));
% 
