params = loadjson(['JSON',filesep,'realizations.json']);
sigparams = params;
nRealizations = sigparams.nRealizations;
addpath(sigparams.path2sensitivitydata)
inFileData = struct(...
    'NSamples',sigparams.sampFreq*sigparams.datalen,...
    'dataY',zeros(nRealizations,sigparams.sampFreq*sigparams.datalen),...
    'dsstPSD',[],...
    'dsstfreqVec',[],...
    'tIntrvl',1/sigparams.sampFreq);

for lpruns = 1:nRealizations
    [inFileData.dataY(lpruns,:),inFileData.dsstPSD,inFileData.dsstfreqVec] = LIGOnoise(inFileData.NSamples,sigparams.sampFreq,1,'sample');
end
%% 
psoParams = loadjson('JSON\pso.json');
signalParams = loadjson('JSON\signal.json');
dsstPSD = inFileData.dsstPSD;

params = gwpsoparams(psoParams,signalParams);
negFStrt = 1-mod(params.N,2);
kNyq = floor(params.N/2)+1;
PSDtotal = [dsstPSD, dsstPSD((kNyq-negFStrt):-1:2)];
dsstTFtotal = 1./sqrt(PSDtotal);
% inData = struct(...
%     'design',cell(nRealizations,1),...
%     'cond',cell(nRealizations,1),...
%     'params',cell(nRealizations,1));

AbysqrtPSD = params.A.*dsstTFtotal;
innProd = (1/params.N)*(AbysqrtPSD)*AbysqrtPSD';
params.normfac = 1/sqrt(real(innProd));
signalParams = params.signal;
m1 = signalParams.masses(1);
m2 = signalParams.masses(2);

%Create Fourier Phase vector
% wavephase = gen2PNwaveform(fpos, ta, phase, fmin, fmax, m1,m2,datalen, initial_phase, snr, N, avec, normfac);
wavephase = gen2PNwaveform(params,signalParams.ta,0,[m1,m2,1],30);

%Create waveform vector in time domain
waveVec = ifft(params.A.*wavephase);

% %Normalized to unit 1
% N = length(PSDtotal);
% % normfac = 1/sqrt((1/N)*sum((fft(waveVec)./PSDtotal).*conj(fft(waveVec))));
% normfac = 1/sqrt((1/(N*signalParams.sampling_freq))*sum((fft(waveVec)./PSDtotal).*conj(fft(waveVec))));
% % normfac = 1/sqrt(innerproduct(waveVec,waveVec,PSDtotal));

% Create final signal
% signal = params.signal.snr*normfac*waveVec;


dataY = inFileData.dataY(1,:) + waveVec*sqrt(params.Fs);

rolloff = 0.5; %Roll-off in seconds
winfiltdata = dataY.*tukeywin(length(dataY), rolloff*params.Fs/params.N)';

fftfiltdata = fft(winfiltdata);

whtndfftfiltdata = fftfiltdata.*dsstTFtotal;

whtndfiltdata = ifft(whtndfftfiltdata);

whtndfiltdata = whtndfiltdata/sqrt(params.Fs);

fftdataYbyPSD = fft(whtndfiltdata).*dsstTFtotal.*params.A;

fwavepos = waveform(params,0,0,params.gwCoefs);

% fwavepos = waveform_tau(fpos,ta,phase,fmin,fmax,tau0,tau1p5,datalen,initial_phase, avec);

if mod(params.N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end

phaseq0 = [fwavepos, fwaveneg];
phaseq1 = phaseq0.*params.phaseDiff;

mf1 = ifft(fftdataYbyPSD.*conj(params.normfac*phaseq0));
mf2 = ifft(fftdataYbyPSD.*conj(params.normfac*phaseq1));

figure; 
plot(sqrt(mf1.^2+mf2.^2))
% wave = snr*inParams.normfac*fwave;

% [inData(:).design] = deal(inFileData);
% % [inData(:).params] = deal(params);
% % mfData = zeros(nRealizations,2);
% tempinData = inData;
% for lpruns = 1:nRealizations
%     tempinData(lpruns).design.dataY = inFileData.dataY(lpruns,:);
%     tempinData(lpruns).cond = load_mfdata(tempinData(lpruns).design,[],params);
%     % tempinData(lpruns).cond.interpPSD = createPSD(tempinData(lpruns).cond.PSD,...
%     %     tempinData(lpruns).cond.freqVec,...
%     %     tempinData(lpruns).cond.tlen,...
%     %     tempinData(lpruns).cond.sampFreq);
%     tempinData(lpruns).params = cond_mfdata(tempinData(lpruns).cond,'mfgwparam_test',1);
% 
%     fftq0 = gen2PNwaveform(tempinData(lpruns).params,0,0,tempinData(lpruns).params.gwCoefs,1);
%     fftq1 = fftq0.*tempinData(lpruns).params.phaseDiff;
% 
%     %Compute fitness value after maximizing by matched filtering
%     mf1 = matchedfiltering(tempinData(lpruns).params.fftdataYbyPSD, fftq0);
%     mf2 = matchedfiltering(tempinData(lpruns).params.fftdataYbyPSD, fftq1);
% end
