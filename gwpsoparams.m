% rungwpso setup script
function varargout = gwpsoparams(psoParams,signalParams,outFileName)
% Read JSON Files 
% jobParams = loadjson(jobParamsFile);
% psoParams = loadjson(jobParams.psoParamsjson);

% if jobParams.injSig == 1
%     signalParams = loadjson(jobParams.signalParamsjson);
% end

% Parameter Structure
T_sig_len = signalParams.signal.T_sig_len; %Length of signal in seconds
sampFreq = signalParams.sampling_freq; %Sampling Frequency, samples per second
nSamples = T_sig_len*sampFreq; %Number of samples
% nRuns = psoParams.nRuns; %Number of independent PSO runs

if psoParams.type == 2 % Search range of phase coefficients
    rmin = [signalParams.rmin_tau(1), signalParams.rmin_tau(2)];
    rmax = [signalParams.rmax_tau(1), signalParams.rmax_tau(2)];
        disp("Tau Space PSO");
else
    rmin = [signalParams.rmin(1), signalParams.rmin(2)];
    rmax = [signalParams.rmax(1), signalParams.rmax(2)];
        disp("Mass Space PSO");
end

% Constants
c = 299792458; %Speed of light in m/s
Msolar = 1.989*10^30; %Solar mass in kg
G = 6.6743*10^-11; %Gravitional constant in m^3/(kg*s^2)

cg = c^3/G;
% Tau coeffs as phase parameters
fmin = signalParams.freq(1); %Min cutoff Frequency 
m1 = signalParams.masses(1); %Mass of object 1 from signal in solar masses
m2 = signalParams.masses(2); %Mass of object 2 from signal in solar masses
if psoParams.type == 2
    %Mass conversion to kg
    m1 = m1*Msolar;
    m2 = m2*Msolar;

    M = m1 + m2;
    mu = m1*m2/(m1 + m2);
    eta = mu/M;
    gwCoefs(1) = (5/(256*pi))*(1/fmin)*((M*pi*fmin/cg)^(-5/3))*(1/eta);
    gwCoefs(2) = (1/8)*(1/fmin)*((M*pi*fmin/cg)^(-2/3))*(1/eta);
end
%% Pre-processing
%Make Frequency Magnitude and Phase difference vectors
phaseParams = phaseCalc(signalParams);

% %% Data Load
% % load(jobParams.outFilePSD,"PSD"); %Either SHAPES or pwelch
% load(jobParams.outFileData,'fftdataYbyPSD','TFtotal','dataY'); 
% 
% %% Create entire PSD vector 
% % negFStrt = 1-mod(nSamples,2);
% % kNyq = floor(nSamples/2)+1;
% % 
% % TFtotal = [TF, TF((kNyq-negFStrt):-1:2)];
% AbysqrtPSD = phaseParams.A.*TFtotal;
% 
% %% Create General Normalization Factor
% % Scalar factor of 1/nSamples is due to Parseval's theorem
% innProd = (1/nSamples)*(AbysqrtPSD)*AbysqrtPSD';
% genNormfacSqr = real(innProd);
% genNormfac = 1/sqrt(genNormfacSqr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optional: Generate and Inject custom CBC signal
% dataY = cbcsiginj(psoParams,signalParams);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timeVec = (0:nSamples-1)*(1/sampFreq);
%% Input Parameters:
params = struct('dataX', timeVec,...
                  'fpos', (0:floor(nSamples/2))*(1/T_sig_len),... %Positive frequency vector
                  'dataY', [],... %must be row vector
                  'fftdataYbyPSD', [],... %must be row vector
                  'frange', [fmin,signalParams.freq(2);],...
                  'datalen',T_sig_len,...,
                  'initial_phase', 0,...
                  'N', nSamples,...
                  'A', phaseParams.A,...
                  'phaseDiff', phaseParams.phaseDiff,...
                  'normfac', [],...
                  'avec', phaseParams.avec,...
		          'T_sig', signalParams.signal.T_sig,...
                  'rmin',rmin,...
                  'rmax',rmax,...
                  'Fs',sampFreq,...
                  'cgFac',cg,...
                  'Msolar',Msolar,...
                  'signal',signalParams,...
                  'pso',psoParams,...
                  'gwCoefs',gwCoefs);
if nargout>0
    varargout{1} = params;
end
save(outFileName,'params','psoParams','signalParams','gwCoefs')