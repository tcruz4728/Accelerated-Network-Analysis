function []=rungwpso(jobParamsFile,outFileName)
%% Script to run Chirp-time or Mass Space PSO on a user specified datafile and PSD. 

%Raghav Girgaonkar, Apr 2023
%Thomas Cruz, May 2023, updated for use with ALINE

% Read JSON Files 
jobParams = loadjson(jobParamsFile);
psoParams = loadjson(jobParams.psoParamsjson);
signalParams = loadjson(jobParams.signalParamsjson);

% Parameter Structure
T_sig_len = signalParams.signal.T_sig_len; %Length of signal in seconds
sampFreq = signalParams.sampling_freq; %Sampling Frequency, samples per second
nSamples = T_sig_len*sampFreq; %Number of samples
nRuns = psoParams.nRuns; %Number of independent PSO runs

if psoParams.type == "tau" % Search range of phase coefficients
    rmin = [signalParams.rmin_tau(1), signalParams.rmin_tau(2)];
    rmax = [signalParams.rmax_tau(1), signalParams.rmax_tau(2)];
        disp("Tau Space PSO");
else
    rmin = [signalParams.rmin(1), signalParams.rmin(2)];
    rmax = [signalParams.rmax(1), signalParams.rmax(2)];
        disp("Mass Space PSO");
end

% Constants
c = 3*10^8; %Speed of light in m/s
Msolar = 1.989*10^30; %Solar mass in kg
G = 6.6743*10^-11; %Gravitional constant in m^3/(kg*s^2)

% Tau coeffs as phase parameters
fmin = signalParams.freq(1); %Min cutoff Frequency 
m1 = signalParams.masses(1); %Mass of object 1 from signal in solar masses
m2 = signalParams.masses(2); %Mass of object 2 from signal in solar masses
if psoParams.type == "tau"
    %Mass conversion to kg
    m1 = m1*Msolar;
    m2 = m2*Msolar;

    M = m1 + m2;
    mu = m1*m2/(m1 + m2);
    eta = mu/M;
    tau0 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-5/3))*(1/eta);
    tau1p5 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3)^(-2/3))*(1/eta);
end
save([outFileName,'params'])
%% Pre-processing
%Make Frequency Magnitude and Phase difference vectors
phaseParams = phaseCalc(signalParams);

%% Data Load
% load(jobParams.outFilePSD,"PSD"); %Either SHAPES or pwelch
load(jobParams.outFileData,'fftdataYbyPSD','TFtotal','dataY'); 

%% Create entire PSD vector 
% negFStrt = 1-mod(nSamples,2);
% kNyq = floor(nSamples/2)+1;
% 
% TFtotal = [TF, TF((kNyq-negFStrt):-1:2)];
AbysqrtPSD = phaseParams.A.*TFtotal;

%% Create General Normalization Factor
% Scalar factor of 1/nSamples is due to Parseval's theorem
innProd = (1/nSamples)*(AbysqrtPSD)*AbysqrtPSD';
genNormfacSqr = real(innProd);
genNormfac = 1/sqrt(genNormfacSqr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optional: Generate and Inject custom CBC signal
% dataY = cbcsiginj(psoParams,signalParams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timeVec = (0:nSamples-1)*(1/sampFreq);
%% Input Parameters:
tauParams = struct('dataX', timeVec,...
                  'fpos', (0:floor(nSamples/2))*(1/T_sig_len),... %Positive frequency vector
                  'dataY', dataY',... %must be row vector
                  'fftdataYbyPSD', fftdataYbyPSD',... %must be row vector
                  'frange', [fmin,signalParams.freq(2);],...
                  'datalen',T_sig_len,...,
                  'initial_phase', 0,...
                  'N', nSamples,...
                  'A', phaseParams.A,...
                  'phaseDiff', phaseParams.phaseDiff,...
                  'normfac', genNormfac,...
                  'avec', phaseParams.avec,...
		          'T_sig', signalParams.signal.T_sig,...
                  'rmin',rmin,...
                  'rmax',rmax,...
                  'Fs',sampFreq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maxSteps = psoParams.maxSteps;
if psoParams.type == "tau"
    original_fitVal = -1*mfqc_tau([tau0, tau1p5], tauParams);
    outStruct = crcbqcpso_tau(tauParams,psoParams,nRuns);
else
    original_fitVal = -1*mfqc([m1, m2], tauParams);
    outStruct = crcbqcpso_mass(tauParams,psoParams,nRuns,sampFreq);
end

bestFitVal = -1*outStruct.bestFitness;
%% Optional: Uncomment to save output of PSO in a .mat file
save([outFileName,'C'],'outStruct','bestFitVal',"original_fitVal",'tauParams');
end
