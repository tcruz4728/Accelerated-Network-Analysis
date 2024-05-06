function varargout = gwpsoparams(psoParams,signalParams,varargin)
% GWPSOPARAMS(P,S)- computes parameters for use in rungwpso and saves
% params.mat file.
% P- PSO parameter structure containing the following fields
%   {
%   "type": "<Parameter space controller, =1 for mass space, =2 for tau space>",
%   "maxSteps": <Number of PSO iterations used in rungwpso's pso call>,
%   "nRuns": <Number of computing cores used in rungwpso's pso call>
%       }
% S- Signal parameter structure defining the parameters of the time series
% data or to define the parameters of an injected signal. Contains the
% following fields
%   {{
%   "sampling_freq": <sampling frequency of signal in Hz>,
%   "ta": <signal's time of arrival in seconds>,
%   "snr": <signal-to-noise ratio of time series data>,
%   "r":<distance between gw sources in parsecs>,
%   "phase": <phase difference of signal>,
%   "masses": <1x2 array defining the masses of gw sources>,
%   },
%   "signal": 
%     {
%       "T_sig_len": <total time series length in seconds>,
%       "T_sig": <length of signal in seconds>,
%       "num": 1,
%       "noise": 1
%       }
%   }
%Optional arguments
% D = GWPSOPARAMS(P,S,O)
% 
% O- Output filename for file containing 'params' structure, if unspecified
% file name will be 'params.mat'.
% 
% D- Returns 'params' structure.
% 
% 'params' structure has the following fields
% 
% 'dataX'-          time series integer vector of x data
% 'fpos'-           Positive frequency vector
% 'dataY'-          time series vector of y data, must be row vector.
% 'fftdataYbyPSD'-  FFT of data Y by PSD, computed by cond_mtchdfltrdata,
%                   must be row vector
% 'frange'-         Ranged of frequency values
% 'datalen'-        length of data in seconds
% 'initial_phase'-  initial phase of signal
% 'N'-              Total number of data samples in segment
% 'A'-              Computed amplitude magnitude frequency vector
% 'phaseDiff'-      Computed phase difference vector for quadrature
%                   templates
% 'normfac'-        Normalization factor
% 'avec'-           Matrix containing alpha terms used for waveform
%                   generation
% 'T_sig'-          Length of signal in seconds
% 'rmin'-           GW source minimum separation distance (from theory)
% 'rmax'-           GW source maximum separation distance (from theory)
% 'Fs'-             Sampling Frequency
% 'cgFac'-          Speed of light cubed and gravitational constant factor
% 'Msolar'-         Mass of Sun in kg
% 'signal'-         signalParams parameter structure
% 'pso'-            psoParams parameter structure
% 'gwCoefs'-        Gravitional wave coefficients, either masses or chirp
%                   times
%
% See also rungwpso, cond_mtchdfltrdata.

% Parameter Structure
T_sig_len = signalParams.signal.T_sig_len; %Length of signal in seconds
sampFreq = signalParams.sampling_freq; %Sampling Frequency, samples per second
nSamples = T_sig_len*sampFreq; %Number of samples

if psoParams.type == 2 % Search range of phase coefficients
    rmin = [signalParams.rmin_tau(1), signalParams.rmin_tau(2)];
    rmax = [signalParams.rmax_tau(1), signalParams.rmax_tau(2)];
        disp("gwpsoparams- Tau Space PSO");
else
    rmin = [signalParams.rmin(1), signalParams.rmin(2)];
    rmax = [signalParams.rmax(1), signalParams.rmax(2)];
        disp("gwpsoparams- Mass Space PSO");
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
                  'normfac', 1,...
                  'avec', phaseParams.avec,...
		          'T_sig', signalParams.signal.T_sig,...
                  'rmin',rmin,...
                  'rmax',rmax,...
                  'Fs',sampFreq,...
                  'cgFac',cg,...
                  'Msolar',Msolar,...
                  'signal',signalParams,...
                  'pso',psoParams,...
                  'gwCoefs',[gwCoefs,psoParams.type]);
%% Output
if nargout>0
    varargout{1} = params;
end
if nargin >2
    if ischar(varargin{1})
    save(varargin{1},'params','psoParams','signalParams','gwCoefs')
    else
        disp('gwpsoparams- No parameter file saved')
        return
    end
else 
    save('params.mat','params','psoParams','signalParams','gwCoefs')
end