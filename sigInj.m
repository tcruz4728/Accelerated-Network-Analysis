function [signal] = sigInj(params,PSDtotal)
% function [signal] = sigInj(fpos, ta, phase, fmin, fmax, m1,m2,r,datalen, initial_phase, N, avec, A)
%SIGINJ Returns time domain vector of waveform for injection with custom
%strain amplitude
%   This function creates a time domain vector of a custom injected
%   waveform normalized to a strain amplitude which is a function of the
%   component masses and the distance to the source. 
% This factor is given in Cutler and Flanagan 1994 for the newtonian waveform 
% (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.49.2658)
%  Input: fpos: Positive frequency vector, 
%         ta: time of arrival, 
%         phase: coalescence phase, 
%         [fmin, fmax]: waveform frequency bounds, 
%         {m1,m2,r}: the component masses and the distance in Mpc respectively, 
%         datalen: data length in seconds, 
%         initial_phase: initial phase, 
%         N: total number of samples, 
%         avec: Precalculated vectors for phase generation, 
%         A: frequency magnitude vector
% Output: signal: time-domain vector of waveform with normalized strain amplitude 

% load(paramsFile);
signalParams = params.signal;
m1 = signalParams.masses(1);
m2 = signalParams.masses(2);

%Create Fourier Phase vector
% wavephase = gen2PNwaveform(fpos, ta, phase, fmin, fmax, m1,m2,datalen, initial_phase, snr, N, avec, normfac);
wavephase = gen2PNwaveform(params,signalParams.ta,0,[m1,m2,1],1);

%Create waveform vector in time domain
waveVec = ifft(params.A.*wavephase);

%Normalized to unit 1
N = length(PSDtotal);
normfac = 1/sqrt((1/N)*sum((fft(waveVec)./PSDtotal).*conj(fft(waveVec))));
% normfac = 1/sqrt(innerproduct(waveVec,waveVec,PSDtotal));

% Create final signal
signal = params.signal.snr*normfac*waveVec;

% %Create waveform vector in time domain
% signal = ifft(wavefourier);
%Constants
% fmin = 30;
% c = 3*10^8;
% G = 6.6743*10^-11;
% Msolar = 1.989*10^30;
% Mpc = 3.8057*10^22; %1 Megaparsec in meters

% m1_val = m1*Msolar;
% m2_val = m2*Msolar;
% M = m1_val + m2_val; %Should be in solar mass
% u = m1_val*m2_val/M; %Should be in solar mass
% chirpmass = (u^3*M^2)^(1/5);
% pzi = 3*((chirpmass/Msolar)^(-5/3))*(fmin/100)^(-8/3);
% Nfac = (1/sqrt(50))*(fmin/100)^(2/3);
end

