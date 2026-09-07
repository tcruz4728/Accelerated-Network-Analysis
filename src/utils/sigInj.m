function [signal] = sigInj(params,PSDtotal)
% function [signal] = sigInj(fpos, ta, phase, fmin, fmax, m1,m2,r,datalen, initial_phase, N, avec, A)
% [S] = sigInj(P,T)
%SIGINJ Returns time domain vector of waveform for injection with custom
%strain amplitude
%   This function creates a time domain vector of a custom injected
%   waveform normalized to a strain amplitude which is a function of the
%   component masses and the distance to the source. 
% This factor is given in Cutler and Flanagan 1994 for the newtonian waveform 
% (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.49.2658)
%  Input: P, a parameter structure containing the following fields:
%         fpos: Positive frequency vector, 
%         ta: time of arrival, 
%         phase: coalescence phase, 
%         [fmin, fmax]: waveform frequency bounds, 
%         {m1,m2,r}: the component masses and the distance in Mpc respectively, 
%         datalen: data length in seconds, 
%         initial_phase: initial phase, 
%         N: total number of samples, 
%         avec: Precalculated vectors for phase generation, 
%         A: frequency magnitude vector.
% T, the complete two-sided design sensitivity PSD used for
%         computing the normalization factor
% Output: signal: time-domain vector of waveform with normalized strain
% amplitude

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
normfac = 1/sqrt((1/(N*signalParams.sampling_freq))*sum((fft(waveVec)./PSDtotal).*conj(fft(waveVec))));
% normfac = sqrt(N/sum((fft(waveVec)./PSDtotal).*conj(fft(waveVec))));
% normfac = 1/sqrt(innerproduct(waveVec,waveVec,PSDtotal));

% Create final signal
signal = params.signal.snr*normfac*waveVec;
% signal = sqrt(signalParams.sampling_freq)*waveVec;
end

