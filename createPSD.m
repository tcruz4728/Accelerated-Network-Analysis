%% Script to create interpolated PSD from SHAPES and PWELCH Estimates
function varargout = createPSD(inPSD,freqVec,T_sig_len,sampFreq,outFileName)
% O = createPSD(I,F,L)
%Takes in PSD I (a hand-picked segment of low noise), sampling frequency F
%and window length L to interpolate a new PSD. The output PSD will be used
%to compute the transfer function in cond_mtchdfltrdata.m
%Thomas Cruz, May 2023, derived from Raghav's
%% Data Parameters
nSamples = sampFreq*T_sig_len; %Total number of samples

kNyq = floor(nSamples/2);
fvec = (0:(kNyq))*sampFreq/nSamples;

%% 1-D Interpolation
% n = winLen*sampFreq; %Number of samples per window
% freqs = (0:n)*(1/(2*winLen));
logPSD = log10(inPSD);
loginterpPSD = interp1(freqVec, logPSD, fvec);

% %% Antilog
interpPSD = (10.^loginterpPSD);

varargout{1} = interpPSD;
save(outFileName,"interpPSD",'-append')
end