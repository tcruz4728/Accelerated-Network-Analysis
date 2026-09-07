function outParams = phaseCalc(signalParams)
%Function to compute frequency magnitudes and phase difference vectors
%Thomas Cruz, May 2023 extracted from rungwpso.m
%% Parameter Load
T_sig_len = signalParams.signal.T_sig_len; %Length of signal in seconds
Fs = signalParams.sampling_freq; %Sampling Frequency
nSamples = floor(T_sig_len*Fs); %Number of samples
fmin = signalParams.freq(1);
fmax = signalParams.freq(2);
fpos = (0:floor(nSamples/2))*(1/T_sig_len); %Positive frequency vector
%% Pre-processing
%Make Frequency Magnitude and Phase difference vectors
Apos = zeros(size(fpos));
phaseDiffpos = -1j*ones(size(fpos));

Apos(2:end) = fpos(2:end).^(-7/6);
%Modify Apos
min_index  = floor(T_sig_len*fmin) + 1;
max_index = floor(T_sig_len*fmax) + 1;
Apos(1:min_index-1) = 0;
Apos(max_index + 1: end) = 0;
%Make Aneg and Phasediffneg
if mod(nSamples,2) == 0
    Aneg = conj(Apos(end-1:-1:2));
    phaseDiffneg = conj(phaseDiffpos(end-1:-1:2));
else
    Aneg = conj(Apos(end:-1:2));
    phaseDiffneg = conj(phaseDiffpos(end:-1:2));
end
%Make full A
A = [Apos, Aneg];
phaseDiff = [phaseDiffpos, phaseDiffneg];

a0fvec = ((fpos(2:end)./fmin).^(-5/3));
a1fvec = ((fpos(2:end)./fmin).^(-4/3));
a2fvec = ((fpos(2:end)./fmin).^(-3/3));
a3fvec = ((fpos(2:end)./fmin).^(-2/3));
a4fvec = ((fpos(2:end)./fmin).^(-1/3));

avec = [a0fvec ; a1fvec; a2fvec; a3fvec; a4fvec];

outParams = struct('A',A,...
                  'phaseDiff', phaseDiff,...
                  'avec', avec);
end