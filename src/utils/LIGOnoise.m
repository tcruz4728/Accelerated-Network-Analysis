function [outNoise, PSD,varargout] = LIGOnoise(N, Fs, noise_num, noisefile,varargin)
%Function to create colored noise using LIGO Design Sensitivities 
% Design PSD is modified between 15 Hz and 700Hz.
% Input: N = Total number of samples,
%        Fs = Sampling Frequency,
%        (Optional) noise_num = noise realization number from a pre-created noise realizations file
%        (Optional) noisefile = pre-created noise realizations filename
% Output: outNoise = colored noise vector,
%         PSD = two-sided PSD vector for positive DFT frequencies

% Raghav Girgaonkar, April 2023

%Optional Input arguments
freqBnds = [30 700];
if nargin > 4
    if ~isempty(varargin{1})
        freqBnds = varargin{1};
    end
end

%Load PSD 
y = load('iLIGOSensitivity.txt','-ascii');
% freqs = y(:,1);
% sqrtPSD = y(:,2);

%Turn one-sided sensitivity to two-sided
y(:,2) = (1/sqrt(2))*y(:,2);

% Interpolate sensitivity curve to positive DFT frequencies
minF = min(y(:,1));
maxF = max(y(:,1));
if minF ~= 0
% f=0 does not exist, put it in
y = [0, y(1,2);...
                  y];
end
if maxF < Fs/2
    error('High frequency limit requested is higher than supplied');
end


%Positive DFT frequencies
kNyq = floor(N/2)+1;
fvec = (0:(kNyq-1))*Fs/N;

%% Interpolation
interPSD = interp1(y(:,1),y(:,2), fvec);

%% Modifications, change cutoff frequencies as needed 
minidx = find(fvec<=freqBnds(1), 1, 'last' );
maxidx = find(fvec<=freqBnds(end), 1, 'last' );

SnBndStrt = interPSD(minidx);
SnBndEnd = interPSD(maxidx);
 
interPSD(1:minidx) = SnBndStrt;
interPSD(maxidx:end) = SnBndEnd;

PSD = interPSD.^2;
varargout{1} = fvec;
%% Make colored Noise
fltrOrdr = 20000;

outNoise_t = statgaussnoisegen(N,[fvec(:),PSD(:)],fltrOrdr,Fs, noise_num, noisefile);

outNoise = outNoise_t(fltrOrdr+1:end - fltrOrdr);