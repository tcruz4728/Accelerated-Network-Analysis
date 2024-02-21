function varargout = createPSD(varargin)
% File input
% O = createPSD(N)
%Loads in data file N with 4 variables; PSD, freqVec, tlen, and sampFreq as
%described under value inputs.
% The output PSD will be used to compute the transfer function in cond_mtchdfltrdata.m
%Value Inputs
% O = createPSD(P,F,L,S)
%Takes in PSD P (a hand-picked segment of low noise),frequency vector F,
%length of signal in seconds, and the sampling frequency S, to interpolate
%the PSD to match the frequency vector's size created from the signal length and
%sampling frequency.
%Thomas Cruz, May 2023, derived from Raghav's

%%Optional Input arguments
if nargin<2 %input file
    outFileName = varargin{1};
    load(outFileName,'PSD','freqVec','tlen','sampFreq')
else %input values
    for largs = 1:5
        switch largs
            case 1
                PSD = varargin{largs};
            case 2
                freqVec= varargin{largs};
            case 3
                tlen = varargin{largs};
            case 4
                sampFreq = varargin{largs};
            case 5
                outFileName = varargin{largs};
        end
    end

end

%% Data Parameters
nSamples = sampFreq*tlen; %Total number of samples

kNyq = floor(nSamples/2);
fvec = (0:(kNyq))*sampFreq/nSamples;

%% 1-D Interpolation
% n = winLen*sampFreq; %Number of samples per window
% freqs = (0:n)*(1/(2*winLen));
logPSD = log10(PSD);
loginterpPSD = interp1(freqVec, logPSD, fvec);

% %% Antilog
interpPSD = (10.^loginterpPSD);

varargout{1} = interpPSD;
save(outFileName,"interpPSD",'-append')
end