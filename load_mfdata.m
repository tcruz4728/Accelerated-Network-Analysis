function varargout = load_mfdata(inFileData,inFilePSD,varargin)
% Function to load LIGO time series for use with matched filtering.
% LOAD_MTCHDFLTRDATA(D,P)
% Loads data file D which can either be a .hdf5 or a .mat file. Applies a
% Tukey Window (see tukeywin.m) which windows the edges of the time series
% and then applies a high-pass FFT IIR filter with frequency cutoff 30 Hz
% on the time series data. Then computes the Power Spectral Density (PSD)
% using the first 32 or 60 seconds of the high-passed time series data
% using a Tukey window with window length 4*sampFreq.
% Saves the following varaibles to P designated in D as a .mat file:
%   'PSD'-      Power Spectral Density of high-passed unwhitened time-
%               series data using pwelch.m.
%   'freqVec'-  Frequency Vector to PSD in Hz.
%   'dataY'-    High-passed unwhitened time-series data, row vector
%   'sampFreq'- Sampling Frequency (Hz).
%   'freqBnd'-  Frequency bounds of interest, used later to filter.
%   'tlen'-     Length of the time series in seconds.
%
% Optional Arguments
% O = LOAD_MTCHDFLTRDATA(D,P,S,E)
% S specifies the starting cutoff frequency in Hz, E is the ending
% frequency cutoff.
% O is the optional output argument which returns the data saved in P as a
% structure with fields of the same name. In addition to those variables,
% the time series segment used to compute the pwelch PSD is returned.
%
% Thomas Cruz, May 2023
%
% See also pwelch, tukeywin, highpass

%% Default Parameters
strtFreq = 30; % Low Frequency cutoff, or starting frequency
endFreq = 700; % High Frequency cutoff, or ending frequency

%Override default parameters if given
nreqArgs = 2;
sigInjChk = [];
window_time = 4;
for lpargs = 1:(nargin-nreqArgs)
    if ~isempty(varargin{lpargs})
        switch lpargs
            case 1
                sigInjChk = varargin{lpargs};
            case 2
                strtFreq = varargin{lpargs};
            case 3
                endFreq = varargin{lpargs};
        end
    end
end

%% Load LIGO .hdf5 File strain data or generated .mat data
if isstruct(inFileData)
    dataY = inFileData.dataY;
    tIntrvl = inFileData.tIntrvl;
    nSamples = length(dataY);
else
    [~,~,fileExt] = fileparts(inFileData);
    switch lower(fileExt)
        case '.hdf5' % Data is assumed to be unwhitened time series strain data
            % from a GW detector.
            % Extracting the strain data, sampling frequency and time interval between
            % data points.
            dataY = h5read(inFileData,'/strain/Strain')';
            nSamples = double(h5readatt(inFileData,'/strain/Strain','Npoints'));
            tIntrvl = double(h5readatt(inFileData,'/strain/Strain','Xspacing')); %Time Interval
        case '.mat'
            load(inFileData,'dataY','tIntrvl')
            if ~isrow(dataY)
                dataY = dataY';
            end
            disp(['load_mfdata- ',inFileData])
            nSamples = length(dataY);
    end
end

tlen = nSamples*tIntrvl;
sampFreq = 1/tIntrvl;

%Downsampling to 4 kHz
sampFreq = double(sampFreq); %Converts to double from int64
if sampFreq == 16384 
    sampFreq = 4096;
    dataY = resample(dataY,1,4);
end
%% Signal Injection
if ~isempty(sigInjChk)
    load(inFileData,'injSigparams')
    dataY = dataY + injSigparams.signal.data;
    disp(['load_mfdata- Injected signal with snr: ', num2str(injSigparams.signal.snr)])
end
%% Band-pass filter on unwhitened time series data
dataYwin = tukeywin(size(dataY,2),window_time*sampFreq/(size(dataY,2)));
dataY = dataY.*dataYwin';
dataY_highpass = highpass(dataY, strtFreq, sampFreq, ImpulseResponse="iir",Steepness=0.95);
dataY = dataY_highpass;

%% Time series training segment generation
strtTime = floor(10*sampFreq);
if  tlen <= 60 && tlen > 32
    endTime = strtTime + floor(32/tIntrvl);
elseif tlen <= 32
    endTime = strtTime + floor(tlen/(4*tIntrvl));
else
    endTime = strtTime + floor(60/tIntrvl);
end
tseriestrainSeg = dataY_highpass(strtTime:endTime);

%% Pwelch PSD estimation
winVec = tukeywin(window_time*sampFreq);
[PSD,freqVec] = pwelch(tseriestrainSeg,winVec,[],[],sampFreq);

PSD = PSD/2; % /2 to convert to 2-sided psd

%outputs need to be column vectors
if ~iscolumn(PSD)
    PSD = PSD.';
end

if ~iscolumn(freqVec)
    freqVec = freqVec.';
end

% plot(freqVec,log10(PSD),'b')
freqBnd = [strtFreq,endFreq]; 

%% Outputs
% if genSig == 1
%     outData = struct('PSD',PSD,'freqVec',freqVec, ...
%         'dataY',dataY,'sampFreq', sampFreq, ...
%         'freqBnd',freqBnd,'tlen',tlen,...
%         'tseriestrainSeg',tseriestrainSeg,...
%         'dsstPSD',dsstPSD,'dsstfreqVec',dsstfreqVec,...
%         'injectedSignal',injSigparams.signal.data,'genSig',genSig);
% else

outData = struct('PSD',PSD,'freqVec',freqVec, ...
    'dataY',dataY,'sampFreq', sampFreq, ...
    'freqBnd',freqBnd,'tlen',tlen,...
    'tseriestrainSeg',tseriestrainSeg);
if ~isempty(sigInjChk)
    outData.injectedSignal = injSigparams.signal.data;
end

varargout{1} = outData;

if ~isempty(inFilePSD)
    whos('PSD','freqVec','dataY','sampFreq') % debugging
    save(inFilePSD,'PSD','freqVec','dataY','sampFreq',...
        'freqBnd','tlen', "-v7")
    if ~isempty(sigInjChk)
        save(inFilePSD,'injSigparams',"-append")
    end
    disp(['load_mfdata- pwelch PSD and time series data saved to: ',inFilePSD])
end


