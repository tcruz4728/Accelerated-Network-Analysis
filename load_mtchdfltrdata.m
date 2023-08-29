function varargout = load_mtchdfltrdata(jobParamsFile,varargin)
% Function to load LIGO time series (read a hdf5 file) for use with matched
% filtering.
% LOAD_MTCHDFLTRDATA(F)
% Loads json file F which contains directory paths to LIGO .HDF5 time
% series and PSD .mat files for loading and saving files, respectively. 
% Applies a high-pass FFT FIR filter of order 30 with frequency cutoff 15 Hz 
% on the time series data. Then computes Power Spectral Density (PSD)
% of high-passed time series data using a Tukey window with window length
% 4*sampFreq. Takes log10 of PSD and zeroes its strain values with
% corresponding frequencies below cutoff 15 Hz.
% Lastly, saves the following to inFilepwPSD designated in F as .mat:
%   'PSD' -  Log base 10 of Power Spectral Density of high-passed
%   unwhitened time series data with zeroed power for frequencies below
%   high-pass cutoff
%   'freqVec' - Frequency Vector to PSD in Hz
%   'dataY' - High-passed unwhitened time series data
%   'sampFreq' - Sampling Frequency (Hz)
% Optional Arguments
% LOAD_MTCHDFLTRDATA(F,C,W,O)
% C specifies the starting cutoff frequency in Hz, W controls the window length
% factor, and O is the FFT FIR filter order.

%Thomas Cruz, May 2023

%% Default Parameters
jobParams = loadjson(jobParamsFile);
strtFreq = 30; % Low Frequency cutoff, or starting frequency
endFreq = 900; % High Frequency cutoff, or ending frequency
% winLenFac = 4; %Window length factor, winLen = winLenFac*sampFreq;
% filtordr = 30; %Filter order for highpass filter

%Override default parameters if given
nreqArgs = 1;
for lpargs = 1:(nargin-nreqArgs)
    if ~isempty(varargin{lpargs})
        switch lpargs
            case 1
                strtFreq = varargin{lpargs};
            case 2
                endFreq = varargin{lpargs};
            case 3
                winLenFac = varargin{lpargs};
            case 4
                filtordr = varargin{lpargs};
        end
    end
end

%% Load LIGO HDF5 File strain data
%Data is assumed to be unwhitened time series strain data from a GW
%detector in .hdf5 file format. extracting the strain data, sampling
%frequency and time interval between data points.
dataY = h5read(jobParams.inFileData,'/strain/Strain')';
% fileInfo = h5info(jobParams.inFileData);
% nSamples = fileInfo.Groups(3).Datasets.Attributes(1).Value;
nSamples = double(h5readatt(jobParams.inFileData,'/strain/Strain','Npoints'));
tIntrvl = double(h5readatt(jobParams.inFileData,'/strain/Strain','Xspacing')); %Time Interval
tlen = nSamples*tIntrvl;
%%%%%%% TEMP RUN WITH NOISE SAMPLE
% load(jobParams.inFileData,'strain');
% dataY = strain;
% nSamples = length(dataY);
% tIntrvl = 1/4096;
% tlen = nSamples*tIntrvl;
%%%%%%%
sampFreq = 1/tIntrvl;

sampFreq = double(sampFreq); %Converts to double from int64
if sampFreq ==16384 %resamples if data has the 16kHz sampling freq
    sampFreq = 4096;
    dataY = resample(dataY,1,4);
end

%% Band-pass filter on unwhitened time series data
dataYwin = tukeywin(size(dataY,2),4*sampFreq/(size(dataY,2)));
dataY = dataY.*dataYwin';

dataY_highpass = highpass(dataY, strtFreq, sampFreq, ImpulseResponse="iir",Steepness=0.95);

dataY = dataY_highpass;

%% Time series training segment generation
strtTime = floor(10*sampFreq);
if  tlen < 60
    endTime = strtTime + floor(32/tIntrvl);
else
    endTime = strtTime + floor(60/tIntrvl);
end
tseriestrainSeg = dataY_highpass(strtTime:endTime);

winVec = tukeywin(4*sampFreq);
[PSD,freqVec] = pwelch(tseriestrainSeg,winVec,[],[],sampFreq);
PSD = PSD.';
freqVec = freqVec.';

freqBnd = [strtFreq,endFreq]; 
outData = struct('PSD',PSD,'freqVec',freqVec, ...
    'dataY',dataY_highpass,'sampFreq', sampFreq, ...
    'FreqBnd',freqBnd,'tlen',tlen,...'PSDog',PSDog,'freqVecog',freqVecog,...
    'tseriestrainSeg',tseriestrainSeg);
varargout{1} = outData;

save(jobParams.inFilePSD,"PSD","freqVec","dataY","sampFreq",'freqBnd','tlen')

