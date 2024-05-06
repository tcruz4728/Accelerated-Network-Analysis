function [varargout] = cond_mfdata(inData,paramsData)
% Function that handles data reading and conditioning before being called
% by rungwpso.m. Works by reading a flag specified in the jobParamsFile
% which specifies the kind of data being used, either unified under a
% single hdf5 (containing whitened time series data and transfer function)
% or separated under a  .hdf5 file and .mat file (containing unwhitened
% time series data and a PSD from pwelch or shapes). It loads the files and
% computes the FFT on both time series. It has outputs of the
% FFT(timeseries) and Transfer function (1/sqrt(PSD)).

% jobParams = loadjson(jobParamsFile);
if isstruct(paramsData)
    params = paramsData;
else
    load(paramsData,"params")
end
%% Create entire PSD vector

% negFStrt = 1-mod(params.N,2);
% kNyq = floor(params.N/2)+1;

% Load data either via structure or .mat file
if isstruct(inData)
    interpPSD = inData.interpPSD;
    sampFreq = inData.sampFreq;
    dataY = inData.dataY;
    freqBnd = inData.freqBnd;
    % dsstPSD = inData.dsstPSD;
else %.mat load
    % load(inData,"interpPSD","sampFreq","dataY",'freqBnd','dsstPSD'); %PSD
    load(inData,"interpPSD","sampFreq","dataY",'freqBnd');
end

%Windowing before fft
dataYwin = tukeywin(size(dataY,2),0.5*sampFreq/(size(dataY,2)));
dataY = dataY.*dataYwin';

%Whiten data and create transfer function
[whtndfiltdata, TFtotal]=segdatacond(dataY, interpPSD, sampFreq,...
    [1,32*sampFreq],freqBnd(1,1));

% %Perform signal injection (optional)
% if ~isempty(varargin{1})
%     %Compute two-sided PSD from design sensitivity PSD for signal injection
%     dsstPSDtotal = [dsstPSD, dsstPSD((kNyq-negFStrt):-1:2)];
%     signal = sigInj(params,dsstPSDtotal);
%     %(Alternative)
%     % signal = gen2PNtemplate_mass(params,params.signal.ta,0,...
%     %     [params.gwCoefs,1],params.signal.snr,dsstPSDtotal);
% 
%     %q0 & q1 are phases for testing if the signal can be detected by mf
%     % q0 = gen2PNtemplate_mass(params,0,0,[params.gwCoefs,1],1,dsstPSDtotal);
%     % q1 = gen2PNtemplate_mass(params,0,pi/2,[params.gwCoefs,1],1,dsstPSDtotal);
% 
%     %Whiten the signal using the transfer function from PSD (pwelch or
%     %shps)
%     % whtndsignal = (1/sqrt(params.signal.sampling_freq))*ifft(fft(signal).*(TFtotal));
%     whtndsignal = ifft(fft(signal).*(TFtotal));
%     whtndfiltdata = whtndfiltdata + whtndsignal;
%     disp(['cond_mfdata- injected signal with snr: ', num2str(params.signal.snr)])
% end
fftdataYbyPSD = fft(whtndfiltdata).*TFtotal.*params.A;

% Optional plotting to check for proper signal injection
% mf1 = matchedfiltering(fftdataYbyPSD,q0);
% mf2 = matchedfiltering(fftdataYbyPSD,q1);
% mfts = sqrt(mf1.^2 + mf2.^2);
% figure; plot(mfts)
%% Create General Normalization Factor
AbysqrtPSD = params.A.*TFtotal;
innProd = (1/params.N)*(AbysqrtPSD)*AbysqrtPSD';
params.normfac = 1/sqrt(real(innProd));

%% Saving Files
params.dataY = dataY;
params.fftdataYbyPSD = fftdataYbyPSD;
% params.normfac = genNormfac;
if nargout > 0
    varargout{1} = params;
    if nargout > 1
        varargout{2} = fftdataYbyPSD;
        if naragout > 2
            varargout{3} = TF;
        end
    end
end
% save(jobParams.outFileData,"fftdataYbyPSD",'TFtotal','dataY');
if ~isempty(paramsData) && ~isstruct(paramsData)
    save(paramsData,'params','-append')
    disp(['cond_mfdata- Parameters MaJ with conditioned data and saved to: ',paramsData])
end
end