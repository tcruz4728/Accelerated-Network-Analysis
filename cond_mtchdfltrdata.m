function [varargout] = cond_mtchdfltrdata(inputFile,paramsFileName,varargin)
% Function that handles data reading and conditioning before being called
% by rungwpso.m. Works by reading a flag specified in the jobParamsFile
% which specifies the kind of data being used, either unified under a
% single hdf5 (containing whitened time series data and transfer function)
% or separated under a  .hdf5 file and .mat file (containing unwhitened
% time series data and a PSD from pwelch or shapes). It loads the files and
% computes the FFT on both time series. It has outputs of the
% FFT(timeseries) and Transfer function (1/sqrt(PSD)).

% jobParams = loadjson(jobParamsFile);
load(paramsFileName,"params")

%% Create entire PSD vector

negFStrt = 1-mod(params.N,2);
kNyq = floor(params.N/2)+1;
% Load data
% switch datatype
%     case 1 %'welch'
load(inputFile,"interpPSD","sampFreq","dataY","freqBnd",'tlen'); %Welchs PSD
% %     case 2 %'shps'
% %         load(jobParams.inFileshpsPSD,'interpPSD','sampFreq',"dataY","freqBnd",'tlen') %SHAPES estimated Welchs PSD
% % end
% load(jobParams.inFilePSD,"dataY","freqBnd",'tlen'); %highpassed dataY

%Perform signal injection (optional)
if ~isempty(varargin{1})
    signal = sigInj(params);
    dataY = dataY + 1*signal;
end

%Windowing before fft
dataYwin = tukeywin(size(dataY,2),0.5*sampFreq/(size(dataY,2)));
dataY = dataY.*dataYwin';

%Compute transfer function and FFT
TF = 1./sqrt((interpPSD./2));
TF(1:(freqBnd(1)*tlen)) = 0; %Zeroing frequencies below min frequency
TFtotal = [TF, TF((kNyq-negFStrt):-1:2)]; %Computed from training segment
fftdataYbyPSD = fft(dataY).*TFtotal.*TFtotal.*params.A;

%% Create General Normalization Factor
AbysqrtPSD = params.A.*TFtotal;
innProd = (1/params.N)*(AbysqrtPSD)*AbysqrtPSD';
genNormfac = 1/sqrt(real(innProd));

%% Saving Files
params.dataY = dataY;
params.fftdataYbyPSD = fftdataYbyPSD;
params.normfac = genNormfac;
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
save(paramsFileName,'params','-append')
end