% Generate and Inject custom CBC signal
function outData = cbcsiginj(inData,psoParams,signalParams)
%This function reads input signal parameters from signal.json and
%creates a custom CBC signal in the time domain that can be injected into
%the data realization. 

if psoParams.type == "tau"
    %     wavephase = gen2PNwaveform_tau(fpos, ta, phase, fmin, fmax,tau0,tau1p5,datalen, initial_phase, snr, N, avec, genNormfac);
    wavephase = gen2PNwaveform_tau(signalParams,signalParams.ta,qcCoefs,signalParams.snr);
else
    wavephase = gen2PNwaveform(fpos, ta, phase, fmin, fmax, m1,m2,datalen, initial_phase, snr, N, avec, genNormfac);
end

wavefourier = psoParams.A.*wavephase;
%% Whitening the injected CBC signal to be consistent with strain data
%% Uncomment following line in the case of Conditioned LIGO HDF5 file, leave commented otherwise
wavefourier = wavefourier.*TFtotal;

%% Uncomment following lines in the case of SHAPES and WELCH PSD, leave commented otherwise
wavefourier = wavefourier.*psdVec4Norm;

wave = ifft(wavefourier);

%% Inject CBC signal into strain data
outData = inData + wave;
end