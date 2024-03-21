function [mfVal, max_arg,varargout] = mfgw(gwCoefs,params)

%MatchedFiltering for Chirp time space PSO
%Raghav Girgaonkar, April 2023
%Updated to include optional outputs from getparamestimates.m - TC, Jan 2024

%Generate normalized quadrature templates
% tau0 = x(1);
% tau1p5 = x(2);
fftq0 = gen2PNwaveform(params,0,0,gwCoefs,1);
disp(num2str(size(fftq0)))
disp(num2str(length(params.fftdataYbyPSD)))
% phaseq0 = gen2PNwaveform_tau(params.fpos, 0, 0, params.frange(1), params.frange(2), tau0,...
%     tau1p5,params.datalen,0,1,params.N,params.avec, params.normfac);

fftq1 = fftq0.*params.phaseDiff;

prod = params.T_sig*params.Fs;

%Compute fitness value after maximizing by matched filtering
mf1 = matchedfiltering(params.fftdataYbyPSD, fftq0);
mf2 = matchedfiltering(params.fftdataYbyPSD, fftq1);
[max_val, max_arg] = max(mf1(1:end - prod).^2 + mf2(1:end - prod).^2);
mfVal = -1*max_val;
varargout{1} = [mf1;mf2];
%Estimated SNR
varargout{2} = sqrt(max_val); %estAmp
%Estimated ToA:
varargout{3} = max_arg/params.Fs; %estTa

%Estimated Phase
yq0 = mf1(max_arg);
yq1 = mf2(max_arg);

varargout{4} = atan2(yq1,yq0); %estPhase