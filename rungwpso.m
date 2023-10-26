function []=rungwpso(paramsFile,outFileName)
%% Script to run Chirp-time or Mass Space PSO on a user specified datafile and PSD. 

%Raghav Girgaonkar, Apr 2023
%Thomas Cruz, May 2023, updated for use with ALINE

load(paramsFile,'gwCoefs','params','psoParams'); % Loads in params structure

original_fitVal = -1*mfqc([gwCoefs, 2], params);
if isempty(gcp('nocreate'))
    parpool('Processes',psoParams.nRuns)
end
outStruct = crcbgwpso(params,psoParams,2);
bestFitVal = -1*outStruct.bestFitness;
save(outFileName,'outStruct','bestFitVal',"original_fitVal",'params');
end
