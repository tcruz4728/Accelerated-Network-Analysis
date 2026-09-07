function []=rungwpso(paramsFile,outFileName)
%RUNGWPSO(P,O)
% This function reads a parameter file (.mat) and attempts a matched
% filtering of a signal in a PSD file within by applying Particle Swarm
% Optimization (PSO) efficiently with parallel computing and then saves the
% estimates of the signal, amplitude, and phase, while detailing the
% signal's time of arrival and SNR, all within a file named O. 
%
% Arguments
% P- a .mat file which contains the following:
% 'gwCoefs'- 1x2 array with gravitational wave coefficients (either chirp
% times or masses). Calculated by gwpsoparams.m using initial parameters,
% used by mfgw.m to compute the fitness value for comparison to matched
% filtering estimate.
% 'params'- parameter structure which will be used to compute the new
% estimates. See gwpsoparams for more information on fields.
% 'psoParams'- paramaters for the use of Particle Swarm Optimization (PSO)
% indicating the number of PSO iterations, maxSteps, and the number of
% computing cores to be called, nRuns.
%
% O- a character string which sets the name of the output file that will be
% saved containing the following variables.
% 'outStruct'- output structure from crcbgwpso.m containing 6 fields:
% 'allRunsOutput'-  Nested structure with each of the 'nRuns' estimates
%                   saved.
%     'fitVal'-     Fitness value which is minimized for the best fit   
%     'gwCoefs'-    Estimated gravitational wave coefficients
%     'estTa'-      Estimated time of arrival for detected signal
%     'estSig'-     Estimated signal isolated from matched filtering
%                   analysis
%     'totalFuncEvals'- Total number of fitness evaluations done
%     'allBestFit'- Structure of all the best fitness values 
%     'allBestLoc'- Structure of all the best location values 
%     'bestRun'-    Best run's index number
%     'bestFitness'-Best run's fitness value estimate
%     'bestSig'-    Best run's estimated signal 
%     'bestQcCoefs'-Best run's GW coefficients 
%     'estAmp'-     Estimated signal's amplitude
%     'estPhase'-   Estimated signal's phase
%     'bestAmp'-    Best run's estimated amplitude
%     'bestPhase'-  Best run's estimated phase
%     'bestTime'-   Best run's estimated time of arrival
%
%Raghav Girgaonkar, Apr 2023
%Thomas Cruz, May 2023, updated and streamlined for use with ALINE
%
% See also crcbgwpso, mfgw, crcbpso, gwpsoparams

load(paramsFile,'gwCoefs','params','psoParams'); % Loads in params structure

[original_fitVal,~,mf] = mfgw([gwCoefs, psoParams.type], params);
original_fitVal = -1*original_fitVal;
tic;
if isempty(gcp('nocreate'))
    parpool('Processes',psoParams.nRuns);
end
toc;
%Parameter that controls if the all fitness values and locations are
%returned, 1 if values, 2 if values and location, empty for none. 
allCtrl = 2;
outStruct = crcbgwpso(params,psoParams,allCtrl);
bestFitVal = -1*outStruct.bestFitness;
save(outFileName,'outStruct','bestFitVal',"original_fitVal",'params','mf');
disp(['rungwpso- Saved matched filtering outputs to ', outFileName]);
end
