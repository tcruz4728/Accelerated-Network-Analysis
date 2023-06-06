% function outResults = crcbqcpso_tau(inParams,psoParams,nRuns, sampling_freq)
function outResults = crcbqcpso_tau(inParams,psoParams,nRuns)
%Regression of 2PNWaveform using Chirp-Time Space PSO
% inParams: Struct containing data and signal parameters
% psoParams: Struct containing PSO parameters
% nRuns: Number of PSO iterations
% sampling_freq: Sampling Frequency of data

%The fields of O are:
% 'allRunsOutput': An N element struct array containing results from each PSO
%              run. The fields of this struct are:
%                 'fitVal': The fitness value.
%                 'qcCoefs': The coefficients [tau0, tau1.5].
%                 'estSig': The estimated signal.
%                 'totalFuncEvals': The total number of fitness
%                                   evaluations.
% 'bestRun': The best run.
% 'bestFitness': best fitness from the best run.
% 'bestSig' : The signal estimated in the best run.
% 'bestQcCoefs' : [tau0, tau1.5] found in the best run.

%Raghav Girgaonkar, April 2023

nSamples = length(inParams.dataX);

fHandle = @(x) psofitfunc_tau(x,inParams);

params = inParams;

nDim = 2;
outStruct = struct('bestLocation',[],...
                   'bestFitness', [],...
                   'totalFuncEvals',[],...
                   'allBestFit',[],...
                   'allBestLoc',[]);
                    
outResults = struct('allRunsOutput',struct('fitVal', [],...
                                           'qcCoefs',zeros(1,2),...
                                           'estTa',[],...
                                           'estSig',zeros(1,nSamples),...
                                           'totalFuncEvals',[],...
                                           'allBestFit',zeros(1,psoParams.maxSteps),...
                                           'allBestLoc',zeros(nDim,psoParams.maxSteps)),...
                    'bestRun',[],...
                    'bestFitness',[],...
                    'bestSig', zeros(1,nSamples),...
                    'bestQcCoefs',zeros(1,2),...
                    'estAmp',[],...
                    'estPhase',[],...
                    'bestAmp',[],...
                    'bestPhase',[],...
                    'bestTime',[]);

%Allocate storage for outputs: results from all runs are stored
for lpruns = 1:nRuns
    outStruct(lpruns) = outStruct(1);
    outResults.allRunsOutput(lpruns) = outResults.allRunsOutput(1);
end
%Independent runs of PSO in parallel. Change 'parfor' to 'for' if the
%parallel computing toolbox is not available.
% fprintf("Running PSO\n");
% parpool(nRuns);
for lpruns = 1:nRuns
    disp(['Run ',num2str(lpruns),' of ', num2str(nRuns)])
    %Reset random number generator for each worker
    rng(lpruns);
    tic;
    outStruct(lpruns)=crcbpso(fHandle,nDim,psoParams,2);
    toc;
end

%Prepare output
fitVal = zeros(1,nRuns);
% sampling_interval = 1/sampling_freq;
%Time vectors for signal and total series
% timeVecSig = (0:(sampling_freq*T_sig - 1))/sampling_freq;
for lpruns = 1:nRuns
    tic;
    outResults.allRunsOutput(lpruns).allBestFit = outStruct(lpruns).allBestFit;
    outResults.allRunsOutput(lpruns).allBestLoc = outStruct(lpruns).allBestLoc;
    fitVal(lpruns) = outStruct(lpruns).bestFitness;
    outResults.allRunsOutput(lpruns).fitVal = fitVal(lpruns);

    [~,qcCoefs,ta_index] = fHandle(outStruct(lpruns).bestLocation);

    outResults.allRunsOutput(lpruns).qcCoefs = qcCoefs;
    outResults.allRunsOutput(lpruns).ta_index = ta_index;
    
    %Calculate time using sampling freq and ta_index
    estTa = ta_index/inParams.Fs;
    outResults.allRunsOutput(lpruns).estTa = estTa;

    phaseq0 = gen2PNwaveform_tau(params,ta_index,0,qcCoefs,1);
    fftq0 = phaseq0;
    fftq1 = phaseq0.*params.phaseDiff;
    %Estimated Phase
    yq0 = innerprodpsd(fftq0, params.fftdataYbyPSD);
    yq1 = innerprodpsd(fftq1, params.fftdataYbyPSD);
    estPhase = atan2(yq1,yq0);
    outResults.allRunsOutput(lpruns).estPhase = estPhase;

    %Estimated Amplitude
    estAmp = cos(estPhase)*yq0 + sin(estPhase)*yq1;
    outResults.allRunsOutput(lpruns).estAmp = estAmp;

    %Estimated Signal
    estSigphase = gen2PNwaveform_tau(params,ta_index,estPhase,qcCoefs,estAmp);
    estSigfourier = (params.A).*estSigphase;
    estSig = ifft(estSigfourier);

    outResults.allRunsOutput(lpruns).estSig = estSig;
    outResults.allRunsOutput(lpruns).totalFuncEvals = outStruct(lpruns).totalFuncEvals;
    toc;
end
%Find the best run
[~,bestRun] = min(fitVal(:));
outResults.bestRun = bestRun;
outResults.bestFitness = outResults.allRunsOutput(bestRun).fitVal;
outResults.bestSig = outResults.allRunsOutput(bestRun).estSig;
outResults.bestAmp = outResults.allRunsOutput(bestRun).estAmp;
outResults.bestPhase = outResults.allRunsOutput(bestRun).estPhase;
outResults.bestQcCoefs = outResults.allRunsOutput(bestRun).qcCoefs;
outResults.bestTaIndex = outResults.allRunsOutput(bestRun).ta_index;
outResults.bestTime = outResults.allRunsOutput(bestRun).estTa;

