% rungwpso test script

clearvars
close all
%Defines paths and creates folders for day's runs
addpath("References\")
anadpfc;

%JSON and filetag creation
psoParams = loadjson(jobParams.psoParamsjson);
signalParams = loadjson(jobParams.signalParamsjson);

% File Naming Conventions
filetagstr = filetagana(psoParams,signalParams);
%Creates a descriptive file name with important data parameters
outdataFilePrfx = [outDir,jobParams.jobName,'_',filetagstr];
%% mtchdfltr Data Load - pwelch PSD filtering
outData = load_mtchdfltrdata(jobParamsFile);

%% Glitch Checking
figure; 
plot(outData.dataY)
strtTime = 6182100;
endTime = strtTime+floor(60/outData.tIntrvl);
tseriestrainSeg = outData.dataY(strtTime:endTime);
[S,F,T] = spectrogram(tseriestrainSeg,8196,8000,[],outData.sampFreq);
S = abs(S);
imagesc(F,T,log10(S)) %Checking spectrogram image for glitches or high noise
title('Training Segment Spectrogram')
winVec = tukeywin(4*sampFreq);
[PSDtrainSeg,freqVectrainSeg] = pwelch(tseriestrainSeg,winVec,[],[],outData.sampFreq);
%% Zeroing beyond frequency cutoffs
hpcfreqInd = find(freqVectrainSeg<=outData.FreqBnd(1),1,"last"); 
PSDtrainSeg(1:hpcfreqInd) = 0;
lpcfreqInd = find(freqVectrainSeg>=outData.FreqBnd(2),1,"first");
PSDtrainSeg(lpcfreqInd:end) = 0;

figure;
plot(outData.freqVecog,outData.PSDog); axis tight
title('Welch"s PSD')
figure;
plot(outData.freqVec,outData.PSD); axis tight
title('Log10 PSD with Zeroed Frequencies')
figure;
semilogy(freqVectrainSeg,PSDtrainSeg); axis tight
% plot(freqVectrainSeg,PSDtrainSeg); axis tight
title('Training Segment PSD')

%% Run SHAPES on pwelch PSD
run_test_drase4lines(userUID);

%% Interpolation
PSD = createPSD(PSDtrainSeg,freqVectrainSeg,4096,4096,jobParams.outFilePSD);

%% Condition Data and Compute FFTs
[fftdataY,TF] = cond_mtchdfltrdata(jobParamsFile, ...
    signalParams.sampling_freq*signalParams.signal.T_sig_len);

%% Run PSO on GW data 
tic;
rungwpso('PCmtchdfltrtest_Job_params.json',outdataFilePrfx)
toc;
%% Post Processing
load([outdataFilePrfx,'C'])
load([outdataFilePrfx,'params'])

figure;
hold on;
plot(tauParams.dataX,tauParams.dataY,'.');
% plot(dataX,wave,'r');
% for lpruns = 1:nRuns
%       plot(dataX,outStruct.allRunsOutput(lpruns).estSig,'Color',[51,255,153]/255,'LineWidth',4.0);
% end
plot(tauParams.dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
% legend('Data','Signal',...
%         ['Estimated signal: ',num2str(nRuns),' runs'],...
%         'Estimated signal: Best run');
legend('Data','Signal',...
        'Estimated signal: Best run');
% saveas(gcf,files.psoresultplot);
hold off;

figure;
iterVec = linspace(1,psoParams.maxSteps,psoParams.maxSteps);
hold on;
for lpruns = 1:nRuns
      plot(iterVec,outStruct.allRunsOutput(lpruns).allBestFit, 'DisplayName',num2str(lpruns));
end
title("Best Fitness Values for All Runs");
xlabel("Iteration");
ylabel("Best Fitness Value");
legend;
% saveas(gcf,files.bestfitplot);
hold off;

if psoParams.type == "tau"
    figure;
    hold on;
    for lpruns = 1:nRuns
          rVec = s2rv(outStruct.allRunsOutput(lpruns).allBestLoc,tauParams);
          plot(rVec(:,1),rVec(:,2),'DisplayName',num2str(lpruns));
    end
    scatter(tau0,tau1p5,140,'red','filled','D','DisplayName','Original Parameters');
    title("Best Parameter Values for All Runs");
    xlabel("\tau_0");
    ylabel("\tau_{1.5}");
    legend;
    boundary_plot;
%     saveas(gcf,files.bestlocplot);
    hold off;

    t0 = outStruct.bestQcCoefs(1);
    t1p5 = outStruct.bestQcCoefs(2);
    est_M = (5/(32*tauParams.frange(1)))*(t1p5/(pi*pi*t0))*(c^3/G);
    est_u = (1/(16*tauParams.frange(1)^2))*(5/(4*pi^4*t0*t1p5^2))^(1/3)*(c^3/G);
    
    est_m1 = (est_M - sqrt(est_M^2 - 4*est_u*est_M))/2;
    est_m2 = (est_M + sqrt(est_M^2 - 4*est_u*est_M))/2;
    
    %% This will display parameters given through signal.json and PSO-estimated parameters
    %% Uncomment Original parameter display command if needed
%     disp(['Original parameters: tau0= ',num2str(tau0),...
%                                   '; tau1p5= ',num2str(tau1p5),...
%                                   '; m1= ', num2str(m1/Msolar),...
%                                   '; m2= ', num2str(m2/Msolar),...
%                                   '; A = ',num2str(signalParams.snr),...
%                                  '; phi = ',num2str(signalParams.phase),...
%                                   '; t_a = ',num2str(signalParams.ta),...
%                                   '; FitVal = ',num2str(original_fitVal)]);
    
    disp(['Estimated parameters: tau0=',num2str(outStruct.bestQcCoefs(1)),...
                                  '; tau1p5=',num2str(outStruct.bestQcCoefs(2)),...
                                  '; m1= ', num2str(est_m1/Msolar),...
                                  '; m2= ', num2str(est_m2/Msolar),...
                                  '; A = ',num2str(outStruct.bestAmp),...
                                 '; phi = ',num2str(outStruct.bestPhase),...
                                  '; t_a = ',num2str(outStruct.bestTime),...
                                  '; FitVal = ',num2str(bestFitVal)]);
else
    figure;
    hold on;
    for lpruns = 1:nRuns
          rVec = s2rv(outStruct.allRunsOutput(lpruns).allBestLoc,tauParams);
          plot(rVec(:,1),rVec(:,2),'DisplayName',num2str(lpruns));
    end
    scatter(m1,m2,140,'red','filled','D','DisplayName','Original Parameters');
    title("Best Parameter Values for All Runs");
    xlabel("m_1");
    ylabel("m_2");
    legend;
%     saveas(gcf,files.bestlocplot);
    hold off;
    
    %% This will display parameters given through signal.json and PSO-estimated parameters
    %% Uncomment Original parameter display command if needed
%     disp(['Original parameters:  m1= ',num2str(m1),...
%                                   '; m2= ',num2str(m2),...
%                                   '; A = ',num2str(signalParams.snr),...
%                                  '; phi = ',num2str(signalParams.phase),...
%                                   '; t_a = ',num2str(signalParams.ta),...
%                                   '; FitVal = ',num2str(original_fitVal)]);

    disp(['Estimated parameters: m1=',num2str(outStruct.bestQcCoefs(1)),...
                              '; m2=',num2str(outStruct.bestQcCoefs(2)),...
                              '; A = ',num2str(outStruct.bestAmp),...
                             '; phi = ',num2str(outStruct.bestPhase),...
                              '; t_a = ',num2str(outStruct.bestTime),...
                              '; FitVal = ',num2str(bestFitVal)]);
end


