% rungwpso test script
%Set UID for particular run or to overwrite a previous run
% if isempty(userUID)
    userUID = input("Type UID value: ");
% end
%Job File 
jobParamsFile = 'C:\Users\Thomas Cruz\Documents\GitHub\Accelerated-Network-Analysis\JSON\PCmtchdfltrtest_Job_params.json';
jobParams = loadjson(jobParamsFile);
%Defines paths and creates folders for day's runs
addpath(jobParams.path2drase)
setpath(jobParamsFile)
[outDir,filepaths] = dpfc(jobParams,userUID);
paramsFile = [outDir,'params'];
paramsFileshps = [paramsFile,'shps'];

fidprog = fopen([outDir,'progress.txt'],'a');
fprintf(fidprog,'%s\n',datestr(datetime('now')));

%JSON and filetag creation
psoParams = loadjson(jobParams.psoParamsjson);
signalParams = loadjson(jobParams.signalParamsjson);

% File Naming Conventions
filetagstr = filetagana(psoParams,signalParams);
%Creates a descriptive file name with important data parameters
outdataFilePrfx = [outDir,jobParams.jobName,'_',filetagstr];

%% Sets up params.m file containing matched filtering parameters that do 
% NOT need dataY and initial tau0 and tau1p5 values
params = gwpsoparams(jobParamsFile,paramsFile);
copyfile([paramsFile,'.mat'],[paramsFileshps,'.mat']);
%% mtchdfltr Data Load - pwelch PSD filtering
outData = load_mtchdfltrdata(jobParamsFile);
fprintf(fidprog,'\nLoad Data...done');
sampFreq = outData.sampFreq;
%% Glitch Checking
figure; 
plot(outData.dataY)
[S,F,T] = spectrogram(outData.tseriestrainSeg,8192,8000,[],sampFreq);
S = abs(S);
imagesc(T,F,log10(S)); axis xy; %Checking spectrogram image for glitches or high noise
title('Training Segment Spectrogram')
saveas(gcf,[filepaths.figs,'Training_Spectrogram']);
fprintf(fidprog,'\nGlitch Check...done');
%% Plots
figure;
semilogy(outData.freqVec,outData.PSD); axis tight
title('PSD of Training Segment')
saveas(gcf,[filepaths.figs,'Training_PSD']);

%% Run SHAPES on pwelch PSD
run_test_drase4lines(jobParams,userUID);
fprintf(fidprog,'\nDrase Run...done');

%% Interpolation
createPSD(outData.PSD,outData.freqVec,outData.tlen,sampFreq,jobParams.outFilePSD);
load(jobParams.inFileshpsPSD,"PSD")
createPSD(PSD,outData.freqVec,outData.tlen,sampFreq,jobParams.outFileshpsPSD);

%% Condition Data and Compute FFTs
paramsPW = cond_mtchdfltrdata(jobParamsFile,params,paramsFile,...
    sampFreq*outData.tlen,1);
paramsSHPS = cond_mtchdfltrdata(jobParamsFile,params, ...
    paramsFileshps,sampFreq*outData.tlen,2);
fprintf(fidprog,'Data Conditioning on Welch...done\n');

%% Run PSO on GW data 
rungwpso(paramsFile,outdataFilePrfx) %pwelch 
fprintf(fidprog,'Data Conditioning on SHAPES...done\n');
rungwpso(paramsFileshps,[outdataFilePrfx,'shps_']) %shapes estimate
fprintf(fidprog,'\nGW PSO Runs...done');

%% Post Processing
load([outdataFilePrfx,'C'])
% load(paramsFile,"gwCoefs")
gwCoefs = params.gwCoefs;
psoParams = params.pso;
signalParams = params.signal;

figure;
hold on;
plot(params.dataX,params.dataY,'.');

plot(params.dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
legend('Data','Signal',...
        'Estimated signal: Best run');
% saveas(gcf,[filepaths.figs,'PSO_Results']);
hold off;

figure;
iterVec = linspace(1,psoParams.maxSteps,psoParams.maxSteps);
hold on;
for lpruns = 1:psoParams.nRuns
      plot(iterVec,outStruct.allRunsOutput(lpruns).allBestFit, 'DisplayName',num2str(lpruns));
end
title("Best Fitness Values for All Runs");
xlabel("Iteration");
ylabel("Best Fitness Value");
legend;
% saveas(gcf,[filepaths.figs,'Best_fit']);
hold off;

figure;
hold on;
scatter(gwCoefs(1),gwCoefs(2),140,'red','filled','D','DisplayName','Original Parameters');
title("Best Parameter Values for All Runs");

if psoParams.type == 2
    xlabel("\tau_0");
    ylabel("\tau_{1.5}");
    legend;
    boundary_plot;
    hold off;

    t0 = outStruct.bestGwCoefs(1);
    t1p5 = outStruct.bestGwCoefs(2);
    est_M = (5/(32*params.frange(1)))*(t1p5/(pi*pi*t0))*(params.cgFac);
    est_u = (1/(16*params.frange(1)^2))*(5/(4*pi^4*t0*t1p5^2))^(1/3)*(params.cgFac);
    
    est_m1 = (est_M - sqrt(est_M^2 - 4*est_u*est_M))/2;
    est_m2 = (est_M + sqrt(est_M^2 - 4*est_u*est_M))/2;
    if jobParams.injSig == 1
        disp(['Original GW Coefficients: tau0= ',num2str(gwCoefs(1)),...
            '; tau1p5= ',num2str(gwCoefs(2)),...
            '; m1= ', num2str(signalParams.masses(1)),...
            '; m2= ', num2str(signalParams.masses(2))]);

    end
Msolar = 1.989*10^30; %Solar mass in kg
    % This will display parameters given through signal.json and PSO-estimated parameters
    disp(['Estimated GW Coefficients: tau0=',num2str(outStruct.bestGwCoefs(1)),...
                                  '; tau1p5=',num2str(outStruct.bestGwCoefs(2)),...
                                  '; m1= ', num2str(est_m1/Msolar),...
                                  '; m2= ', num2str(est_m2/Msolar)]);
else
    xlabel("m_1");
    ylabel("m_2");
    legend;
    hold off;  
    % This will display parameters given through signal.json and PSO-estimated parameters
    disp(['Estimated parameters: m1=',num2str(outStruct.bestGwCoefs(1)),...
                              '; m2=',num2str(outStruct.bestGwCoefs(2))]);
end
if jobParams.injSig == 1
    disp(['Original parameters: A = ',num2str(signalParams.snr),...
        '; phi = ',num2str(signalParams.phase),...
        '; t_a = ',num2str(signalParams.ta),...
        '; FitVal = ',num2str(original_fitVal)]);
end
disp(['Estimated parameters: A = ',num2str(outStruct.bestAmp),...
                                 '; phi = ',num2str(outStruct.bestPhase),...
                                  '; t_a = ',num2str(outStruct.bestTime),...
                                  '; FitVal = ',num2str(bestFitVal)]);
fprintf(fidprog,'\nPost-processing...done');
fclose(fidprog);