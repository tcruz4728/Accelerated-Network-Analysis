% function postprocessing(jobParams,dataFile,shpsDataFile,filepaths)
function postprocessing(inDataFile,filepaths,injSig)
% Loads in a file which contains 5 data structures.
% 'inputData' -  time series data and time interval
% 'psdData' - contains pwelch data and interpolated PSD
% 'estpsdData' - contains shapes estimated PSD
% 'outData' - output data from matched filtering run
% 'estoutData' - output data from matched filtering run on SHAPES estimated
% data

load(inDataFile,'inputData','psdData','estpsdData','outData','estoutData')

%PSDs from each
% load(jobParams.inFilePSD,'PSD'); %Welchs PSD
% pwelchPSD = PSD;
% load(jobParams.inFileshpsPSD,'PSD'); %SHAPES estimated Welchs PSD
% shpsPSD = PSD;
%% Time Series
figure; 
plot(outData.params.dataX,inputData.dataY)
title('Time series data')
%% PSD plots
figure;
semilogy(psdData.freqVec,psdData.PSD)
hold on
semilogy(estpsdData.freqVec,estpsdData.PSD)
axis tight
title('Pwelch and SHAPES Estimated Pwelch data')
legend('Pwelch','SHAPES Est')
saveas(gcf,[filepaths.figs,'PSD']);
hold off

freqVecinterp = linspace(psdData.freqVec(1),psdData.freqVec(end),length(psdData.interpPSD));
figure;
semilogy(freqVecinterp,psdData.interpPSD)
hold on
semilogy(freqVecinterp,estpsdData.interpPSD)
axis tight
title('Interpolated Pwelch and SHAPES Estimated Pwelch data')
legend('Pwelch','SHAPES Est')
saveas(gcf,[filepaths.figs,'PSD_Interpolated']);
hold off

%Records parameters in txt file
fidparams = fopen([filepaths.end,'parameters.txt'],'a');
fprintf(fidparams,'%s\n',datestr(datetime('now')));
for i = 1:2
    switch i
        case 1
            % load(dataFile,"params","outStruct","bestFitVal","original_fitVal")
            dataStruct = outData;
            disp('Running on Welch Data')
            titlestr = 'Pwelch Data';

        case 2
            % load(shpsDataFile,"params","outStruct","bestFitVal","original_fitVal")
            dataStruct = estoutData;
            disp('Running on Shapes Estimated Welch Data')
            titlestr = 'Shapes Estimated Data';
    end
    fprintf(fidparams,'%s\n',titlestr);
    params = dataStruct.params;
    outStruct = dataStruct.outStruct;
    bestFitVal = dataStruct.bestFitVal;
    original_fitVal = dataStruct.original_fitVal;
    % load(paramsFile,"gwCoefs")
    gwCoefs = params.gwCoefs;
    %     psoParams = params.pso;
    signalParams = params.signal;
    %% Time Series plot with Matched Filtered signals overlayed
    figure;
    hold on;
    plot(params.dataX,params.dataY,'.');
    title(titlestr)
    plot(params.dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
    legend('Data','Signal',...
        'Estimated signal: Best run');
    axis tight
    saveas(gcf,[filepaths.figs,'PSO_Results']);
    hold off;

    %% GW Coefficients Iteration Optimization
    if i == 1
        figure;
        hold on;
        scatter(gwCoefs(1),gwCoefs(2),140,'red','filled','D','DisplayName','Original Parameters');
        title("Best Parameter Values for All Runs");
        %     if psoParams.type == 2
        xlabel("\tau_0");
        ylabel("\tau_{1.5}");
        legend;
        boundary_plot;
        hold off;
    end

    t0 = outStruct.bestGwCoefs(1);
    t1p5 = outStruct.bestGwCoefs(2);
    est_M = (5/(32*params.frange(1)))*(t1p5/(pi*pi*t0))*(params.cgFac);
    est_u = (1/(16*params.frange(1)^2))*(5/(4*pi^4*t0*t1p5^2))^(1/3)*(params.cgFac);
    est_m1 = (est_M - sqrt(est_M^2 - 4*est_u*est_M))/2;
    est_m2 = (est_M + sqrt(est_M^2 - 4*est_u*est_M))/2;

    if injSig == 1
        injsigGWcoefs = ['Injected Signal GW Coefficients: tau0= ',num2str(gwCoefs(1)),...
            '; tau1p5= ',num2str(gwCoefs(2)),...
            '; m1= ', num2str(signalParams.masses(1)),...
            '; m2= ', num2str(signalParams.masses(2))];
        fprintf(fidparams,'%s\n',injsigGWcoefs);
    end
    Msolar = 1.989*10^30; %Solar mass in kg
    % This will display parameters given through signal.json and PSO-estimated parameters
    estGWcoefs = ['Estimated GW Coefficients: tau0=',num2str(outStruct.bestGwCoefs(1)),...
        '; tau1p5=',num2str(outStruct.bestGwCoefs(2)),...
        '; m1= ', num2str(est_m1/Msolar),...
        '; m2= ', num2str(est_m2/Msolar)];
    fprintf(fidparams,'%s\n',estGWcoefs);
    %     else
    %         xlabel("m_1");
    %         ylabel("m_2");
    %         legend;
    %         hold off;
    %         % This will display parameters given through signal.json and PSO-estimated parameters
    %         estmParams = ['Estimated mass parameters: m1=',num2str(outStruct.bestGwCoefs(1)),...
    %             '; m2=',num2str(outStruct.bestGwCoefs(2))];
    %         fprintf(fidparams,'%s\n',estmParams);
    %     end
    if injSig == 1
        injsigParams = ['Injected Signal parameters: A = ',num2str(signalParams.snr),...
            '; phi = ',num2str(signalParams.phase),...
            '; t_a = ',num2str(signalParams.ta),...
            '; FitVal = ',num2str(original_fitVal)];
        fprintf(fidparams,'%s\n',injsigParams);
    end
    estParams = ['Estimated parameters: A = ',num2str(outStruct.bestAmp),...
        '; phi = ',num2str(outStruct.bestPhase),...
        '; t_a = ',num2str(outStruct.bestTime),...
        '; FitVal = ',num2str(bestFitVal)];
    fprintf(fidparams,'%s\n',estParams);
end
fclose(fidparams);
end