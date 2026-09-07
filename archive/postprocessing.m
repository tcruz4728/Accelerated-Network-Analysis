function postprocessing(inDataFile,filepaths,varargin)
% POSTPROCESSING(F,P,I,C)
% Loads in a file F which contains 5 data structures corresponding to the
% data as it transitioned from 
% 'inputData' -  time series data and time interval
% 'psdData' - contains pwelch data and interpolated PSD
% 'estpsdData' - contains shapes estimated PSD
% 'outData' - output data from matched filtering run
% 'estoutData' - output data from matched filtering run on SHAPES estimated
% data
% P is a filepath structure; the parameter text file containing various
% matched filtering values (amplitude, phase, etc.) is saved to P.tables.
% Alongside it a .mat file is saved with the same information for ease of
% manipulation in a structure array with these additional fields.


% Optional input arguments
injSig = [];
realizationCount = 1;
nreqArgs = 2;
dsstFileName = [];
for lpargs = 1:(nargin-nreqArgs)
    if ~isempty(varargin{lpargs})
        switch lpargs
            case 1
                injSig = varargin{lpargs};
            case 2
                realizationCount = varargin{lpargs};
            case 3
                dsstFileName = varargin{lpargs};
                load(dsstFileName,'dsstPSD','dsstfreqVec')
        end
    end
end

%Variable Loading
if iscell(inDataFile)
    inputData = inDataFile{:,1};
    psdData = inDataFile{:,2};
    estpsdData = inDataFile{:,3};
    outData = inDataFile{:,4};
    estoutData = inDataFile{:,5};
else
    load(inDataFile,'inputData','psdData','estpsdData','outData','estoutData')
end

%% Time-Series Plots
figure(1);
if realizationCount == 1
    tsd = tiledlayout('flow');
    title(tsd,'Time series data')

    nexttile
else
    nexttile
end
hold on
plot(outData.params.dataX,inputData.dataY)
title(['Realization: ', num2str(realizationCount)])
xlabel('Time(s)')
ylabel('Amplitude')
axis tight
hold off
%% PSD plots
%Cutoff Frequency index update
initialCut = psdData.freqBnd(1,1);
freqVecIndxco = find(psdData.freqVec<=initialCut,1,'last');
freqVecCut = psdData.freqVec(freqVecIndxco:end);
figure(2);
if realizationCount == 1
    psdPlot = tiledlayout('flow');
    title(psdPlot,'Logarithmic PSD')
    nexttile
else
    nexttile
end
hold on
plot(freqVecCut,log10(psdData.PSD(freqVecIndxco:end)),'g',...
    'DisplayName',['Pwelch ',num2str(realizationCount)]);
plot(freqVecCut,log10(estpsdData.PSD(freqVecIndxco:end)),'b',...
    'DisplayName',['SHAPES Est ',num2str(realizationCount)]);
if ~isempty(dsstFileName)
    plot(dsstfreqVec,log10(dsstPSD),'k--',...
        'DisplayName','Design Sensitivity')
end
axis tight
title(['Realization: ', num2str(realizationCount)])
xlabel('Frequency (Hz)')
ylabel('Amplitude Spectrum (1/sqrt(Hz))')
legend([],'Location','northeast')
% saveas(gcf,[filepaths.figure,'PSD']);
hold off

freqVecinterp = linspace(psdData.freqVec(1),psdData.freqVec(end),length(psdData.interpPSD));
freqVecinterpIndxco = find(freqVecinterp<=initialCut,1,'last');
freqVecinterpCut = freqVecinterp(freqVecinterpIndxco:end);
figure(3);
if realizationCount == 1
    interppsdPlot = tiledlayout('flow');
    title(interppsdPlot,...
        'Interpolated Logarithmic PSD')
    nexttile
else
    nexttile
end
hold on
plot(freqVecinterpCut,log10(psdData.interpPSD(freqVecinterpIndxco:end)),'g',...
    'DisplayName',['Pwelch ',num2str(realizationCount)]);
plot(freqVecinterpCut,log10(estpsdData.interpPSD(freqVecinterpIndxco:end)),'b',...
    'DisplayName',['SHAPES Est ',num2str(realizationCount)]);
if ~isempty(dsstFileName)
    plot(dsstfreqVec,log10(dsstPSD),'k--',...
        'DisplayName','Design Sensitivity')
end
axis tight
title(['Realization: ', num2str(realizationCount)])
xlabel('Frequency (Hz)')
ylabel('Amplitude Spectrum (1/sqrt(Hz))')
legend([],'Location','northeast')
hold off
% saveas(gcf,[filepaths.figure,'PSD_Interpolated']);

%% Data type loop
%Records parameters in txt file
fidparams = fopen([filepaths.tables,filesep,'parameters.txt'],'a');
fprintf(fidparams,'%s\n',char(datetime("today")));
fprintf(fidparams,'%s\n',['Realization: ', num2str(realizationCount)]);
for datatype = 1:2
    switch datatype
        case 1
            % load(dataFile,"params","outStruct","bestFitVal","original_fitVal")
            dataStruct = outData;
            % disp('Running on Welch Data')
            titlestr = 'Pwelch Data';
            legendstr = ['Pwelch ',num2str(realizationCount)];
            linespecstr = '-';
        case 2
            % load(shpsDataFile,"params","outStruct","bestFitVal","original_fitVal")
            dataStruct = estoutData;
            % disp('Running on Shapes Estimated Welch Data')
            titlestr = 'Shapes Estimated Data';
            legendstr = ['SHAPES Est ',num2str(realizationCount)];
            linespecstr = '--';
    end
    fprintf(fidparams,'%s\n',titlestr);
    params = dataStruct.params;
    outStruct = dataStruct.outStruct;
    bestFitVal = dataStruct.bestFitVal;
    original_fitVal = dataStruct.original_fitVal;
    mf1 = dataStruct.mf(1,:);
    mf2 = dataStruct.mf(2,:);
    gwCoefs = params.gwCoefs;
    signalParams = params.signal;
    %% Matched Filtering quadrature
    figure(4);
    if realizationCount == 1 && datatype == 1
        mfq = tiledlayout('flow');
        title(mfq,'Matched Filtering Quadrature')
        nexttile
    elseif datatype == 1
        nexttile
    end
    title(['Realization: ', num2str(realizationCount)])
    hold on
    plot(params.dataX,sqrt(mf1.^2 + mf2.^2),'DisplayName',legendstr)
        legend([])
        xlabel('Time(s)')
    axis tight
    hold off
    
    %% Time Series plot with Matched Filtered signals overlayed
    figure(5);
    if realizationCount == 1 && datatype == 1
        sts = tiledlayout('flow');
        title(sts,'Signal Time Series')
        nexttile
    elseif datatype == 1
        nexttile
    end
   
    hold on;
    title(['Realization: ', num2str(realizationCount)])
    plot(params.dataX,outStruct.bestSig,'DisplayName',legendstr)
    %     plot(params.dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
    xlabel('Time(s)')
    ylabel('Amplitude')
    legend([],'Location','southwest')
    axis tight
        % saveas(gcf,[filepaths.figure,'PSO_Results']);
    hold off
    %% Allbest Plots
    figure(6);
    hold on
    plot(outStruct.allBestFitness,linespecstr,'DisplayName',legendstr)
    title('Best Fitness and Locations vs. Iterations')
    legend([])
    ylabel('Fitness Value')
    xlabel('Iteration Number')
    hold off

    figure(7);
    hold on
    plot(outStruct.allBestLocation(:,1),linespecstr,'DisplayName',legendstr)
    title('Best x-Location vs. Iterations')
    legend([])
    ylabel('X-Coordinate')
    xlabel('Iteration Number')
    hold off

    figure(8);
    hold on
    plot(outStruct.allBestLocation(:,2),linespecstr,'DisplayName',legendstr)
    title('Best y-Location vs. Iterations')
    legend
    ylabel('Y-Coordinate')
    xlabel('Iteration Number')
    hold off
    %% GW Coefficients Iteration Optimization
    if datatype == 1
        figure(9);
        boundary_plot(gwCoefs);
        hold off;
    end

    t0 = outStruct.bestGwCoefs(1);
    t1p5 = outStruct.bestGwCoefs(2);
    est_M = (5/(32*params.frange(1)))*(t1p5/(pi*pi*t0))*(params.cgFac);
    est_u = (1/(16*params.frange(1)^2))*(5/(4*pi^4*t0*t1p5^2))^(1/3)*(params.cgFac);
    est_m1 = (est_M - sqrt(est_M^2 - 4*est_u*est_M))/2;
    est_m2 = (est_M + sqrt(est_M^2 - 4*est_u*est_M))/2;
%GW coefficients
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
%Parameters
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
    fprintf(fidparams,'%s\n\n',estParams);

    % Create accompanying structure with parameters
    gwResults = struct();
    % Injected (true) signal parameters if available
    if injSig == 1
        gwResults.injected = struct(...
            'tau0', gwCoefs(1), ...
            'tau1p5', gwCoefs(2), ...
            'm1_kg', signalParams.masses(1), ...
            'm2_kg', signalParams.masses(2), ...
            'm1_Msun', signalParams.masses(1)/Msolar, ...
            'm2_Msun', signalParams.masses(2)/Msolar, ...
            'A', signalParams.snr, ...
            'phase', signalParams.phase, ...
            'ta', signalParams.ta, ...
            'fitVal', original_fitVal);
    else
        gwResults.injected = [];
    end
    % Estimated results from PSO/SHAPES
    gwResults.estimated = struct(...
        'tau0', outStruct.bestGwCoefs(1), ...
        'tau1p5', outStruct.bestGwCoefs(2), ...
        'm1_kg', est_m1, ...
        'm2_kg', est_m2, ...
        'm1_Msun', est_m1/Msolar, ...
        'm2_Msun', est_m2/Msolar, ...
        'A', outStruct.bestAmp, ...
        'phase', outStruct.bestPhase, ...
        'ta', outStruct.bestTime, ...
        'fitVal', bestFitVal, ...
        'gwCoefs_raw', outStruct.bestGwCoefs);

    % Attach metadata
    gwResults.metadata = struct(...
        'realization', realizationCount, ...
        'dataFile', inDataFile, ...
        'estimationType', titlestr, ...
        'timestamp', char(datetime('now')));

    % Save into results structure for output
    if ~exist('results','var') || ~isstruct(results)
        results = struct();
    end
    if ~isfield(results,'GW')
        results.GW = {};
    end
    results.GW{end+1} = gwResults;

end
fclose(fidparams);
% Save GW coefficients and parameters to a structure
save([filepaths.tables,'parameters.mat'],"results")
disp(['postprocessing- Saved parameters to .txt and .mat files: ',filepaths.tables,'parameters.txt and parameters.mat'])
end