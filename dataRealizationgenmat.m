params = loadjson(['JSON',filesep,'realizations.json']);
psoParams = loadjson(['JSON',filesep,'pso.json']);
signalParams = loadjson(['JSON',filesep,'signal.json']);
%% Data Load
switch params.dataGenType
    case 1 %Data Loading from externally generated files
        load(params.extGenFilename)  
        dsstPSD = sigparams.PSD;
        frange = sigparams.frange;
    case 2 %Data Created from sensitivity data
        % userNrealizations = input("Set the number of data realizations to be generated: ");
        sigparams = params;
        addpath(sigparams.path2sensitivitydata)
        Nsamples = sigparams.sampFreq*sigparams.datalen;
        data_realizations = zeros(sigparams.nRealizations,Nsamples);
        for lpruns = 1:sigparams.nRealizations
            [data_realizations(lpruns,:),dsstPSD,dsstfreqVec] = LIGOnoise(Nsamples,sigparams.sampFreq,1,'sample');
        end
end
tIntrvl = 1/sigparams.sampFreq;

%% Injected Signal Creation
negFStrt = 1-mod(Nsamples,2);
kNyq = floor(Nsamples/2)+1;
% Compute two-sided PSD from design sensitivity PSD for signal injection
dsstPSDtotal = [dsstPSD, dsstPSD((kNyq-negFStrt):-1:2)];
injSigparams = gwpsoparams(psoParams,signalParams,0);
injSigparams.signal.data = sigInj(injSigparams,dsstPSDtotal);

for i = 1:size(data_realizations,1)
    dataY = data_realizations(i,:);
    save(['TMPPSDDATA',filesep,'realization','_',num2str(sigparams.datalen),...
        's_inj',num2str(sigparams.ta),'_fs',num2str(sigparams.sampFreq),'_n',num2str(i),'.mat'],...
        "dataY","tIntrvl","dsstPSD","dsstfreqVec","injSigparams")
end
disp(['dataRealizationgenmat- Saved ',num2str(size(data_realizations,1)),...
    ' data files to Accelerated-Network-Analysis',filesep,'TMPPSDDATA',...
    filesep,...
    ' as realization_',num2str(sigparams.datalen),...
        's_inj',num2str(sigparams.ta),...
        '_fs',num2str(sigparams.sampFreq),'_n<#>','.mat'])
