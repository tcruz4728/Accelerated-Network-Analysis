params = loadjson(['JSON',filesep,'realizations.json']);

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

for i = 1:size(data_realizations,1)
    dataY = data_realizations(i,:);
    save(['TMPPSDDATA',filesep,'realization','_',num2str(sigparams.datalen),...
        's_inj',num2str(sigparams.ta),'_fs',num2str(sigparams.sampFreq),'_n',num2str(i),'.mat'],...
        "dataY","tIntrvl",'dsstPSD',"dsstfreqVec")
end
disp(['dataRealizationgenmat- Saved ',num2str(size(data_realizations,1)),...
    ' data files to Accelerated-Network-Analysis',filesep,'TMPPSDDATA',...
    filesep,...
    ' as realization_',num2str(sigparams.datalen),...
        's_inj',num2str(sigparams.ta),...
        '_fs',num2str(sigparams.sampFreq),'_n<#>','.mat'])
