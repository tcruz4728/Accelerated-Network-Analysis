params = loadjson(['JSON',filesep,'realizations.json']);

switch params.dataGenType
    case 1 %Data Loading from externally generated files
        load(params.extGenFilename)  
        PSD = sigparams.PSD;
        frange = sigparams.frange;
    case 2 %Data Created from sensitivity data
        % userNrealizations = input("Set the number of data realizations to be generated: ");
        sigparams = params;
        addpath(sigparams.path2sensitivitydata)
        Nsamples = sigparams.sampFreq*sigparams.datalen;
        data_realizations = zeros(sigparams.nRealizations,Nsamples);
        for lpruns = 1:sigparams.nRealizations
            [data_realizations(lpruns,:),PSD] = LIGOnoise(Nsamples,sigparams.sampFreq,1,'sample');
        end
end
tIntrvl = 1/sigparams.sampFreq;

for i = 1:size(data_realizations,1)
    dataY = data_realizations(i,:);
    save(['TMPPSDDATA',filesep,'test','_',num2str(sigparams.datalen),...
        's_inj',num2str(sigparams.ta),'_fs',num2str(sigparams.sampFreq),'_n',num2str(i),'.mat'],...
        "dataY","tIntrvl",'PSD')
end
disp(['dataRealizationgenmat- Saved ',num2str(size(data_realizations,1)),...
    ' data files to Accelerated-Network-Analysis',filesep,'TMPPSDDATA',...
    filesep,...
    ' as test_',num2str(sigparams.datalen),...
        's_inj',num2str(sigparams.ta),...
        '_fs',num2str(sigparams.sampFreq),'_n<#>','.mat'])
