dataGenType = input("Data Generation Type: ");
switch dataGenType
    case 1 %Data Loading from externally generated files
        load("data_realizations_10Nov2023.mat")  
        PSD = sigparams.PSD;
        frange = sigparams.frange;
    case 2 %Data Created from sensitivity data
        % userNrealizations = input("Set the number of data realizations to be generated: ");
        sigparams = loadjson('JSON\realizations.json');
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
    save(['TMPPSDDATA\test',num2str(i),'_',num2str(sigparams.datalen),...
        's_inj',num2str(sigparams.ta),'_fs',num2str(sigparams.sampFreq),'.mat'],...
        "dataY","tIntrvl",'PSD')
end
disp(['Saved ',num2str(size(data_realizations,1)),' data files to Accelerated-Network-Analysis\TMPPSDDATA\ as test<#>_',num2str(sigparams.datalen),...
        's_inj',num2str(sigparams.ta),'_fs',num2str(sigparams.sampFreq),'.mat'])
