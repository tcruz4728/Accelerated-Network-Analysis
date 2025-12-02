function genanarealization(path2json,outDir)
% Documentation for the genanarealization function
% 
% This function generates data realizations based on specified parameters 
% and saves the results to a designated output directory. The function 
% supports two types of data generation: loading from externally generated 
% files or creating data from sensitivity data.
%
% Inputs:
%   inFileParams - A JSON file path containing parameters for the 
%                  generation process, including paths to other JSON files 
%                  for PSO and signal parameters.
%   varargin - Optional argument to specify the output directory. If not 
%              provided, the output directory is taken from the input 
%              parameters.
%
% Outputs:
%   The function saves generated data realizations as .mat files in the 
%   specified output directory. Each file contains the generated data, 
%   time interval, power spectral density (PSD), frequency vector, and 
%   injected signal parameters.
%
% Example usage:
%   genanarealization('params.json', 'output_directory')
addpath(path2json)
params = loadjson('realizations.json');
psoParams = loadjson('pso.json');
signalParams = loadjson('signal.json');

%% Data Load
switch params.dataGenType
    case 1 %Data Loading from externally generated files
        load(params.extGenFilename,"data_realizations","sigparams")  
        dsstPSD = sigparams.PSD;
        kNyq = floor(length(dsstPSD)/2)+1;
        dsstfreqVec = (0:(kNyq-1))*sigparams.sampFreq/length(dsstPSD);
    case 2 %Data Created from sensitivity data
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

%% File Saving
for i = 1:size(data_realizations,1)
    dataY = data_realizations(i,:);
    save([outDir,'TMPPSDDATA',filesep,...
        'realization','_',num2str(sigparams.datalen),...
        's_inj',num2str(sigparams.ta),'_fs',num2str(sigparams.sampFreq),...
        '_n',num2str(i),'.mat'],...
        "dataY","tIntrvl","dsstPSD","dsstfreqVec","injSigparams")
end
disp(['dataRealizationgenmat- Saved ',num2str(size(data_realizations,1)),...
    ' data files to ',filesep,outDir,filesep,'TMPPSDDATA',...
    filesep,...
    ' as realization_',num2str(sigparams.datalen),...
        's_inj',num2str(sigparams.ta),...
        '_fs',num2str(sigparams.sampFreq),'_n<#>','.mat'])
end

