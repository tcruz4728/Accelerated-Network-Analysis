function combmultimfshps(path2jsonlab,jobParamsFile,userUID,varargin)
%COMBMULTIMFSHPS(J,F,U)
%Post-job function to combine data for PC download.
%J is the path to the jsonlab package. If set to '', the jsonlab package is
%assumed to be in the Matlab search path. Note that this search path must
%be accessible from every compute node, which cannot always be guaranteed.
%Hence it is safer to specify the path explicitly.
%
%F is a JSON file containing the parameters of the job whose output is to
%be post-processed.
%
%The file F is read to obtain the name of the folder containing all the
%output files as well as a .txt file containing the list of the output file
%names. The convention followed for the naming of the required files is as
%defined in GENMULTISHPSLNCHRJB. Each output file contains the processed
%data matrices from Drase. The separate processed data matrices are
%combined to produce the unified SHAPES processed data matrices of the
%original input data file. The combined data matrices are stored out in a
%file in the folder containing all the output files.
%
%U is the userUID that was set for the GENMULTISHPSLNCHRJB.m


% Modified from combsplitlinesshps to work with ANA data

%Add path to JSONLAB
addpath(path2jsonlab);

%Load job description
jobParams = loadjson(jobParamsFile);



datad = date;
%Override the file name if optional input given
nreqArgs = 3;
for lpargs = 1:(nargin-nreqArgs)
    if ~isempty(varargin{lpargs})
        switch lpargs
            case 1
                datad = varargin{lpargs};
        end
    end
end

[~,outdataFilePrfx] = ana_prep(jobParams,userUID,datad,[],1);
% paramsFileshps = [paramsFile,'shps'];

%Nx3 Cell Array with inFileData, inFilePSD, and inFileshpsPSD paths
inFilesList = readcell([outdataFilePrfx,'_inFilesList.txt'],"Delimiter","  ");
%Nx2 Cell Array with dataFile and shpsdataFile paths
outFilesList = readcell([outdataFilePrfx,'_outFilesList.txt'],"Delimiter","  ");
%Nx2 Cell Array with params and paramsshps paths
% paramsFilesList = readcell([outdataFilePrfx,'_paramsFilesList.txt'],"Delimiter","  ");

%Time series initial file
inFileNameList = inFilesList{:,1};
%Pwelch data file (inFilePSD-goes into SHAPES)
interFileNameList = inFilesList{:,2};
%SHAPES estimated PSD data file (inFileshpsPSD)
outFileNameList = inFilesList{:,3};

% %Pwelch params files
% paramsFileList = paramsFilesList{:,1};
% %SHAPES params files
% paramsFileshpsList = paramsFilesList{:,2};

%Outgoing data files
dataFileList = outFilesList{:,1};
shpsDataFileList = outFilesList{:,2};

%% Data Loop 
combData = cell(jobParams.inFileDataRange(2),5);
for fileCount = 1:jobParams.inFileDataRange(2)
    %Input Data
    inputData = load(inFileNameList{fileCount});

    %PSD Data and interpolated PSD data
    psdData = load(interFileNameList{fileCount}); %pwelch psd and highpassed time series
    estpsdData = load(outFileNameList{fileCount}); %estimated PSD and sampFreq only

    %Rungwpso Data
    outData = load(dataFileList{fileCount});
    estoutData = load(shpsDataFileList{fileCount});
    combData{fileCount} = {inputData,psdData,estpsdData,outData,estoutData};
    save([outdataFilePrfx,'_n',num2str(fileCount),'F.mat'],...
        "inputData","psdData","estpsdData","outData","estoutData");
    disp(['combmultimfshps- Saved individual run to ',...
        outdataFilePrfx,'_n',num2str(fileCount),'F.mat'])
end
save([outdataFilePrfx,'F.mat'],"combData");
disp(['combmultimfshps- Saved combined runs to ',...
        outdataFilePrfx,'F.mat'])
end