function comb_anashpsjb(path2jsonlab,jobParamsFile,userUID,varargin)
%COMBANASHPSJB(J,F,U)
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
%defined in GENANASHPSLNCHRJB. Each output file contains the processed
%data matrices from Drase. The separate processed data matrices are
%combined to produce the unified SHAPES processed data matrices of the
%original input data file. The combined data matrices are stored out in a
%file in the folder containing all the output files.
%
%U is the userUID that was set for the GENANASHPSLNCHRJB.m


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
                combFileName = varargin{lpargs};
            case 2
                datad = varargin{lpargs};
        end
    end
end

[~,outdataFilePrfx] = ana_basics(jobParams,userUID,datad,1);
dataFile = [outdataFilePrfx,'C'];
shpsDataFile = [outdataFilePrfx,'shps_C'];
% paramsFileshps = [paramsFile,'shps'];
%Construct name of the .txt file containing the list of output file names
% outFilesListFile = [outdataFilePrfx,'outFilesList.txt'];
% fidFileList = fopen(outFilesListFile,'r');
% if fidFileList == -1
%     error(['Error opening ',outFilesListFile]);
% end
% nFiles = 0;
% outFilesList = {};
% % while ~feof(fidFileList)
% %     outFileName = fgetl(fidFileList);
% %     outFilesList = [outFilesList, outFileName];
% %     nFiles = nFiles+1;
% % end
% fclose(fidFileList);

%Construct name of the file to store the combined data
combFileName = [outdataFilePrfx,'F.mat'];



%% Data 
%Input Data
inputData = load(jobParams.inFileData);

%PSD Data and interpolated PSD data
psdData = load(jobParams.inFilePSD); %pwelch psd and highpassed time series
estpsdData = load(jobParams.inFileshpsPSD); %estimated PSD and sampFreq only

%Rungwpso Data
outData = load(dataFile);
estoutData = load(shpsDataFile);

save(combFileName,"inputData","psdData","estpsdData","outData","estoutData");
[pathstr,filename,~] = fileparts(combFileName);
disp(['Saved ',filename,' to ',pathstr])
end