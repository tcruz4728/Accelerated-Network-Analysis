function combmultimfshps(outdataFilePrfx,inFileList,outFileList)
%COMBMULTIMFSHPS(D,I,O)
%Post-job function to combine data for PC download.
% D- the output file data prefix given to the files that are saved.
% I- the input file list used in the 1st ls6 job with SHAPES, contains N
% lines corresponding to the N files with data realizations. Needs to match
% the file created by genmultishpslnchrjb.m
% O- the output file list used in the 2nd ls6 job with rungwpso, contains N
% lines corresponding to the N files with PSD estimates. Needs to match the
% file created by genmultipsolnchrjb.m
%
% Modified from combsplitlinesshps to work with ANA data
%
% See also 

%Nx3 Cell Array with inFileData, inFilePSD, and inFileshpsPSD paths
inFilesList = readcell(inFileList,"Delimiter","  ");
%Nx2 Cell Array with dataFile and shpsdataFile paths
outFilesList = readcell(outFileList,"Delimiter","  ");

nFiles = min(size(outFilesList,1),size(inFilesList,1));

%Time series initial file
inFileNameList = inFilesList(nFiles,1);
%Pwelch data file (inFilePSD-goes into SHAPES)
interFileNameList = inFilesList(nFiles,2);
%SHAPES estimated PSD data file (inFileshpsPSD)
outFileNameList = inFilesList(nFiles,3);

%Outgoing data files
dataFileList = outFilesList(nFiles,1);
shpsDataFileList = outFilesList(nFiles,2);

%% Data Loop 
combData = cell(nFiles,5);

for fileCount = 1:nFiles
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