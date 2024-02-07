% anabasics - bare bones from initial setup and before call to rungwpso
function [paramsFile,outdataFilePrfx,varargout] = ana_basics(jobParams,userUID,varargin)
% Optional input Arguments 
% [P,Y] = ana_basics(J,U,D,E,C,T)
% Setting D allows the user to set a specific date for folder creation
% This could be used to pull data from a specific date if E is not empty,
% or overwride pre-existing data if it is. E, if specified, returns the
% function after only writing file strings for P and Y.
% C controls whether or not a progress.txt file is generated for the run, 1 is yes, else is no.
% Optional outputs. T is the plotting control option, default is no plot, 
% setting '1' generates a spectrogram of the training series, '2' generates 
% the PSD curve and '12' generates both.
% Optional output Arguments
% [P,Y,F,N] = ana_basics(J,U,D,E,C,T)
% F returns a structure of file paths directories for figures, PSDs, and
% end files. N is the full filename of the progress.txt file. 
%% jobParams datatype check - checks for whether jobParams is a structure, else assumes a file
if ~isstruct(jobParams)
    jobParams = loadjson(jobParams);
end

%% Optional Argument
nreqArgin = 2;
datad = [];
anabreak = [];
progCtrl = [];
pltCtrl = 0;
for largs = 1:(nargin-nreqArgin)
    if ~isempty(varargin{largs})
        switch largs
            case 1
                datad = varargin{largs};
            case 2
                anabreak = varargin{largs};
            case 3
                progCtrl = varargin{largs};
            case 4
                pltCtrl = varargin{largs};
        end
    end
end

%% Initial Setup: Parameters - job settings
%Project and DRASE function and JSON load
addpath(genpath(jobParams.path2drase));

% Defining File Paths/Names and Folder Creation
[outDir,filepaths] = dpfc(jobParams,userUID,datad);
varargout{1} = filepaths;
paramsFile = [outDir,'params']; %rungwpso params file
paramsFileshps = [paramsFile,'shps']; %rungwpso params file for shapes data

%Project Parameters
psoParams = loadjson(jobParams.psoParamsjson); %matched filtering PSO params
signalParams = loadjson(jobParams.signalParamsjson); %signal injection params
%%%NOTE%%% - additional parameters for drase.m are called in
%%% via jobParams.drasejobParamsFile

%File Prefix Generation - Names files with project-specific parameters
filetagstr = filetagana(psoParams,signalParams); 
outdataFilePrfx = [outDir,jobParams.jobName,'_',filetagstr];

% Quick stop for partial runs, namely for setting file conventions
if ~isempty(anabreak)
    return
end

%% Progress Text file - Monitor code progress and completion (optional)
if ~isempty(progCtrl) && progCtrl == 1
    progressFile = [outDir,'progress.txt'];
    varargout{2} = progressFile;
    fidprog = fopen(progressFile,'a');
    disp(['Progress File created in ',outDir])
    proglines = struct('nd','done.',...
    't',datestr(datetime('now')),... %fprintf(fidprog,'%s\n',datestr(datetime('now')));
    'c','Computing Parameters...',...
    'l','Loading Data...',...
    'd','Running Drase...',...
    'w','Data Conditioning on Welch...',...
    's','Data Conditioning on SHAPES...',...
    'g','Glitch Check...');
    progstatus(proglines.t,fidprog,progCtrl)
end
%% GW Parameter settings - combines relevant settings and performs necessary 
% computations for rungwpso. Does not require time-series or PSD data
    progstatus(proglines.c,fidprog,progCtrl)
gwpsoparams(psoParams,signalParams,paramsFile);
copyfile([paramsFile,'.mat'],[paramsFileshps,'.mat']); %identical but separate parameter settings
    progstatus(proglines.nd,fidprog,progCtrl)
%% Data Load - loads time-series performs bandpass, computes training segment PSD
% input: inFileData - time series data from LIGO or simulations
% output: inFilePSD - training segment PSD
    progstatus(proglines.l,fidprog,progCtrl)
outData = load_mtchdfltrdata(jobParams.inFileData,jobParams.inFilePSD);
    progstatus(proglines.nd,fidprog,progCtrl)
%% Glitch Checking - Visual clarification with spectrogram
if pltCtrl == 1 || pltCtrl == 12
    progstatus(proglines.g,fidprog,progCtrl)
    figure;
    plot(outData.dataY)
    [S,F,T] = spectrogram(outData.tseriestrainSeg,8192,8000,[],outData.sampFreq);
    S = abs(S);
    imagesc(T,F,log10(S)); axis xy; %Checking spectrogram image for glitches or high noise
    title('Training Segment Spectrogram')
    saveas(gcf,[filepaths.figs,'Training_Spectrogram']);
    progstatus(proglines.nd,fidprog,progCtrl)
end
%% PSDs training segment plot
if pltCtrl == 2 || pltCtrl == 12
    figure;
    semilogy(outData.freqVec,outData.PSD); axis tight
    title('PSD of Training Segment')
    saveas(gcf,[filepaths.figs,'Training_PSD']);
end
%% SHAPES PSD estimate - takes the pwelch PSD (not the log version!) and returns in the
%same form
% input: inFilePSD - training segment PSD from load_mtchdfltrdata.m
% output: inFileshpsPSD - shapes estimation of training segment PSD
    progstatus(proglines.d,fidprog,progCtrl)
run_drase4lines(jobParams,userUID);
    progstatus(proglines.nd,fidprog,progCtrl)

%% Interpolation - Takes log10 of PSDs, interpolates and inverses the log 
%input: outData/inFilePSD - data structure from load_mtchdfltrdata.m
%output: inFilePSD - appending the interpolated PSD
createPSD(outData.PSD,outData.freqVec,outData.tlen,outData.sampFreq,jobParams.inFilePSD);
%input: inFileshpsPSD - estimated PSD from drase
%output: inFileshpsPSD - appending the interpolated shpsPSD
load(jobParams.inFileshpsPSD,"PSD")
createPSD(PSD,outData.freqVec,outData.tlen,outData.sampFreq,jobParams.inFileshpsPSD);

%% Condition Data and Compute FFTs
%inputs: inFilePSD - interpolated PSD highpassed time series
%output: paramsFile - updated from gwpsoparams with fft
    progstatus(proglines.w,fidprog,progCtrl)
cond_mtchdfltrdata(jobParams.inFilePSD,paramsFile,jobParams.injSig);
    progstatus(proglines.nd,fidprog,progCtrl)
%inputs: inFileshpsPSD - interpolated shps PSD & inFilePSD - highpassed time series
%output: paramsFileshps - updated from gwpsoparams with fft
    progstatus(proglines.s,fidprog,progCtrl)
cond_mtchdfltrdata(jobParams.inFileshpsPSD,paramsFileshps,jobParams.injSig);
    progstatus(proglines.nd,fidprog,progCtrl)
[paramsFilepath,paramsFileName,~] = fileparts([paramsFile,'.mat']);
[~,paramsFileshpsName,~] = fileparts([paramsFileshps,'.mat']);
disp(['Parameter files ',paramsFileName, ' and ',paramsFileshpsName,' saved into ',paramsFilepath])
if ~isempty(progCtrl) && progCtrl == 1, fclose(fidprog); end
end