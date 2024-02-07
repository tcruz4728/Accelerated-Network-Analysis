function [paramsFile,outdataFilePrfx,varargout] = ana_prep(jobParams,userUID,varargin)
% Input Arguments
% [P,Y] ana_prep(J,U)
% J is the job parameters json which contains the following fields:
% {
% "jobName":"<Name of job>",
% "path2project":"<Path to Accelerated-Network-Analysis>",
% "path2drase":"<path to drase>",
% "path2pso":"<path to pso>",
% "path2shapes":"<path to shapes>",
% "drasejobParamsFile":"<path to json file for drase job parameters>",
% "inFileData":"<file containing raw data>",
% "inFilePSD":"/work/08302/tcruz/Accelerated-Network-Analysis/TMPPSDDATA/pwPSD.mat",
% "inFileshpsPSD":"/work/08302/tcruz/Accelerated-Network-Analysis/TMPPSDDATA/estPSD.mat",
% "outDir":"/scratch/08302/tcruz/Accelerated-Network-Analysis",
% "scrtchDir":"/scratch/08302/tcruz/Accelerated-Network-Analysis",
% "psoParamsjson":"pso.json",
% "signalParamsjson":"signal.json",
% "injSig":"",
% "accessctrl":"sep",
% "jbTime": 1.5,
% "nNodes": 1,
% "nJbsPerNode": 2,
% "qType":"normal",
% "email":"thomas.cruz03@utrgv.edu"
% }
% U is the user-specified identification number used for file keeping.
% Optional Input Arguments 
% [P,Y] = ana_basics(J,U,D,E,C)
% Setting D allows the user to set a specific date for folder creation
% This could be used to pull data from a specific date if E is not empty,
% or overwride pre-existing data if it is. E, if specified, returns the
% function after only writing file strings for P and Y.
% C controls whether or not a progress.txt file is generated for the run, 1 is yes, else is no.
%
% Optional Output Arguments
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
for largs = 1:(nargin-nreqArgin)
    if ~isempty(varargin{largs})
        switch largs
            case 1
                datad = varargin{largs};
            case 2
                anabreak = varargin{largs};
            case 3
                progCtrl = varargin{largs};
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
