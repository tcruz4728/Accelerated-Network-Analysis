function files = ana_basics_multi(params_jb,filepaths,opts)
% This function uses the job parameters defined in J and the UID U to
% produce the parameter files P and defines the output file prefix which
% files will be saved with. The program loads in time series data with
% load_mtchdfltrdata and saves a PSD training segment to be estimated by
% SHAPES and then conditions the PSD and SHAPES estimated PSD for use in
% matched filtering analysis with rungwpso.m.
%Inputs
% J- Job's parameter json file path or MATLAB structure which contains
% pathing for data and JSON files. Required fields are defined below:
%   {
% "jobName":"<a name for this job which will become the prefix for file
%            names>",
% "path2project":"<path to the project directory>",
% "path2drase":"<path to the DRASE directory>",
% "path2pso":"<path to the PSO directory>",
% "path2shapes":"<path to the SHAPES directory>",
% "inFile":"<path to the file containing pwelch PSD training data>",
% "outFile":"<path to the file containing SHAPES estimated PSD
% training data>",
% "outFilePSD":"<path to the file containing pwelch PSD training data to be
% run by rungwpso>",
% "outFileshpsPSD":"<path to the file containing SHAPES estimated PSD training
% data to be run by rungwpso>",
% "inFileData":"<path to the file containing detector strain data>",
% "outDir":"<path to directory where all output files will be stored, a
% subdirectory for the date is created under this directory. A subdirectory
% under outDir/<date> with a UID>",
% "scrtchDir": <path to directory where all temporary files, e.g., SLURM
%               file, will be stored>",
% (REQUIRED FOR RUNNING SHAPES)
% "genDataParamsjson":"<path to JSON file containing parameters for
%                generating line data>",
% "draseParamsFile":"<path to JSON file containing parameters for each
%                SHAPES run using drase>",
% }
% U- Specifies the user identification number for file indexing beyond
% date. See DRASE\mkdirUID.m for more information. 
%Outputs
% P- Parameter file containing necessary fields for rungwpso. Contains
% structure 'params' defined in gwpsoparams.m.
% Y- Output File Naming Prefix defined by J.outDir,J.jobName, and output
% from gwpsoparams.m
%Optional input Arguments 
% [P,Y] = ana_basics(J,U,D,E,C,T)
% D- Allows the user to specify date for folder indexing.
%   If =<non-existing date>, creates new folder for specified date.
%   If =<pre-existing date>, overwrites pre-existing data in folder.
% E- Returns the function earlier after completing certain functions.     
%   If =1, returns after only writing file strings for P and Y,
%   If =2, returns after =1's response and creating parameter files.
% C- Controls whether or not a progress.txt file is generated for the run.
%   If =1, progress file is generated.
%   Else, (DEFAULT) no progress file is generated.
% T- Plotting control option,
%   If =1, Glitch spectrogram of the training series is plotted
%   Elseif =2, Pwelch PSD training segment is plotted
%   Elseif =12, Both 1's and 2's are plotted,
%   Else, (DEFAULT) No plots are generated. 
%
%Optional output Arguments
% [P,Y,F,N] = ana_basics(J,U,D,E,C,T)
% F- Returns a structure of file paths directories for figures, PSDs, and
% end files as defined in DRASE\dpfc.m 
% N- Returns the full filename of the progress.txt file (if generated). 
%
% See also GWPSOPARAMS, FILETAGANA, DRASE\MKDIRUID, DRASE\DRASE.

%%
arguments
    params_jb           {mustBeNonempty}
    filepaths               (1,1) struct = struct()
    opts.RunID              (1,1) string  = ""
    opts.EarlyReturn        (2,1) logical = [false,false]
    opts.ProgressMonitoring (1,1) logical = false
    opts.IncludePlotting    (2,1) logical = [false,false]
end
%% jobParams datatype check - checks for whether jobParams is a structure, else assumes a file
if ~isstruct(params_jb)
    params_jb = loadjson(params_jb);
end

%% Initial Setup: Parameters - job settings
paramsFile = [filepaths.intermediate,'params']; %rungwpso params file
paramsFileshps = [paramsFile,'shps']; %rungwpso params file for shapes data

%Project Parameters
psoParams = loadjson(params_jb.configs.psoParamsfile); %matched filtering PSO params
signalParams = loadjson(params_jb.configs.signalParamsfile); %signal injection params

%% GW Parameter settings - combines relevant settings and performs necessary 
% computations for rungwpso. Does not require time-series or PSD data
        
gwpsoparams(psoParams,signalParams,paramsFile);
copyfile([paramsFile,'.mat'],[paramsFileshps,'.mat']); %identical but separate parameter settings

switch params_jb.inFileFormat
    case "singleton"
        fileName = params_jb.inFile;
        runs = 1;
    case "manifest"
        fileName = loadjson(params_jb.inFile,"fileNames");
        runs = length(fileName);
        fileName()
end
outFilePSDs = string(zeros(runs,1));
for j = 1:runs
    outFilePSDs(j) = fullfile(filepaths.intermediate,sprintf('PSD_%03.0f',j));
    %% Data Load - loads time-series performs bandpass, computes training segment PSD
    % input: inFileData - time series data from LIGO or simulations
    % output: inFile - training segment PSD
    outData = load_mfdata(fileName(j),outFilePSDs(j),...
        jobParams.injSig);

    %% SHAPES PSD estimate - takes pwelch linear PSD and returns in same form
    % input: inFile - training segment PSD from load_mfdata.m
    % output: outFile - shapes estimation of training segment PSD
    results = drase4lines(params_jb,filepaths);
    %% Interpolation - Takes log10 of PSDs, interpolates and inverses the log
    %input: outData/inFile - data structure from load_mfdata.m
    %output: inFile - appending the interpolated PSD

    createPSD(outData.PSD,outData.freqVec,outData.tlen,outData.sampFreq,params_jb.inFile);
    %input: outFile - estimated PSD from drase
    %output: outFile - appending the interpolated shpsPSD
    createPSD(results.estimate,outData.freqVec,outData.tlen,outData.sampFreq,params_jb.outFile);

    %% Condition Data and Compute FFTs
    %inputs: inFile - interpolated PSD highpassed time series
    %output: paramsFile - updated from gwpsoparams with fft
    cond_mfdata(params_jb.inFile,paramsFile);
    %inputs: outFile - interpolated shps PSD & inFile - highpassed time series
    %output: paramsFileshps - updated from gwpsoparams with fft
    cond_mfdata(params_jb.outFile,paramsFileshps);
end

%Display file names
disp(['ana_basics- Parameter files saved: ',paramsFile, '.mat and ',paramsFileshps,'.mat'])
end