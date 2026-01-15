function gencombmultilnchrjb(path2jsonlab,jobParamsFile,inFileList,outFileList,varargin)
%Post-matched filtering function to combine data for PC download across a
% cluster
% GENCOMBMULTILNCHRJB(P,F,I,O)
% P- Path to JSON lab dependency
% F- JSON file for creating file name prefixes and slurm parameters.
% I- the input file list used in the 1st ls6 job with SHAPES, contains N
% lines corresponding to the N files with data realizations. Needs to match
% the file created by genmultipremflnchrjb.m
% O- the output file list used in the 2nd ls6 job with rungwpso, contains N
% lines corresponding to the N files with PSD estimates. Needs to match the
% file created by genmultimflnchrjb.m
%
% Modified from combmultimfshps to work on massive data sets when runtimes
% exceed idev session time.
%
% See also genmultipremflnchrjb, genmultimflnchrjb.

%Add path to jsonlab
addpath(path2jsonlab);

%Load job parameters
jobParams = loadjson(jobParamsFile);

%Optionals for ana_basics setup
userUID = 1;
datad = [];

%Override the file name if optional input given
nreqArgs = 2;
for lpargs = 1:(nargin-nreqArgs)
    if ~isempty(varargin{lpargs})
        switch lpargs
            case 1
                userUID = varargin{lpargs};
            case 2
                datad = varargin{lpargs};
        end
    end
end

% progCtrl = 1; % Always set to on for debugging

%% ANA prep package
[~,outdataFilePrfx] = ana_basics(jobParams,userUID,datad,1);

inFilesList  = readcell(inFileList,  "Delimiter","  ");
outFilesList = readcell(outFileList, "Delimiter","  ");
nFiles = min(size(outFilesList,1), size(inFilesList,1));

inFileNameList     = inFilesList(1:nFiles,1);
interFileNameList  = inFilesList(1:nFiles,2);
outFileNameList    = inFilesList(1:nFiles,3);

dataFileList       = outFilesList(1:nFiles,1);
shpsDataFileList   = outFilesList(1:nFiles,2);

%% Constuct job file for Launcher
fidJbFile = fopen([jobParams.scrtchDir,filesep,jobParams.jobName,'_jbfile.txt'],'w');
disp(['Job File: ',jobParams.jobName,'_jbfile.txt',' file created in ',jobParams.scrtchDir,filesep])
%Store list of output files in .txt file for post-processing codes
fidOutFileList = fopen([outdataFilePrfx,'_outFilesList.txt'],'w');
disp(['Output File list created: ',outdataFilePrfx,'_outFilesList.txt'])
nJobs = 1;

for fileCount = 1:nFiles
    % Output MAT file for this realization
    outMat = sprintf('%s_n%dF.mat', outdataFilePrfx, fileCount);
    fprintf(fidOutFileList, '%s\n', outMat);

    in1 = inFileNameList{fileCount};
    in2 = interFileNameList{fileCount};
    in3 = outFileNameList{fileCount};
    in4 = dataFileList{fileCount};
    in5 = shpsDataFileList{fileCount};

    % Escape single quotes in paths
    in1 = strrep(in1,'''','''''');
    in2 = strrep(in2,'''','''''');
    in3 = strrep(in3,'''','''''');
    in4 = strrep(in4,'''','''''');
    in5 = strrep(in5,'''','''''');
    outMatEsc = strrep(outMat,'''','''''');

    % MATLAB command to execute for this realization
    matlabCmd = sprintf([ ...
        "inputData   = load('%s'); " ...
        "psdData     = load('%s'); " ...
        "estpsdData  = load('%s'); " ...
        "outData     = load('%s'); " ...
        "estoutData  = load('%s'); " ...
        "save('%s','inputData','psdData','estpsdData','outData','estoutData');" ...
        "exit"], ...
        in1,in2,in3,in4,in5,outMatEsc);

    % Wrap in matlab -batch "<cmd>"
    % Use single quotes around the whole string in the shell, double quotes inside
    shellLine = sprintf('matlab -batch "%s"\n', matlabCmd);

    fprintf(fidJbFile, '%s', shellLine);
    nJobs = nJobs + 1;
end

fclose(fidJbFile);
fclose(fidOutFileList);

%% Slurm file generation
genslurm(jobParams,nJobs)
end

