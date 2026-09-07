function [] = genmultimflnchrjb(path2jsonlab,jobParamsFile,varargin)
%GENMULTIMFLNCHRJB(J,P)
% Generates a .slurm file containing a LAUNCHER job for matched filtering
% of both pwelch PSDs and (optional) SHAPES estimated PSDs. It is necessary
% to 1st run genmultipremflnchrjb.m before this function to have the proper
% parameter files created and to condition time-series data for MF.
% J is the path to the jsonlab package. If set to '', the jsonlab package
% is assumed to be in the Matlab search path. Note that this search path
% must be accessible from every compute node, which cannot always be
% guaranteed. Hence it is safer to specify the path explicitly.
%
% P is a JSON file containing the following job parameters. Text in <>
% should be replaced by an appropriate value. Examples are shown for some
% of the parameters.
%   {
% "jobName":"<a name for this job which will become the prefix for file
%            names>",
% "path2project":"<path to the project directory>",
% "path2drase":"<path to the DRASE directory>",
% "path2pso":"<path to the PSO directory>",
% "path2shapes":"<path to the SHAPES directory>",
% "inFileDataPrfx":"<file prefix and location of data realizaitons>",
% "inFileDataRange":"<1x2 array giving the range of values from which the
%       full file names of data realizations, i.e. [1, 5] will have 5 data
%       realizations loaded in and estimated.>",
% "inFile":"<path to the file containing the pwelch estimated PSD
%       training data>",
% "outFile":"<path to the file containing SHAPES estimated PSD
%       training data>",
% "outDir":"<path to directory where all output files will be stored, a
%       subdirectory for the current date is created under this directory. 
%       Another subdirectory under outDir/<date> with a UID is also set.>",
% "scrtchDir": <path to directory where all temporary files, e.g., SLURM
%       files, job text files, will be stored>",
% "psoParamsjson":"<path to JSON file containing parameters for rungwpso's
%       PSO run.>",
% "signalParamsjson":"<path to JSON file containing parameters for a
%       realized signal or an injected signal>",
% "injSig":"<Injection signal control parameter, if empty no signal
%       injection is performed, else the injected signal will have 
%       parameters set by signalParams.json>",
% "jbTime":"<Time per Matlab job in hours>",
% "nNodes":"<Number of compute nodes for this job>",
% "qType":"<job queue: leave empty for default/ls6,skx-normal for 
%       stampede2>",
% "email":"<email address>"
%   }
% }
% Optional Input Arguments
%GENMULTIMFLNCHRJB(J,P,U,D,S)
% U is an unique ID number that is used when creating a folder under the
% outfile directory. If specified, it uses that value with the format %05d,
% i.e. U = 1, out directory is <outDir>/00001. Can be used to overwrite
% preexisiting data(provided the date is the same) if an already used UID
% is specified.
% D is a date specifying variable which if not set, uses the current date
% for folder indexing. Date formats are written as '01-Jan-2000'
% S is the SHAPES control parameter which determines whether matched
% filtering Jobs will be created for a SHAPES estimated PSD in addition to
% the norm. S = 1 is the default which creates the SHAPES jobs, set to 0
% to disable SHAPES jobs.

% Modified from genlineshpslnchrjb, Aug 2023 for ANA use
% Modified name from genmultipsolnchrjb.m to better match function purpose,
% Jul 2025

%Add path to jsonlab
addpath(path2jsonlab);

%Load job parameters
jobParams = loadjson(jobParamsFile);

%Optionals for ana_basics setup
userUID = 1;
datad = [];
shpsCtrl = 1;
nBatches = 1; % number of matlab -batch lines per launcher job file

%Override the file name if optional input given
nreqArgs = 2;
for lpargs = 1:(nargin-nreqArgs)
    if ~isempty(varargin{lpargs})
        switch lpargs
            case 1
                userUID = varargin{lpargs};
            case 2
                datad = varargin{lpargs};
            case 3
                shpsCtrl = varargin{lpargs};
            case 4
                nBatches = varargin{lpargs};
        end
    end
end

% progCtrl = 1; % Always set to on for debugging

%% ANA prep package
[paramsFile,outdataFilePrfx] = ana_basics(jobParams,userUID,datad,1);

paramsFileList = cell(jobParams.inFileDataRange(2),1);
dataFileList = cell(jobParams.inFileDataRange(2),1);

if shpsCtrl == 1
    paramsFileshps = [paramsFile,'shps'];
    paramsFileshpsList = cell(jobParams.inFileDataRange(2),1);
    shpsDataFileList = cell(jobParams.inFileDataRange(2),1);
end

for fileCount = jobParams.inFileDataRange(1):jobParams.inFileDataRange(2)
    paramsFileList{fileCount} = [paramsFile,'_n',...
        num2str(fileCount),'.mat'];
    dataFileList{fileCount} = [outdataFilePrfx,'_n',...
        num2str(fileCount),'C'];
    if shpsCtrl == 1
        paramsFileshpsList{fileCount} = [paramsFileshps,'_n',...
            num2str(fileCount),'.mat'];
        shpsDataFileList{fileCount} = [outdataFilePrfx,'_n',...
            num2str(fileCount),'shps_C'];
    end
end
% fidShpsOutFileList = fopen([outdataFilePrfx,'shpsoutFilesList.txt'],'r');
%% Constuct job file for Launcher
% fidJbFile = fopen([jobParams.scrtchDir,filesep,jobParams.jobName,'_jbfile.txt'],'w');
% disp(['Job File: ',jobParams.jobName,'_jbfile.txt',' file created in ',jobParams.scrtchDir,filesep])
% %Store list of output files in .txt file for post-processing codes
% fidOutFileList = fopen([outdataFilePrfx,'_outFilesList.txt'],'w');
% disp(['Output File list created: ',outdataFilePrfx,'_outFilesList.txt'])
% 
% nJobs = 1;
% batchIndex = 1;            % which batch we are on
% jobCountInBatch = 0;       % how many jobs currently in this batch
% nJobsTotal = 0;            % overall job counter
% 
% % Helper to open a new batch job file and out‑list
% openNewBatch = @(bIdx) deal( ...
%     fopen(fullfile(jobParams.scrtchDir, ...
%            sprintf('%s_batch%03d_jbfile.txt', jobParams.jobName, bIdx)), 'w'), ...
%     fopen(sprintf('%s_batch%03d_outFilesList.txt', outdataFilePrfx, bIdx), 'w') );
% 
% [fidJbFile, fidOutFileList] = openNewBatch(batchIndex);
% fprintf('Job File: %s_batch%03d_jbfile.txt created in %s', ...
%     jobParams.jobName, batchIndex, jobParams.scrtchDir);
% fprintf('Output File list created: %s_batch%03d_outFilesList.txt\n', ...
%     outdataFilePrfx, batchIndex);
% 
% for nCount = jobParams.inFileDataRange(1):jobParams.inFileDataRange(2)
%     for runType = 1:(1+shpsCtrl)
%         % If this batch is full, close files and start a new batch
%         if jobCountInBatch >= batchSize
%             fclose(fidJbFile);
%             fclose(fidOutFileList);
%             batchIndex       = batchIndex + 1;
%             jobCountInBatch  = 0;
%             [fidJbFile, fidOutFileList] = openNewBatch(batchIndex);
%             fprintf('Job File: %s_batch%03d_jbfile.txt created in %s', ...
%                 jobParams.jobName, batchIndex, jobParams.scrtchDir);
%             fprintf('Output File list created: %s_batch%03d_outFilesList.txt\n', ...
%                 outdataFilePrfx, batchIndex);
%         end
% 
%         %  PSO and drase command on input file.
%         fprintf(fidJbFile,'matlab -batch ');
%         %path to jsonlab,
%         fprintf(fidJbFile,' "addpath ''%s''; ', path2jsonlab);
%         %path to DRASE
%         fprintf(fidJbFile,' addpath ''%s''; ', jobParams.path2drase);
%         %path to SHAPES, PSO, and project
%         fprintf(fidJbFile,' setpath(''%s''); ', jobParamsFile);
%         %Call to matched filtering code
%         switch runType
%             case 1 %pwelch run
%                 fprintf(fidJbFile, ' rungwpso(''%s'',''%s'');" \n', ...
%                     paramsFileList{nCount},dataFileList{nCount});
%             case 2 %shapes run
%                 fprintf(fidJbFile, ' rungwpso(''%s'',''%s'');" \n', ...
%                     paramsFileshpsList{nCount},shpsDataFileList{nCount});           
%         end
%         jobCountInBatch = jobCountInBatch + 1;
%         nJobsTotal      = nJobsTotal + 1;
% 
%         % Write to this batch's out‑file list
%         fprintf(fidOutFileList,'%s', dataFileList{nCount});
%         if runType == 2
%             fprintf(fidOutFileList,'  %s\n', shpsDataFileList{nCount});
%         else
%             fprintf(fidOutFileList,'\n');
%         end
%     end
% end
% fclose(fidOutFileList);
% fclose(fidJbFile);
% %% Slurm file generation
% genslurm(jobParams,nJobs)
% end

%% Job Count
nJobsTotal = 0;
for nCount = jobParams.inFileDataRange(1):jobParams.inFileDataRange(2)
    for runType = 1:(1+shpsCtrl)
        nJobsTotal = nJobsTotal + 1;
    end
end
jobsPerBatch    = ceil(nJobsTotal / nBatches);

%% Second pass: construct batched job files and output lists
batchIndex      = 1;
jobCountInBatch = 0;

batchJobCounts  = [];    % jobs in each batch
batchNames      = {};    % jobName for each batch

openNewBatch = @(bIdx, baseJobName) deal( ...
    fopen(fullfile(jobParams.scrtchDir, ...
          sprintf('%s_batch%03d_jbfile.txt', baseJobName, bIdx)), 'w'), ...
    fopen(sprintf('%s_batch%03d_outFilesList.txt', outdataFilePrfx, bIdx), 'w'), ...
    sprintf('%s_batch%03d', baseJobName, bIdx) );

[fidJbFile, fidOutFileList, currentBatchName] = ...
    openNewBatch(batchIndex, jobParams.jobName);

disp(['Job File: ', currentBatchName, '_jbfile.txt file created in ', ...
      jobParams.scrtchDir, filesep]);
disp(['Output File list created: ', outdataFilePrfx, ...
      sprintf('_batch%03d_outFilesList.txt', batchIndex)]);

% reset total job counter if you still want it
nJobsTotalCheck = 0;

for nCount = jobParams.inFileDataRange(1):jobParams.inFileDataRange(2)

    for runType = 1:(1+shpsCtrl)

        % If this batch reached its quota, close and start a new batch
        if jobCountInBatch >= jobsPerBatch
            fclose(fidJbFile);
            fclose(fidOutFileList);

            batchJobCounts(end+1) = jobCountInBatch;
            batchNames{end+1}     = currentBatchName;

            batchIndex      = batchIndex + 1;
            jobCountInBatch = 0;

            [fidJbFile, fidOutFileList, currentBatchName] = ...
                openNewBatch(batchIndex, jobParams.jobName);

            disp(['Job File: ', currentBatchName, ...
                  '_jbfile.txt file created in ', ...
                  jobParams.scrtchDir, filesep]);
            disp(['Output File list created: ', outdataFilePrfx, ...
                  sprintf('_batch%03d_outFilesList.txt', batchIndex)]);
        end

        fprintf(fidJbFile,'matlab -batch ');
        fprintf(fidJbFile,' "addpath ''%s''; ', path2jsonlab);
        fprintf(fidJbFile,' addpath ''%s''; ', jobParams.path2drase);
        fprintf(fidJbFile,' setpath(''%s''); ', jobParamsFile);

        switch runType
            case 1 % pwelch run
                fprintf(fidJbFile, ' rungwpso(''%s'',''%s'');" \n', ...
                    paramsFileList{nCount}, dataFileList{nCount});
            case 2 % shapes run
                fprintf(fidJbFile, ' rungwpso(''%s'',''%s'');" \n', ...
                    paramsFileshpsList{nCount}, shpsDataFileList{nCount});
        end

        jobCountInBatch = jobCountInBatch + 1;
        nJobsTotalCheck = nJobsTotalCheck + 1;

        % write to this batch's outFiles list

        if runType == 2
            shpsList = shpsDataFileList{nCount};
            % fprintf(fidOutFileList,'  %s\n', shpsDataFileList{nCount});
        else
            % fprintf(fidOutFileList,'\n');
            shpsList = '';
        end
        fprintf(fidOutFileList,'%s %s\n', dataFileList{nCount},shpsList);
    end
end

% close the final batch
fclose(fidJbFile);
fclose(fidOutFileList);

if jobCountInBatch > 0
    batchJobCounts(end+1) = jobCountInBatch;
    batchNames{end+1}     = currentBatchName;
end

disp(['Total jobs (check): ', num2str(nJobsTotalCheck)]);
disp(['Requested batches:  ', num2str(nBatches)]);
disp(['Actual batches:     ', num2str(numel(batchNames))]);

%% Generate one SLURM file per batch
for b = 1:numel(batchNames)
    thisName  = batchNames{b};
    thisNJobs = batchJobCounts(b);

    % genslurm(jobParams, nJobs, jobNameOverride)
    genslurm(jobParams, thisNJobs, thisName);
end

end
