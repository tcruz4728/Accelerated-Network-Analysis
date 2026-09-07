function results = runAnaBatch(config,opts)

arguments
    config                      (1,1) struct = struct()
    opts.ParameterOverride      (1,1) struct = struct()
    opts.UseDelayedSegments    (1,1) logical = false
    opts.nRuns                  (1,1) double = 6
    opts.RunID                 (1,1) string  = ""
end


jobTemplatePath = config.job_config;

if ~isfile(jobTemplatePath)
    error("runAnaBatch:JobTemplateNotFound", ...
        "AANA job config does not exist:\n%s", jobTemplatePath);
end

jobParams = jsondecode(fileread(jobTemplatePath));

validateJobTemplate(jobParams);

% Preserve the input selected by your existing job template initially.
% Later, move this selection into config.inputs if desired.
jobParams.jobName = char(config.batch.id + "_" + config.run.id);

% Redirect all generated output to this exact run.
jobParams.outDir = char(config.run.tables_dir);

if strlength(config.run.intermediate_dir) == 0
    error("runAnaBatch:IntermediateDisabled", ...
        "This batch was commissioned without an intermediate directory.");
end

jobParams.scrtchDir = char(config.run.intermediate_dir);

% Save an immutable record of the exact DRASE/ANA configuration used.
jobConfigPath = fullfile( ...
    config.run.metadata_dir, ...
    "ana_job_" + config.run.id + ".json");

writeJsonFile(jobConfigPath, jobParams);

fprintf("ANA/DRASE job config:   %s\n", jobConfigPath);
fprintf("Input file:             %s\n", string(jobParams.inFile));
fprintf("Output directory:       %s\n", string(jobParams.outDir));
fprintf("Scratch directory:      %s\n", string(jobParams.scrtchDir));

%% SNEAR
%Setup

% Addpaths to dependencies
addpath(fullfile(config.paths.environment.drase_dir,'src/','utils/'))
setpath(jobParams.environment)

% Prepare the file paths for the SNEAR processing
% jobParams.path2projectdata = char(fullfile(config.run.tables_dir, 'project_data'));
filepaths.extracted = char(fullfile(config.run.extracted_dir));
filepaths.intermediate = char(fullfile(config.run.intermediate_dir));
filepaths.estimated = char(fullfile(config.run.estimated_dir));
filepaths.tables = char(fullfile(config.run.tables_dir));
filepaths.figures = char(fullfile(config.run.figures_dir));

% Snear Call
pipelineOutput = ana(jobParams,filepaths,...
    "ParameterOverride",opts.ParameterOverride,...
    "UseDelayedSegments",opts.UseDelayedSegments,...
    "nRuns",opts.nRuns,...
    "RunID",opts.RunID,...
    "ProgressMonitoring",true);

%% Results
result = struct();
% result.status = "completed";
result.job_config_path = string(jobConfigPath);
result.input_file = string(jobParams.inFile);
result.output_dir = string(jobParams.outDir);
result.scratch_dir = string(jobParams.scrtchDir);
result.pipeline_output = pipelineOutput;

end

function validateJobTemplate(job_params)
% Check required fields in the job parameters
requiredFields = {'jobName','dataRoot','environment','inFile','inFileFormat','configs'};
missingFields = setdiff(requiredFields, fieldnames(job_params));

if ~isempty(missingFields)
    error("validateJobTemplate:MissingFields", ...
        "The following required fields are missing: %s", strjoin(missingFields, ', '));
end

% Verify the existence of the data root directory
if ~isfolder(job_params.dataRoot)
    error("validateJobTemplate:DataRootNotFound", ...
        "The specified data root directory does not exist: %s", job_params.dataRoot);
end
% Ensure the input folder exists
if ~isfile(job_params.inFile)
    error("validateJobTemplate:InputFileNotFound", ...
        "The specified input file does not exist: %s", job_params.inFile);
end
% Verify inFile is a valid filetype, either an .hdf5 or .mat file.
[~,~,ext] = fileparts(job_params.inFile);
if ~ismember(ext, {'.hdf5', '.mat','.json'})
    error("validateJobTemplate:InvalidFileType", ...
        "The input file must be either an .hdf5, .mat, or .json file: %s", job_params.inFile);
end

% Verify configs has fields draseParamsfile and genDataParamsfile
requiredConfigFields = {'draseParamsfile', 'genDataParamsfile',...
                        'psoParamsfile','signalParamsfile'};
missingConfigFields = setdiff(requiredConfigFields, fieldnames(job_params.configs));

if ~isempty(missingConfigFields)
    disp(job_params.configs)
    error("validateJobTemplate:MissingConfigFields", ...
        "The following required config fields are missing: %s", strjoin(missingConfigFields, ', '));
end

end