%test_ana_script
ana_parameter_override = struct('variable','nKnts',...
        'type','single');
batch_description = 'testing batch creation with ANA';
batch_focus = "test";
batch_ID = 1;
confirm_check = false;
use_delayed_segments = false;
%% Do not change below
%Load environment config
projectPaths = loadjson('configs/environment.local.json');
addpath(genpath(projectPaths.path2drase))
setpath(projectPaths)

overrides = struct();

overrides.processing = struct();
overrides.processing.algorithm = "ana";

environment = struct();
environment.project_dir = projectPaths.path2project;
environment.drase_dir = projectPaths.path2drase;
environment.pso_dir = projectPaths.path2pso;
environment.shapes_dir = projectPaths.path2shapes;

path2templateFile = fullfile(projectPaths.path2project,'configs','templates','ana.json');
configPath = commissionBatch(...
    "DataRoot",projectPaths.path2projectdata,...
    "Focus",batch_focus,"TemplateFile",path2templateFile, ...
    "Interactive",false,"Overrides",overrides,"Confirm",confirm_check, ...
    "AllowExisting",true,"Environment",environment, ...
    "Description",batch_description,"BatchID",batch_ID);

switch ana_parameter_override.variable
    case 'nKnts'
        switch ana_parameter_override.type
            case 'large_sweep'
                ana_parameter_override.value = [5,7,10,15,17,20,25,30,35,40];
            case 'focused'
                ana_parameter_override.value = [15,16,17,18,19,20];
            case 'single'
                ana_parameter_override.value = 17;
        end
        ana_parameter_override.control = 20;
    case 'rGain'
        switch ana_parameter_override.type
            case 'large_sweep'
                ana_parameter_override.value = [0.55,0.45,0.35,0.25,0.15,0.05,0.025,0.005];
            case 'focused'
                ana_parameter_override.value = [0.075 0.065 0.05 0.04 0.03];
            case 'single'
                ana_parameter_override.value = 0.05;
        end
        ana_parameter_override.control = 0.0;
end

path2jsonlab = 'C:\Users\tcruz\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\JSONLab_ a toolbox to encode_decode JSON files\jsonlab-2.0';
addpath(path2jsonlab)

runInfo = runBatch("ConfigPath",configPath,"Processor",@runAnaBatch,...
    "ConfirmRun",confirm_check,"ParameterOverride",ana_parameter_override,...
    "UseDelayedSegments",use_delayed_segments);

if runInfo.status ~= "cancelled_before_run"
    fidprog = fopen(progressFile,'a');
    proglines = struct('nd','done.',...
        'p','Post-processing...');

    progstatus(proglines.p,fidprog,ProgressMonitoring)
    combFileName = comb_anashpsjb(path2jsonlab,path2templateFile);
    postprocessing(combFileName,filepaths,path2templateFile.injSig);
    progstatus(proglines.nd,fidprog,ProgressMonitoring)
    fclose(fidprog);
end
% ana_post_paramsfile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\SNEAR\configs\templates\snear_post.yml';
% results = snear_post(runInfo.processor_result.pipeline_output,ana_post_paramsfile,ana_parameter_override,"SaveResults",{true,runInfo.tables_dir});
% 
% snear_plots(results,"ParameterOverride",ana_parameter_override,"SaveFigures",{false,runInfo.figures_dir},"SaveOldPlots",true);
