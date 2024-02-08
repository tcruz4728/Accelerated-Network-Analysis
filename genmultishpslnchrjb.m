function [outFileNameList] = genmultishpslnchrjb(path2jsonlab,jobParamsFile,varargin)
% path2jsonlab = 'C:\Users\tcruz\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\JSONLab_ a toolbox to encode_decode JSON files\jsonlab-2.0';
% jobParamsFile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\JSON\multi_mtchdfltr_PC_Job_params.json';
addpath(path2jsonlab)

% userUID = input("UID:");
jobParams = loadjson(jobParamsFile);

%File Naming Convention
[interFilePath,interFileName,~] = fileparts(jobParams.inFilePSD);

inFileNameList = cell(jobParams.inFileDataRange(2),1);
interFileNameList = cell(jobParams.inFileDataRange(2),1);
outFileNameList = cell(jobParams.inFileDataRange(2),1);
for fileCount = 1:jobParams.inFileDataRange(2)
    inFileNameList{fileCount} = [jobParams.inFileDataPrFx,...
        num2str(fileCount),'.mat'];
    interFileNameList{fileCount} = [interFilePath,filesep,...
        interFileName,'_n',num2str(fileCount),'.mat'];
    outFileNameList{fileCount} = [jobParams.inFileshpsPSD,'_n',...
        num2str(fileCount),'.mat'];
end
[~,outdataFilePrfx] = ana_prep(jobParams,varargin{1},[],1);

%% Job Text File
fidJbFile = fopen([jobParams.scrtchDir,filesep,jobParams.jobName,'_jbfile.txt'],'w');
disp(['Job File created: ',jobParams.scrtchDir,filesep,jobParams.jobName,'_jbfile.txt'])
%Store list of output files in .txt file for post-processing codes
fidOutFileList = fopen([outdataFilePrfx,'shpsoutFilesList.txt'],'w');
disp(['Output File list created: ',outdataFilePrfx,'shpsoutFilesList.txt'])

nJobs = 1;
for nCount = 1:fileCount
    fprintf(fidJbFile,'matlab -batch ');
    %path to jsonlab,
    fprintf(fidJbFile,' "addpath ''%s''; ', path2jsonlab);
    %path to DRASE
    fprintf(fidJbFile,' addpath ''%s''; ', jobParams.path2drase);
    %path to SHAPES, PSO, and project
    fprintf(fidJbFile,' setpath(''%s''); ', jobParamsFile);
    fprintf(fidJbFile,' load_mtchdfltrdata(''%s'',''%s''); ',...
        inFileNameList{nCount},interFileNameList{nCount});
    %SHAPES call 
    fprintf(fidJbFile,' run_drase4lines(''%s'',''%s'',''%s'',''%s'');" \n',...
        jobParamsFile,num2str(varargin{1}),interFileNameList{nCount},outFileNameList{nCount}); %%%%%%%%%%FIX USERUID%%%%%%%%
    fprintf(fidOutFileList,'%s\n',outFileNameList{nCount});
    nJobs = nJobs +1;
end
fclose(fidJbFile);
fclose(fidOutFileList);
%% Slurm File Generation
genslurm(jobParams,nJobs)