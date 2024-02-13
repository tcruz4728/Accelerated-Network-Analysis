function [outFileNameList] = genmultishpslnchrjb(path2jsonlab,jobParamsFile,varargin)
% path2jsonlab = 'C:\Users\tcruz\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\JSONLab_ a toolbox to encode_decode JSON files\jsonlab-2.0';
% jobParamsFile = 'C:\Users\tcruz\OneDrive\Onedrive_Documents\GitHub\Accelerated-Network-Analysis\JSON\multi_mtchdfltr_PC_Job_params.json';
addpath(path2jsonlab)

% userUID = input("UID:");
jobParams = loadjson(jobParamsFile);

%File Naming Convention
[paramsFile,outdataFilePrfx] = ana_prep(jobParams,varargin{1});
paramsFileshps = [paramsFile,'shps'];
[interFilePath,interFileName,~] = fileparts(jobParams.inFilePSD);
[outFilePath,outFileName,~] = fileparts(jobParams.inFileshpsPSD);
inFileNameList = cell(jobParams.inFileDataRange(2),1);
interFileNameList = cell(jobParams.inFileDataRange(2),1);
outFileNameList = cell(jobParams.inFileDataRange(2),1);
paramsFileList = cell(jobParams.inFileDataRange(2),1);
paramsFileshpsList = cell(jobParams.inFileDataRange(2),1);
for fileCount = 1:jobParams.inFileDataRange(2)
    %Time series initial file
    inFileNameList{fileCount} = [jobParams.inFileDataPrFx,...
        num2str(fileCount),'.mat'];
    %Pwelch data file (inFilePSD-goes into SHAPES)
    interFileNameList{fileCount} = [interFilePath,filesep,...
        interFileName,'_n',num2str(fileCount),'.mat'];
    %SHAPES estimated PSD data file (inFileshpsPSD)
    outFileNameList{fileCount} = [outFilePath,filesep,...
        outFileName,'_n',num2str(fileCount),'.mat'];
    %Pwelch params files
    paramsFileList{fileCount} = [paramsFile,'_n',...
        num2str(fileCount),'.mat'];
    %SHAPES params files
    paramsFileshpsList{fileCount} = [paramsFileshps,'_n',...
        num2str(fileCount),'.mat'];
    %Copying root file to others, calculated params do not change
    copyfile([paramsFile,'.mat'],paramsFileList{fileCount});
    copyfile([paramsFileshps,'.mat'],paramsFileshpsList{fileCount});
end

%% Job Text File
fidJbFile = fopen([jobParams.scrtchDir,filesep,jobParams.jobName,'_jbfile.txt'],'w');
disp(['Job File created: ',jobParams.scrtchDir,filesep,jobParams.jobName,'_jbfile.txt'])
%Store list of output files in .txt file for post-processing codes
fidOutFileList = fopen([outdataFilePrfx,'_inFilesList.txt'],'w');
disp(['Input File list created: ',outdataFilePrfx,'_inFilesList.txt'])
fidparamsFileList = fopen([outdataFilePrfx,'_paramsFilesList.txt'],'w');
disp(['Parameter File list created: ',outdataFilePrfx,'_paramsFilesList.txt'])

nJobs = 1;
for nCount = 1:fileCount
    fprintf(fidJbFile,'matlab -batch ');
    %path to jsonlab,
    fprintf(fidJbFile,' "addpath ''%s''; ', path2jsonlab);
    %path to DRASE
    fprintf(fidJbFile,' addpath ''%s''; ', jobParams.path2drase);
    %path to SHAPES, PSO, and project
    fprintf(fidJbFile,' setpath(''%s''); ', jobParamsFile);
    %Data Load
    fprintf(fidJbFile,' load_mtchdfltrdata(''%s'',''%s''); ',...
        inFileNameList{nCount},interFileNameList{nCount});
    %SHAPES call 
    fprintf(fidJbFile,' run_drase4lines(''%s'',''%s'',''%s'',''%s''); ',...
        jobParamsFile,outdataFilePrfx,interFileNameList{nCount},outFileNameList{nCount});
    %PSD Interpolation
    fprintf(fidJbFile,' createPSD(''%s''); ',... welch
        interFileNameList{nCount});
    fprintf(fidJbFile,' createPSD(''%s''); ',... shapes
        outFileNameList{nCount});
    %Matched filtering conditioning
    fprintf(fidJbFile,' cond_mtchdfltrdata(''%s'',''%s'',''%s'');', ...
        interFileNameList{nCount},paramsFileList{nCount},jobParams.injSig);
    %SHAPES conditioning
    fprintf(fidJbFile,' cond_mtchdfltrdata(''%s'',''%s'',''%s'');" \n', ...
        outFileNameList{nCount},paramsFileshpsList{nCount},jobParams.injSig);
    fprintf(fidOutFileList,'%s  %s  %s\n',...
        inFileNameList{nCount},interFileNameList{nCount},...
        outFileNameList{nCount});
    fprintf(fidparamsFileList,'%s  %s\n',...
        paramsFileList{nCount},paramsFileshpsList{nCount});
    nJobs = nJobs +1;
end
fclose(fidJbFile);
fclose(fidOutFileList);
%% Slurm File Generation
genslurm(jobParams,nJobs)