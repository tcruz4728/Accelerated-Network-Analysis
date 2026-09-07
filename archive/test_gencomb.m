
path2jsonlab = 'C:\Users\tcruz\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\JSONLab_ a toolbox to encode_decode JSON files\jsonlab-2.0';
jobParamsFile = 'JSON\comb_PC_Job_params.json';
inFileList = '/Users/tcruz/OneDrive/Onedrive_Documents/GitHub/Accelerated-Network-Analysis/SCRATCH\15-Jan-2026\00001\gwtsnr_premf_PC_Job_tau_fs4096stp5tsL512ta138snr30_inFilesList.txt';
outFileList = '/Users/tcruz/OneDrive/Onedrive_Documents/GitHub/Accelerated-Network-Analysis/SCRATCH\15-Jan-2026\00001\gwtsnr_mf_PC_Job_tau_fs4096stp5tsL512ta138snr30_outFilesList.txt';

gencombmultilnchrjb(path2jsonlab,jobParamsFile,inFileList,outFileList)