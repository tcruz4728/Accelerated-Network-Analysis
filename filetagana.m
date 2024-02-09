function filetagstr = filetagana(psoParams,signalParams)
% T = filetag(P,S)
%Filetag creator for ana-related functions using specific parameters from P
%and S. Order of functions should be filetagana->createPSD->rungwpso.
%Typical naming procedure will have outFilePrfx = [outFilePath,filetagstr].
% P is a JSON file with parameters relating to pso's use, S is a JSON file
% with parameters relating to the injected signal.
%
%Created May. 2023 by Thomas Cruz from DRASE/filetag.m
switch psoParams.type
    case 1
        filetagstr = 'mass_';
    case 2
        filetagstr = 'tau_';
end
% filetagstr = [psoParams.type,'_'];

for  lp=1:5
    switch lp
        case 1
            filetagstr = [filetagstr,'fs',num2str(signalParams.sampling_freq)];
        case 2
            filetagstr = [filetagstr,'stp',num2str(psoParams.maxSteps)];
        case 3
            filetagstr = [filetagstr,'tsL',num2str(signalParams.signal.T_sig_len)];
        case 4
            filetagstr = [filetagstr, 'ta',num2str(signalParams.ta)];
        case 5
            filetagstr = [filetagstr, 'snr',num2str(signalParams.snr)];
    end
end

end