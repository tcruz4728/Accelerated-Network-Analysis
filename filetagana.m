function filetagstr = filetagana(psoParams,signalParams)
% T = filetag(P,S)
%Filetag creator for ana-related functions using specific parameters from P
%and S. Order of functions should be filetagana->createPSD->rungwpso.
%Typical naming procedure will have outFilePrfx = [outFilePath,filetagstr].
% P is a JSON file with parameters relating to pso's use, S is a JSON file
% with parameters relating to the injected signal.
%
%Created May. 2023 by Thomas Cruz from DRASE/filetag.m
filetagstr = '';

% Prefix from psoParams.type
if isfield(psoParams,'type') && ~isempty(psoParams.type)
    switch psoParams.type
        case 1
            filetagstr = 'mass_';
        case 2
            filetagstr = 'tau_';
        otherwise
            % optional: default prefix or leave empty
    end
end

% fs (sampling_freq)
if isfield(signalParams,'sampling_freq') && ~isempty(signalParams.sampling_freq)
    filetagstr = [filetagstr, 'fs', num2str(signalParams.sampling_freq)];
end

% stp (maxSteps)
if isfield(psoParams,'maxSteps') && ~isempty(psoParams.maxSteps)
    filetagstr = [filetagstr, 'stp', num2str(psoParams.maxSteps)];
end

% tsL (T_sig_len)
if isfield(signalParams,'signal') && isfield(signalParams.signal,'T_sig_len') && ...
        ~isempty(signalParams.signal.T_sig_len)
    filetagstr = [filetagstr, 'tsL', num2str(signalParams.signal.T_sig_len)];
end

% ta
if isfield(signalParams,'ta') && ~isempty(signalParams.ta)
    filetagstr = [filetagstr, 'ta', num2str(signalParams.ta)];
end

% snr
if isfield(signalParams,'snr') && ~isempty(signalParams.snr)
    filetagstr = [filetagstr, 'snr', num2str(signalParams.snr)];
end


end