function filetagstr = filetagana(psoParams,signalParams)
% T = filetag(P,S)
%Filetag creator for ana-related functions using specific parameters from P
%and S. Order of functions should be filetagana->createPSD->rungwpso.
%Typical naming procedure will have outFilePrfx = [outFilePath,filetagstr].
% P is a JSON file with parameters relating to pso's use, S is a JSON file
% with parameters relating to the injected signal.
%
%Created May. 2023 by Thomas Cruz from DRASE/filetag.m

%% Structure Check
psoParams = ensureStruct(psoParams);
signalParams   = ensureStruct(signalParams);

%% Build filename tag
tags = strings(0, 1);

% Prefix from psoParams.type
% if isfield(psoParams,'type') && ~isempty(psoParams.type)
%     switch psoParams.type
%         case 1
%             filetagstr = 'mass';
%         case 2
%             filetagstr = 'tau';
%     end
% end

% From psoParams
tags(end+1) = getTag(psoParams, 'type',       '',    '%s');
tags(end+1) = getTag(psoParams, 'maxSteps',    'stp',   '%g');

% From signalParams
tags(end+1) = getTag(signalParams, 'sampling_freq', 'fs',  '%g');
tags(end+1) = getTag(signalParams.signal, 'T_sig_len',        'tsL',   '%g');
tags(end+1) = getTag(signalParams, 'ta',        'ta',   '%g');
tags(end+1) = getTag(signalParams, 'snr',        'snr',   '%g');



% fs (sampling_freq)
% if isfield(signalParams,'sampling_freq') && ~isempty(signalParams.sampling_freq)
%     filetagstr = [filetagstr, 'fs', num2str(signalParams.sampling_freq)];
% end

% stp (maxSteps)
% if isfield(psoParams,'maxSteps') && ~isempty(psoParams.maxSteps)
%     filetagstr = [filetagstr, 'stp', num2str(psoParams.maxSteps)];
% end

% % tsL (T_sig_len)
% if isfield(signalParams,'signal') && isfield(signalParams.signal,'T_sig_len') && ...
%         ~isempty(signalParams.signal.T_sig_len)
%     filetagstr = [filetagstr, 'tsL', num2str(signalParams.signal.T_sig_len)];
% end
% 
% % ta
% if isfield(signalParams,'ta') && ~isempty(signalParams.ta)
%     filetagstr = [filetagstr, 'ta', num2str(signalParams.ta)];
% end

% % snr
% if isfield(signalParams,'snr') && ~isempty(signalParams.snr)
%     filetagstr = [filetagstr, 'snr', num2str(signalParams.snr)];
% end

% Drop fields that were missing or empty, then add one trailing underscore.
tags = tags(tags ~= "");

if isempty(tags)
    filetagstr = '';
else
    filetagstr = char(strjoin(tags, '_'));
end

end

function S = ensureStruct(S)
    if ~isstruct(S)
        S = loadjson(S);
    end
end

function tag = getTag(S, fieldName, prefix, formatSpec)
    tag = "";

    if ~isfield(S, fieldName) || isempty(S.(fieldName))
        return
    end

    value = S.(fieldName);

    % String-like fields, such as "name"
    if ischar(value) || isstring(value)
        tag = string(prefix) + string(value);
        return
    end

    % Numeric fields
    tag = string(prefix) + string(sprintf(formatSpec, value));
end

