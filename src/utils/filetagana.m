function filetagstr = filetagana(psoParams,signalParams)
% T = filetag(P,S)
%Filetag creator for ana-related functions using specific parameters from P
%and S. Order of functions should be filetagana->createPSD->rungwpso.
%Typical naming procedure will have outFilePrfx = [outFilePath,filetagstr].
% P is a JSON file with parameters relating to pso's use, S is a JSON file
% with parameters relating to the injected signal.
%
%Created May. 2023 by Thomas Cruz from DRASE/filetag.m
%Updated Sept. 2026

%% Structure Check
psoParams = ensureStruct(psoParams);
signalParams   = ensureStruct(signalParams);

%% Build filename tag
tags = strings(0, 1);

% From psoParams
tags(end+1) = getTag(psoParams, 'type',       '',    '%s');
tags(end+1) = getTag(psoParams, 'maxSteps',    'stp',   '%g');

% From signalParams
tags(end+1) = getTag(signalParams, 'sampling_freq', 'fs',  '%g');
tags(end+1) = getTag(signalParams.signal, 'T_sig_len',        'tsL',   '%g');
tags(end+1) = getTag(signalParams, 'ta',        'ta',   '%g');
tags(end+1) = getTag(signalParams, 'snr',        'snr',   '%g');

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

