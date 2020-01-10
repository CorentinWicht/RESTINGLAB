function EEG=average_ref(EEG,varargin)

% This function corrects the original pop_reref function flaws. Namely,
% pop_reref, replaces the old reference electrode at the end of the
% EEG.chanlocs structure. The current code replaces the old reference
% electrode at its right place (according to EEG.chanlocs.urchan)

% Usage:
%    >> [EEG] = average_ref(EEG,varargin);
%
% Inputs (mandatory):
%   EEG       = EEG dataset structure

% Inputs (optional):
%   CzStruct  = EEG.chanlocs field of the old reference electrode

% Outputs:
%   EEG   - Modified EEG dataset structure

% Author: Corentin Wicht, LCNS, 2018

%% Process Secondary Arguments
CzStruct = [];

if nargin > 1
    for i = 1:2:length(varargin)
        Param = varargin{i};
        Value = varargin{i+1};
        if ~ischar(Param)
          error('Flag arguments must be strings')
        end
        Param = lower(Param);
    
        switch Param
            case 'n_permutes'
                N_Permutes  = Value;
            otherwise
                display (['Unknown parameter setting: ' Param])
        end
    end
end

%% Run
if ~isempty(CzStruct)
    CzStruct.ref=CzStruct.labels;
    CzStruct.datachan=0;
    EEG = pop_reref( EEG, [],'refloc',CzStruct);
else
    EEG = pop_reref( EEG, []);
end
[~,sortIndexes] = sort([EEG.chanlocs.urchan],'ascend');
EEG.chanlocs = EEG.chanlocs(sortIndexes);
EEG.data = EEG.data(sortIndexes,:,:);

end

