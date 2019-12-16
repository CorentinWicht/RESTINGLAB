function EEG=RemovEvents(EEG,varargin)

% This function computes :
% Usage:
%    >> EEG=RemovEvents(EEG,varargin);
%
% Inputs (mandatory):
%   EEG  = EEGLab dataset with .event and .urevent fields

% Inputs (Optional):
%   EventTypes  = cell of string(s) containing the event types to be
%                 removed

% Outputs:
%   EEG = Modified EEGLab dataset 

% Author: Corentin Wicht, LCNS, 2019
% - corentin.wicht@unifr.ch
% - https://github.com/CorentinWicht

%% SETTINGS
EventTypes = {'boundary'};

% Process Secondary Arguments
if nargin > 1
    for i = 1:2:length(varargin)
        Param = varargin{i};
        Value = varargin{i+1};
        if ~ischar(Param)
          error('Flag arguments must be strings')
        end
        Param = lower(Param);
    
        switch Param
            case 'eventtypes'
                EventTypes  = Value;
            otherwise
                display (['Unknown parameter setting: ' Param])
        end
    end
end

%% EVENTS REMOVAL

% Settings
Eventfields = {'event','urevent'};
TEMPEEG = EEG; % Using a temporary structure not to change to original 
%n data type of the event structure

% Changing all event types into string
for p=1:length(Eventfields)
    if ~isempty(TEMPEEG.(Eventfields{p}))
        for k=1:length(TEMPEEG.(Eventfields{p})) 
            if isnumeric(TEMPEEG.(Eventfields{p})(k).type)
                TEMPEEG.(Eventfields{p})(k).type = num2str(TEMPEEG.(Eventfields{p})(k).type);
           end
        end
    end
end

% Removing events
for p=1:length(Eventfields)
    if ~isempty(TEMPEEG.(Eventfields{p}))
        for k = 1:length(EventTypes)
            IdxBound = contains({TEMPEEG.(Eventfields{p}).type},EventTypes{k});
            EEG.(Eventfields{p})(IdxBound)=[];
        end
    end
end
