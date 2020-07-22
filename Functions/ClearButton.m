% Author: Corentin Wicht, LCNS, 2018
% corentin.wicht@unifr.ch

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

function GUI = ClearButton (GUI)
    uicontrol('Style', 'pushbutton', 'String', 'Clear',...
    'Position', [330 20 50 20],'Callback',{@ClearCallback,GUI});

function ClearCallback(~,~,GUI)
    assignin('base','GUI',GUI)
    GUI.Data(:,2) = repmat({false},[length(GUI.Data(:,2)),1]);
end

end
