% eegplugin_EEGBrowser()
% 
% Adds the EEGBrowser plugin to the path.
%
% Author: Alejandro Ojeda, Neural Engineering and Translation Lab, 
% Department of Psychiatry, UC San Diego, 2019

function v = eegplugin_EEGBrowser(fig,try_strings, catch_strings)

v = '1.1';
p = fileparts(which('eegplugin_EEGBrowser'));
addpath(p);
hPlot = findobj(fig,'Label','Plot');
uimenu( hPlot, 'label', 'Channel data (EEG Browser)', 'callback','pop_eegbrowser(EEG);', 'position', 3);
uimenu( hPlot, 'label', 'Component activations (EEG Browser)', 'callback','pop_eegbrowser(EEG,0);', 'position', 10);
