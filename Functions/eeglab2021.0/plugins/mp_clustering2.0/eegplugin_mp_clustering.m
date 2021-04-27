% eegplugin_mp_clustering() - Add both Measure Product pre-clustering AND Clustering commands to
% STUDY menu.
function vers = eegplugin_mp_clustering( fig, try_strings, catch_strings)
 
% title shown when loading into eeglab at the beginning
vers = 'mp_clustering2.0';

% create menu
plotmenu = findobj(fig, 'tag', 'study');

% make a sub-menu under study
submenu = uimenu( plotmenu, 'Label', 'Measure Product clustering', 'position', 6, 'userdata', 'startup:off;study:on');

uimenu( submenu, 'label', 'Build preclustering matrices', 'callback', [ try_strings.no_check '[STUDYTMP ALLEEGTMP LASTCOM]= pop_mpreclust(STUDY, ALLEEG);' catch_strings.update_study], 'userdata', 'startup:off;study:on');
uimenu( submenu, 'label', 'Cluster components', 'callback', [ try_strings.no_check '[STUDYTMP ALLEEGTMP LASTCOM] = pop_mpcluster(STUDY, ALLEEG);' catch_strings.update_study], 'userdata', 'startup:off;study:on');
