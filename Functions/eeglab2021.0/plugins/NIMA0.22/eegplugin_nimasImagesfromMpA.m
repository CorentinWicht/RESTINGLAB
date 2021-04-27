% eegplugin_nimasImagesfromMpA() - This is a plugin for using Nima's Images 
%                                  from Measure-projection Analysis (NIMA).

% Copyright (C) 2018, Makoto Miyakoshi, Nima Bigdely-Shamlo. SCCN,INC,UCSD.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% History
% 12/04/2018 Makoto. Created for Christiane Glatz.


function eegplugin_nimasImagesfromMpA( fig, try_strings, catch_strings);

vers = '0.10';
    if nargin < 3
        error('eegplugin_nimasImagesfromMpA() requires 3 arguments');
    end
    
% create menu
std = findobj(fig, 'tag', 'study');
uimenu( std, 'label', 'NIMA', 'callback', 'pop_nimasImagesfromMpA', 'userdata', 'startup:off;study:on');