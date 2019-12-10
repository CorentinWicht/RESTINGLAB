% eegplugin_fitTwoDipoles()--see help of IC_sym_class.m

% 05/03/2016 Makoto. Created.
%
% Copyright (C) 2016 Makoto Miyakoshi, SCCN,INC,UCSD
%                    Caterina Piazza, Politecnico di Milano/Scientific Institute IRCCS Eugenio Medea
%
% Reference: 
% Piazza, Miyakoshi, Akalin-Acar, Cantiani, Reni, Bianchi, Makeig. (2016). An automated
% function for identifying EEG independent components representing bilateral source
% activity. XIV Mediterranean Conference on Medical and Biological Engineering and
% Computing 2016 IFMBE Proceedings 57:105-109.
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function vers = eegplugin_fitTwoDipoles(fig, trystrs, catchstrs)
    
vers = 'beta';
if nargin < 3
    error('eegplugin_fitTwoDipoles requires 3 arguments');
end;
    
% create menu
toolsmenu = findobj(fig, 'tag', 'tools');
uimenu( toolsmenu, 'label', 'Fit bilateral dipoles', 'separator','on','callback', 'pop_fitTwoDipoles');