% Author: Corentin Wicht, LCNS, 2018
% corentin.wicht@unifr.ch

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

function Fig=MessageBox(Text,Title,FontSize,Width,Height)
% This function is an updated version of msgbox () enabling to resize the
% message box and the text inside

if nargin<3
    FontSize=16;
    Width=800;
    Height=300; 
elseif nargin<4
    Width=800;
    Height=300; 
elseif nargin<5
    Height=300; 
end

Fig=msgbox(Text,Title); 

ScreenDimensions=get(0,'ScreenSize');
% Change the size of the figure
object_handles = findall(Fig);
set( object_handles(4), 'FontSize', FontSize)

set(Fig, 'Units', 'Pixels', 'Position',...
        [ScreenDimensions(3)/2-Width/2 ScreenDimensions(4)/2-Height/2 Width Height]);  

waitfor(Fig)
%uiwait(Fig) 
end