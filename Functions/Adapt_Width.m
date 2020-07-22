function MaxString=Adapt_Width(StringList)
%-------------------------------------------------------------------------%
%  EMPTY EXCEL FILES
%-------------------------------------------------------------------------%

% This function adapts the width of columns in uitable based on the length
% (in pixels) of the longest line.

% Usage:
%    >> MaxString=Adapt_Width(StringList);
%
% Inputs:
%   StringList       = Cell array of strings which will be measured in
%                      pixels and only the maximum value will be saved as 
%                      output. 

% Outputs:
%   MaxString = Variable containing the length of the longest line in
%               pixels (double)

% Author: Corentin Wicht, LCNS, 2018
% corentin.wicht@unifr.ch

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

% Determines the size of the longest string to display in uitable
PixelString=zeros(1,length(StringList));
for i=1:length(StringList)
    TempSize=uicontrol('Style', 'text','String',StringList{i},'Visible','off','Units','pixels');
    PixelString(i)=TempSize.Extent(3);
    close gcf
end
% Selects the maximum width
MaxString=max(PixelString);
end