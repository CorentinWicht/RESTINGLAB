function SaveFigures(fig,Directory,Col,Formats,Close)
%-------------------------------------------------------------------------%
% SAVING FULL SCREEN FIGURES
%-------------------------------------------------------------------------%
% This function enables to export full-screen pictures automically adapted
% to screen size. They will be exported in the directory specified by the
% argument Directory. Moreover, the functions can automatically save the
% figure in 3 different formats: pdf, bmp and eps. This function relies on
% another function called "backslash". 

% Usage:
%    >> [EEG] = SaveFigures(fig,Directory,Col, Formats);
%
% Inputs:
%   fig       = Figure to rescale and to save
%   Directory = Directory where to store the resulting figures 

% Optional input
%   Col       = Defines background color (e.g. 'b' for blue) (default is
%               white)
%   Formats   = Defines the format of export ('pdf', 'bmp', 'eps', 'fig', 'all')
%               (default is all three)
%   Close     = Force to close the current figure (def = 1, i.e. Yes).


% Author: Corentin Wicht, LCNS, 2018
% corentin.wicht@unifr.ch

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

if nargin<3
    Col = 'w';
    Formats='All'; 
elseif nargin<4
    Formats='All'; 
elseif nargin <5
    Close = 1;
end

set(fig,'Units','Inches');
pos = get(fig,'Position');  
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(fig,'color',Col);
fig.InvertHardcopy = 'off';

    if strcmpi(Formats,'All') 
        print(fig,backslash(Directory),'-dpdf','-r0'); 
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(fig,backslash(Directory),'bmp');
        saveas(fig,backslash(Directory),'epsc');
    elseif strcmpi(Formats,'pdf') 
        % Saves the figure as .pdf format 
        print(fig,backslash(Directory),'-dpdf','-r0'); 
    elseif strcmpi(Formats,'bmp')
        % Saves the figure as .bmp format 
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(fig,backslash(Directory),'bmp');
    elseif strcmpi(Formats,'eps') 
        % Saves the figure as .eps format 
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(fig,backslash(Directory),'epsc');
    elseif strcmpi(Formats,'fig') 
        % Saves the figure as MATLAB .fig format 
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        savefig(fig,backslash(Directory));
    end

    % Closes the current figure
    if Close
        close gcf
    end
end
