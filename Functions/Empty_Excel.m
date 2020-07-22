function Empty_Excel(Directory,Extension,ImportParameters)
%-------------------------------------------------------------------------%
%  EMPTY EXCEL FILES
%-------------------------------------------------------------------------%

% This function enables the user to erase the content of all or selected
% excel files. It also saves a matrix (.mat) including deletion
% parameters for further use. 

% Usage:
%    >> Empty_Excel(Directory,Extension);
%    >> Empty_Excel(Directory,Extension,ImportParameters);
%
% Inputs:
%   Directory        = Most upper folder containing all the excel files
%   Extension        = Extension of excel files (e.g. xlsx, csv)

% Optional Input:
%   ImportParameters = .Mat file saved from filling the uitable in this
%                       function. This will save you a lot of time, the
%                       next time you run the function.

% Author: Corentin Wicht, LCNS, 2018
% corentin.wicht@unifr.ch

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)
    

% Include possibility to have more than 1 extensions!!!

% Avoid crashing if no parameters specified
if nargin < 3
    ImportParameters='';
end

% Temporarily stores the current path to reset it at the end
TempDirectory=pwd;

% Current directory = where the excel files are stored
cd(Directory)

% Defines the excel files extension to look for
if ~ismember('.', Extension)
    Extension=['.' lower(Extension)];
end

% Retrieve all files with corresponding extension in the current folder
% (and substructures)
FileList = dir(['**/*' Extension]);

% Uitable to select which excel file to empty
[AllNames, NamesIndex]=unique({FileList.name});

% Prompt for participants selection
PromptPartSelect = questdlg('Would you like to erase the content of i) all excel files, ii) selected one or iii) using loaded parameters ?', ...
    'Files selection','All','Selected','Loaded Parameters','Selected');
if strcmp(PromptPartSelect,'All')
    LogicalSelect={true};
elseif strcmp(PromptPartSelect,'Selected')
    LogicalSelect={false};
end
if ~strcmp(PromptPartSelect,'Loaded Parameters')
    % Matrix to integrate in the following uitable
    to_display=[AllNames', repmat(LogicalSelect,[size(AllNames,2) 1]), repmat({' '},[size(AllNames,2) 1])];
else
    % Importing saved .mat file
    if ~ismember('.',ImportParameters)
        ImportFile=[ImportParameters '.mat'];
    end
    [file,path] = uigetfile(ImportFile);
    File = fullfile(path, file);
    to_display_temp = struct2cell(load(File));
    to_display=to_display_temp{:};
end
% In case nothing changes in the uitable
MarkerList=to_display;

% Determines width of columns
Column1=Adapt_Width(AllNames);
Column2=Adapt_Width({'Empty file ?'});
Column3=Adapt_Width({'Range (leave empty if all content need to be emptied)'}) ;

% Select folders on which to apply analyses 
f=figure('Position', [150 150 800 400]);
Table1=uitable('Parent', f,'Data',to_display,'ColumnEdit',[false true true],'ColumnName',...
    {'Files', 'Empty file ?','Range (leave empty if all content need to be emptied)'},'CellEditCallBack','MarkerList = get(gco,''Data'');');
uicontrol('Style', 'text', 'Position', [0 325 400 50], 'String',...
        {'Files to empty','Click on the box of the files you want to empty and provide the excel range (e.g. B2:C4)'});
Table1.Position = [0 0 800 320];    
set(Table1,'ColumnWidth', {Column1,Column2,Column3});
   
% Wait for t to close until running the rest of the script
waitfor(Table1)

% Build shit
TempFolder={FileList.folder}';
FolderCells=TempFolder(NamesIndex,:);
ToErase=[FolderCells MarkerList];

% Saves the design in a .mat file
uisave('MarkerList','ExcelErasing_Parameters.mat')

% Processes the selected files
for i=1:size(ToErase,1)
    if ToErase{i,3}==true
        if ~strcmp(ToErase{i,4},' ')
            EraseExcelRange(backslash([ToErase{i,1} '\' ToErase{i,2}]),ToErase{i,4});
        else
            EraseExcelSheets(backslash([ToErase{i,1} '\' ToErase{i,2}]));
        end
        disp([ToErase{i,2} '...content has been emptied'])
    end
end

% Sets directory to beginning value
cd(TempDirectory)
end