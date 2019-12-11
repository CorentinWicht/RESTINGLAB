function EraseExcelRange(XlsFile, XlsRange)
%-------------------------------------------------------------------------%
% ERASE EXCEL RANGE
%-------------------------------------------------------------------------%
% This function erases excel files only in specified ranges.

% Usage:
%    >> EraseExcelRange(XlsFile, XlsRange)
%
% Inputs:
%   XlsFile   = Name of the excel file including the path and
%               extension.                    
%   XlsRange  = Range with format similar to Excel structure (e.g. B2:D49)

% Author: Corentin Wicht, LCNS, 2018
% corentin.wicht@unifr.ch

% Name of the excel file
filename = XlsFile;
% Open Excel as a COM Automation server
Excel = actxserver('Excel.Application');
% Open Excel workbook
Workbook = Excel.Workbooks.Open(filename);
% Load Excel sheets name
[~,Sheets] = xlsfinfo(filename); 
for i=1:length(Sheets)
    % Clear the content of the sheet
    Workbook.Worksheets.Item(Sheets{i}).Range(XlsRange).ClearContents
end
% Now save/close/quit/delete
Workbook.Save;
Excel.Workbook.Close;
invoke(Excel, 'Quit');
delete(Excel)
end
