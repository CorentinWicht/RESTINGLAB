% Erase all data in all sheets in the workbook.
function EraseExcelSheets(fileName)
% By ImageAnalyst
% Check whether the file exists
if ~exist(fileName,'file')
error([fileName ' does not exist !']);
else
% Check whether it is an Excel file
typ = xlsfinfo(fileName);
if ~strcmp(typ,'Microsoft Excel Spreadsheet')
error([fileName ' not an Excel sheet !']);
end
end

% If fileName does not contain a "\" the name of the current path is
% added to fileName. The reason for this is that the full path is required
% for the command "excelObj.workbooks.Open(fileName)" to work properly.
if ~contains(fileName,'\')
fileName = [cd '\' fileName];
end

excelObj = actxserver('Excel.Application');
excelWorkbook = excelObj.workbooks.Open(fileName);
worksheets = excelObj.sheets;
numSheets = worksheets.Count;
% Prevent beeps from sounding if we try something not allowed.
excelObj.EnableSound = false;

% Loop over all sheets
    for s = 1 : numSheets
        % worksheets.Item(sheetIndex).UsedRange is the range of used cells
        worksheets.Item(s).UsedRange.Delete;
    end
excelObj.EnableSound = true;
excelWorkbook.Save;
excelWorkbook.Close(false);
excelObj.Quit;
delete(excelObj);
return;