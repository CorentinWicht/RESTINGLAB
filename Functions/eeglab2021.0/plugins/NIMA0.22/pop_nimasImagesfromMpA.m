% pop_nimasImagesfromMpA() - a wrapper function for nimasImagesfromMpA().

% Copyright (C) 2018, Makoto Miyakoshi and Nima Bigdely-Shamlo, SCCN, INC, UCSD.
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

% History
% 12/05/2018 Makoto. Created for Christiane Glatz.

function varargout = pop_nimasImagesfromMpA(varargin)
% POP_NIMASIMAGESFROMMPA MATLAB code for pop_nimasImagesfromMpA.fig
%      POP_NIMASIMAGESFROMMPA, by itself, creates a new POP_NIMASIMAGESFROMMPA or raises the existing
%      singleton*.
%
%      H = POP_NIMASIMAGESFROMMPA returns the handle to a new POP_NIMASIMAGESFROMMPA or the handle to
%      the existing singleton*.
%
%      POP_NIMASIMAGESFROMMPA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_NIMASIMAGESFROMMPA.M with the given input arguments.
%
%      POP_NIMASIMAGESFROMMPA('Property','Value',...) creates a new POP_NIMASIMAGESFROMMPA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pop_nimasImagesfromMpA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pop_nimasImagesfromMpA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pop_nimasImagesfromMpA

% Last Modified by GUIDE v2.5 06-Dec-2018 16:26:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pop_nimasImagesfromMpA_OpeningFcn, ...
                   'gui_OutputFcn',  @pop_nimasImagesfromMpA_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pop_nimasImagesfromMpA is made visible.
function pop_nimasImagesfromMpA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pop_nimasImagesfromMpA (see VARARGIN)

% Choose default command line output for pop_nimasImagesfromMpA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pop_nimasImagesfromMpA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pop_nimasImagesfromMpA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fwhmEdit_Callback(hObject, eventdata, handles)
% hObject    handle to fwhmEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fwhmEdit as text
%        str2double(get(hObject,'String')) returns contents of fwhmEdit as a double


% --- Executes during object creation, after setting all properties.
function fwhmEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fwhmEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in voxelSizePopupmenu.
function voxelSizePopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to voxelSizePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns voxelSizePopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from voxelSizePopupmenu


% --- Executes during object creation, after setting all properties.
function voxelSizePopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voxelSizePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in blobOrVoxelPopupmenu.
function blobOrVoxelPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to blobOrVoxelPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns blobOrVoxelPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from blobOrVoxelPopupmenu


% --- Executes during object creation, after setting all properties.
function blobOrVoxelPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blobOrVoxelPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alphaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to alphaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphaEdit as text
%        str2double(get(hObject,'String')) returns contents of alphaEdit as a double


% --- Executes during object creation, after setting all properties.
function alphaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function clusterIndicesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to clusterIndicesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clusterIndicesEdit as text
%        str2double(get(hObject,'String')) returns contents of clusterIndicesEdit as a double


% --- Executes during object creation, after setting all properties.
function clusterIndicesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusterIndicesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function optionalInputEdit_Callback(hObject, eventdata, handles)
% hObject    handle to optionalInputEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of optionalInputEdit as text
%        str2double(get(hObject,'String')) returns contents of optionalInputEdit as a double


% --- Executes during object creation, after setting all properties.
function optionalInputEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to optionalInputEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in projectOnlyOnXyPlaneCheckbox.
function projectOnlyOnXyPlaneCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to projectOnlyOnXyPlaneCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of projectOnlyOnXyPlaneCheckbox


% --- Executes on button press in helpPushbutton.
function helpPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to helpPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pophelp('nimasImagesfromMpA.m');



% --- Executes on button press in plotPushbutton.
function plotPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plotPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain STUDY and ALLEEG
STUDY = evalin('base', 'STUDY');
ALLEEG = evalin('base', 'ALLEEG');

% Obtain cluster indices.
clusterIndicesToPlot = str2num(get(handles.clusterIndicesEdit, 'String'));

% Obtain cluster dipoles.
clusterDipoles = std_readdipoles(STUDY, ALLEEG, clusterIndicesToPlot);

% Obtain cluster dipole locations.
clusterDipoleLocations = cell(1,length(clusterIndicesToPlot));
for clsIdx = 1:length(clusterIndicesToPlot)
    currentClsDip = clusterDipoles{clsIdx};
    
    % Dual dipoles are counted as two.
    dipList = [];
    for dipIdx = 1:length(currentClsDip)
        currentDipList = round(currentClsDip(dipIdx).posxyz);
        dipList = cat(1, dipList, currentDipList);
    end
    
    % Store as cluster dipole locations.
    clusterDipoleLocations{1,clsIdx} = dipList;
end

% Pass the dipole location info to the visualization function.
voxelSize      = get(handles.voxelSizePopupmenu, 'value')+1;
FWHM           = str2num(get(handles.fwhmEdit, 'string'));
blobOrVoxelIdx = get(handles.blobOrVoxelPopupmenu, 'value');
uniformAlpha   = str2num(get(handles.alphaEdit, 'string'));

% Parse optional input to construct string-value pairs if it is not empty.
optionStr = get(handles.optionalInputEdit, 'string');
if isempty(optionStr)
    optionalInput = {};
else
    commaIdx = strfind(optionStr, ',');
    for optionalInputIdx = 1:length(commaIdx)+1
        
        if optionalInputIdx == 1;
            currentStrings = optionStr(1:commaIdx(1)-1);
        elseif optionalInputIdx == length(commaIdx)+1
            currentStrings = optionStr(commaIdx(end):length(optionStr));
        else
            currentStrings = optionStr(commaIdx(optionalInputIdx-1)+1:commaIdx(optionalInputIdx)-1);
        end
        currentStrings = strtrim(currentStrings);
        
        
        if mod(optionalInputIdx,2) == 0 % Odd numbers are strings, even numbers are values.
            optionalInput{optionalInputIdx} = str2num(currentStrings);
        else
            optionalInput{optionalInputIdx} = currentStrings(2:end-1);
        end
    end
end
    
nimasImagesfromMpA(clusterDipoleLocations, 'voxelSize', voxelSize,...
                                           'FWHM', FWHM,...
                                           'blobOrVoxel', blobOrVoxelIdx,...
                                           'uniformAlpha', uniformAlpha,...
                                           optionalInput);
disp('Done.')