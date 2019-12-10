function varargout = pop_fitTwoDipoles(varargin)
% POP_FITTWODIPOLES MATLAB code for pop_fitTwoDipoles.fig
%      POP_FITTWODIPOLES, by itself, creates a new POP_FITTWODIPOLES or raises the existing
%      singleton*.
%
%      H = POP_FITTWODIPOLES returns the handle to a new POP_FITTWODIPOLES or the handle to
%      the existing singleton*.
%
%      POP_FITTWODIPOLES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_FITTWODIPOLES.M with the given input arguments.
%
%      POP_FITTWODIPOLES('Property','Value',...) creates a new POP_FITTWODIPOLES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pop_fitTwoDipoles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pop_fitTwoDipoles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pop_fitTwoDipoles

% Last Modified by GUIDE v2.5 03-May-2016 16:45:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pop_fitTwoDipoles_OpeningFcn, ...
                   'gui_OutputFcn',  @pop_fitTwoDipoles_OutputFcn, ...
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


% --- Executes just before pop_fitTwoDipoles is made visible.
function pop_fitTwoDipoles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pop_fitTwoDipoles (see VARARGIN)

% Choose default command line output for pop_fitTwoDipoles
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pop_fitTwoDipoles wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pop_fitTwoDipoles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in SymmetryRegionPopupmenu.
function SymmetryRegionPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to SymmetryRegionPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SymmetryRegionPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SymmetryRegionPopupmenu


% --- Executes during object creation, after setting all properties.
function SymmetryRegionPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SymmetryRegionPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresholdEdit_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdEdit as a double


% --- Executes during object creation, after setting all properties.
function thresholdEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in startPushbutton.
function startPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to startPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load EEG from the base workspace.
EEG = evalin('base', 'EEG');

% Fetch User Inputs
symmetryRegionIdx = get(handles.SymmetryRegionPopupmenu, 'value');
switch symmetryRegionIdx
    case 1
        symmetryRegion = 'LRR';
    case 2
        symmetryRegion = 'LSR';
    case 3
        symmetryRegion = 'SRR';
    case 4
        symmetryRegion = 'SSR';
end
threshold = str2num(get(handles.thresholdEdit, 'String'));

% Run Caterina's solution.
EEG = fitTwoDipoles(EEG, symmetryRegion, threshold);

% Generate eegh log.
com = ['EEG = fitTwoDipoles(EEG, ''' symmetryRegion ''', ' num2str(threshold) ');'];

% Store the history to EEG
EEG = eegh(com, EEG);

% Assign-in EEG to the base workspace.
assignin('base', 'EEG', EEG);

% Display end message.
disp('Done.')
