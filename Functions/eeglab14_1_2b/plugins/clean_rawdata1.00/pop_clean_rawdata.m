% History:
% 03/26/2019 Makoto and Chiyuan. Major update. Supported 'availableRAM_GB'. GUI switched to GUIDE-made.

function varargout = pop_clean_rawdata(varargin)
% POP_CLEAN_RAWDATA MATLAB code for pop_clean_rawdata.fig
%      POP_CLEAN_RAWDATA, by itself, creates a new POP_CLEAN_RAWDATA or raises the existing
%      singleton*.
%
%      H = POP_CLEAN_RAWDATA returns the handle to a new POP_CLEAN_RAWDATA or the handle to
%      the existing singleton*.
%
%      POP_CLEAN_RAWDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_CLEAN_RAWDATA.M with the given input arguments.
%
%      POP_CLEAN_RAWDATA('Property','Value',...) creates a new POP_CLEAN_RAWDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pop_clean_rawdata_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pop_clean_rawdata_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pop_clean_rawdata

% Last Modified by GUIDE v2.5 26-Mar-2019 18:54:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pop_clean_rawdata_OpeningFcn, ...
                   'gui_OutputFcn',  @pop_clean_rawdata_OutputFcn, ...
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


% --- Executes just before pop_clean_rawdata is made visible.
function pop_clean_rawdata_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pop_clean_rawdata (see VARARGIN)

% Choose default command line output for pop_clean_rawdata
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pop_clean_rawdata wait for user response (see UIRESUME)
% uiwait(handles.pop_clean_rawdataGUI);


% --- Outputs from this function are returned to the command line.
function varargout = pop_clean_rawdata_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function flatChannelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to flatChannelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flatChannelEdit as text
%        str2double(get(hObject,'String')) returns contents of flatChannelEdit as a double


% --- Executes during object creation, after setting all properties.
function flatChannelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flatChannelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poorCorrChanEdit_Callback(hObject, eventdata, handles)
% hObject    handle to poorCorrChanEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of poorCorrChanEdit as text
%        str2double(get(hObject,'String')) returns contents of poorCorrChanEdit as a double


% --- Executes during object creation, after setting all properties.
function poorCorrChanEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poorCorrChanEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lineNoiseChanEdit_Callback(hObject, eventdata, handles)
% hObject    handle to lineNoiseChanEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lineNoiseChanEdit as text
%        str2double(get(hObject,'String')) returns contents of lineNoiseChanEdit as a double


% --- Executes during object creation, after setting all properties.
function lineNoiseChanEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lineNoiseChanEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function asrEdit_Callback(hObject, eventdata, handles)
% hObject    handle to asrEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of asrEdit as text
%        str2double(get(hObject,'String')) returns contents of asrEdit as a double


% --- Executes during object creation, after setting all properties.
function asrEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to asrEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function removeWindowEdit_Callback(hObject, eventdata, handles)
% hObject    handle to removeWindowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of removeWindowEdit as text
%        str2double(get(hObject,'String')) returns contents of removeWindowEdit as a double


% --- Executes during object creation, after setting all properties.
function removeWindowEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to removeWindowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in showResultsPopupmenu.
function showResultsPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to showResultsPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns showResultsPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from showResultsPopupmenu


% --- Executes during object creation, after setting all properties.
function showResultsPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to showResultsPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in helpPushbutton.
function helpPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to helpPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

doc('clean_rawdata')



% --- Executes on button press in cancelPushbutton.
function cancelPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Close the GUI.
close(handles.pop_clean_rawdataGUI)
disp('pop_clean_rawdata() cancelled by the user.')


function highpassEdit_Callback(hObject, eventdata, handles)
% hObject    handle to highpassEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of highpassEdit as text
%        str2double(get(hObject,'String')) returns contents of highpassEdit as a double


% --- Executes during object creation, after setting all properties.
function highpassEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highpassEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function optionEdit_Callback(hObject, eventdata, handles)
% hObject    handle to optionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of optionEdit as text
%        str2double(get(hObject,'String')) returns contents of optionEdit as a double


% --- Executes during object creation, after setting all properties.
function optionEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to optionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in obtainRamPushbutton.
function obtainRamPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to obtainRamPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

availableRamInGB = hlp_memfree/(2^30);
set(handles.optionEdit, 'String', sprintf('''availableRAM_GB'', %.2f', availableRamInGB))



% --- Executes on button press in okPushbutton.
function okPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to okPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain EEG from the base workspace. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = evalin('base', 'EEG');

% Delete EEG.etc.clean_channel_mask and EEG.etc.clean_sample_mask if present.
if isfield(EEG.etc, 'clean_channel_mask')
    EEG.etc = rmfield(EEG.etc, 'clean_channel_mask');
    disp('EEG.etc.clean_channel_mask present: Deleted.')
end
if isfield(EEG.etc, 'clean_sample_mask')
    EEG.etc = rmfield(EEG.etc, 'clean_sample_mask');
    disp('EEG.etc.clean_sample_mask present: Deleted.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert user inputs into numerical variables. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arg_flatline = get(handles.flatChannelEdit,'String');
if ~strcmp(arg_flatline, 'off')
    arg_flatline = str2num(arg_flatline);
end

arg_highpass = get(handles.highpassEdit,'String');
if ~strcmp(arg_highpass, 'off')
    arg_highpass = str2num(arg_highpass);
end

arg_channel = get(handles.poorCorrChanEdit,'String');
if ~strcmp(arg_channel, 'off')
    arg_channel = str2num(arg_channel);
end

arg_noisy = get(handles.lineNoiseChanEdit,'String');
if ~strcmp(arg_noisy , 'off')
    arg_noisy  = str2num(arg_noisy );
end

arg_burst = get(handles.asrEdit,'String');
if ~strcmp(arg_burst, 'off')
    arg_burst = str2num(arg_burst);
end

arg_window = get(handles.removeWindowEdit,'String');
if ~strcmp(arg_window, 'off')
    arg_window = str2num(arg_window);
end

arg_visartfc = get(handles.showResultsPopupmenu, 'Value');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build a package of optional inputs to make varargin. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optionalInput = hlp_split(get(handles.optionEdit, 'String'), ',');
optionalInputCells = {};
for itemIdx = 1:length(optionalInput)
    currentItem = strtrim(optionalInput{itemIdx}); % Trim the first and last space.
    if mod(itemIdx, 2) == 1
        optionalInputCells{itemIdx} = currentItem(2:end-1); % Trim one too many single quote from the first and the last.
    else
        optionalInputCells{itemIdx} = str2num(currentItem);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call the wrapper for Christian's main function. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cleanEEG = clean_rawdata(EEG, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst, arg_window, optionalInputCells);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use Christian's before and after comparison visualization. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if arg_visartfc == 1;
    try
        vis_artifacts(cleanEEG,EEG);
    catch
        warning('vis_artifacts failed. Skipping visualization.')
    end
end



% Update EEG.
EEG = cleanEEG;

% Output eegh.
com = EEG.etc.clean_rawdata_log;
EEG = eegh(com, EEG);

% Update EEG.
EEG.etc.clean_rawdata_log = com;
assignin('base', 'EEG', EEG);

% Close the GUI.
close(handles.pop_clean_rawdataGUI)
disp('Done.')
