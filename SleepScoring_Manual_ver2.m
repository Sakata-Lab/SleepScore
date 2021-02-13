function varargout = SleepScoring_Manual_ver2(varargin)
%SLEEPSCORING_MANUAL_VER2 M-file for SleepScoring_Manual_ver2.fig
%      SLEEPSCORING_MANUAL_VER2, by itself, creates a new SLEEPSCORING_MANUAL_VER2 or raises the existing
%      singleton*.
%
%      H = SLEEPSCORING_MANUAL_VER2 returns the handle to a new SLEEPSCORING_MANUAL_VER2 or the handle to
%      the existing singleton*.
%
%      SLEEPSCORING_MANUAL_VER2('Property','Value',...) creates a new SLEEPSCORING_MANUAL_VER2 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to SleepScoring_Manual_ver2_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SLEEPSCORING_MANUAL_VER2('CALLBACK') and SLEEPSCORING_MANUAL_VER2('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SLEEPSCORING_MANUAL_VER2.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SleepScoring_Manual_ver2

% Last Modified by GUIDE v2.5 30-Sep-2014 16:48:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SleepScoring_Manual_ver2_OpeningFcn, ...
                   'gui_OutputFcn',  @SleepScoring_Manual_ver2_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%
% --- Executes just before SleepScoring_Manual_ver2 is made visible.
function SleepScoring_Manual_ver2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

%%% default parameters
handles.FileName=[];
handles.FilePath=[];
handles.SR=1000;

handles.sampsPerVolt=double(intmax('int16'))/5;
handles.timeWin=4; % 4 sec window (can be updated later)
handles.nWinsDisp=8; % the number of time window displayed
handles.ScoreFileName='sleepscoring_results.mat'; % contains Score(1D vect) and Log(2D mat:LogID,SegID,Score) 


% Choose default command line output for SleepScoring_Manual_ver2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SleepScoring_Manual_ver2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line.
function varargout = SleepScoring_Manual_ver2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Next.
function Next_Callback(hObject, eventdata, handles)
% hObject    handle to Next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Score and Log update
PrevScore=handles.Score(handles.currentWin);
CurrScore=get(handles.State,'Value');
if CurrScore==4
    CurrScore=0;
end
handles.Score(handles.currentWin)=CurrScore;
MyLogID=size(handles.Log,1)+1;
handles.Log(end+1,:)=[MyLogID handles.currentWin PrevScore CurrScore];

Score=handles.Score;
Log=handles.Log;
save(handles.ScoreFileName,'Score','Log');

%% update handles.currentWin
if handles.currentWin+1<=handles.TotalNumWin
    handles.currentWin=handles.currentWin+1;
    handles.Score(handles.currentWin)=handles.Score(handles.currentWin-1);
    handles=MyUpdate(handles);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in Previous.
function Previous_Callback(hObject, eventdata, handles)
% hObject    handle to Previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Score and Log update
PrevScore=handles.Score(handles.currentWin);
CurrScore=get(handles.State,'Value');
if CurrScore==4
    CurrScore=0;
end
handles.Score(handles.currentWin)=CurrScore;
MyLogID=size(handles.Log,1)+1;
handles.Log(end+1,:)=[MyLogID handles.currentWin PrevScore CurrScore];

Score=handles.Score;
Log=handles.Log;
save(handles.ScoreFileName,'Score','Log');

%% update handles.currentWin
if handles.currentWin-1>0
    handles.currentWin=handles.currentWin-1;
    handles=MyUpdate(handles);
end

% Update handles structure
guidata(hObject, handles);



% --- Executes on selection change in State.
function State_Callback(hObject, eventdata, handles)
% hObject    handle to State (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns State contents as cell array
%        contents{get(hObject,'Value')} returns selected item from State


% --- Executes during object creation, after setting all properties.
function State_CreateFcn(hObject, eventdata, handles)
% hObject    handle to State (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function Where_Callback(hObject, eventdata, handles)
% hObject    handle to Where (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Where as text
%        str2double(get(hObject,'String')) returns contents of Where as a double


% --- Executes during object creation, after setting all properties.
function Where_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Where (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on State and none of its controls.
function State_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to State (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% if get(handles,'CurrentCharacter')=='1'
%     set(handles.State,'Value',1);
% elseif get(handles,'CurrentCharacter')=='2'
%     set(handles.State,'Value',2);
% elseif get(handles,'CurrentCharacter')=='3'
%     set(handles.State,'Value',3);
% end

% if eventdata.Key==1
%     fprintf('pressed\n');
%     set(handles.State,'Value',1);
% elseif eventdata.Key=='2'
%     set(handles.State,'Value',2);
% elseif eventdata.Key=='3'
%     set(handles.State,'Value',3);
% end


% --------------------------------------------------------------------
function MyFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to MyFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.FileName,handles.PathName] = uigetfile('.dat','Open a target dat file');

%% recording info
handles.nCh=str2double(get(handles.NumCh,'String'));
handles.EEch=str2double(get(handles.EEGch,'String'));
handles.EMch1=str2double(get(handles.EMGch1,'String'));
handles.EMch2=str2double(get(handles.EMGch2,'String'));
handles.Ref=get(handles.EMGref,'Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load file and display
cd(handles.PathName);
% get end point
fileInfo = dir(handles.FileName);
fileSize = fileInfo.bytes;
handles.EndTime=fileSize/2/handles.nCh/handles.SR; % rec duration (sec)
handles.TotalNumWin=floor(handles.EndTime/handles.timeWin);

%% load previous job
if fopen(handles.ScoreFileName)~=-1
    load(handles.ScoreFileName); % Score, Log
    handles.Score=Score;
    handles.Log=Log;
    if ~isempty(Log)
        handles.currentWin=Log(end,2);
    else
        handles.currentWin=1;
    end
else
    handles.Score=zeros(handles.TotalNumWin,1);
    handles.Log=[];
    handles.currentWin=1;
end


%% load data
handles=MyUpdate(handles);

Log=handles.Log;
Score=handles.Score;
save(handles.ScoreFileName,'Score','Log');

set(handles.Slide,'SliderStep',[length(handles.StartWin:handles.FinishWin)/handles.TotalNumWin 0.1]);

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FigName=uiputfile('.fig');

% Update handles structure
guidata(hObject, handles);


function NumCh_Callback(hObject, eventdata, handles)
% hObject    handle to NumCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumCh as text
%        str2double(get(hObject,'String')) returns contents of NumCh as a double


% --- Executes during object creation, after setting all properties.
function NumCh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function Slide_Callback(hObject, eventdata, handles)
% hObject    handle to Slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(hObject,'SliderStep',[length(handles.StartWin:handles.FinishWin)/handles.TotalNumWin 0.1]);

MySlidePosition=get(hObject,'Value');
if MySlidePosition==0
    handles.currentWin=1;
else
    handles.currentWin=ceil(handles.TotalNumWin*MySlidePosition);
end

%% update handles.currentWin
handles=MyUpdate(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function EEGx_Callback(hObject, eventdata, handles)
% hObject    handle to EEGx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EEGx as text
%        str2double(get(hObject,'String')) returns contents of EEGx as a double


% --- Executes during object creation, after setting all properties.
function EEGx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EEGx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EMGx_Callback(hObject, eventdata, handles)
% hObject    handle to EMGx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EMGx as text
%        str2double(get(hObject,'String')) returns contents of EMGx as a double


% --- Executes during object creation, after setting all properties.
function EMGx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EMGx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offset_Callback(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset as text
%        str2double(get(hObject,'String')) returns contents of offset as a double


% --- Executes during object creation, after setting all properties.
function offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function EEGch_Callback(hObject, eventdata, handles)
% hObject    handle to EEGch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EEGch as text
%        str2double(get(hObject,'String')) returns contents of EEGch as a double


% --- Executes during object creation, after setting all properties.
function EEGch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EEGch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EMGch1_Callback(hObject, eventdata, handles)
% hObject    handle to EMGch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EMGch1 as text
%        str2double(get(hObject,'String')) returns contents of EMGch1 as a double


% --- Executes during object creation, after setting all properties.
function EMGch1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EMGch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EMGch2_Callback(hObject, eventdata, handles)
% hObject    handle to EMGch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EMGch2 as text
%        str2double(get(hObject,'String')) returns contents of EMGch2 as a double


% --- Executes during object creation, after setting all properties.
function EMGch2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EMGch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EMGref.
function EMGref_Callback(hObject, eventdata, handles)
% hObject    handle to EMGref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EMGref



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% updating function
function handles=MyUpdate(handles)

%% update channels
handles.EEch=str2double(get(handles.EEGch,'String'));
handles.EMch1=str2double(get(handles.EMGch1,'String'));
handles.EMch2=str2double(get(handles.EMGch2,'String'));
handles.Ref=get(handles.EMGref,'Value');

%% load data
%% determine start and end points to display
if handles.currentWin-handles.nWinsDisp<=0
    handles.StartPoint=1;
    handles.StartWin=1;
    x0=0; % in sec
else
    handles.StartPoint=(handles.currentWin-handles.nWinsDisp-1)*handles.timeWin*handles.SR+1;
    handles.StartWin=handles.currentWin-handles.nWinsDisp;
    x0=(handles.currentWin-handles.nWinsDisp-1)*handles.timeWin;
end
if (handles.currentWin+handles.nWinsDisp)>handles.TotalNumWin
    handles.FinishPoint=handles.TotalNumWin*handles.timeWin*handles.SR;
    handles.FinishWin=handles.TotalNumWin;
    xEnd=handles.TotalNumWin*handles.timeWin;
else
    handles.FinishPoint=(handles.currentWin+handles.nWinsDisp)*handles.timeWin*handles.SR;
    handles.FinishWin=handles.currentWin+handles.nWinsDisp;
    xEnd=(handles.currentWin+handles.nWinsDisp)*handles.timeWin;
end

x=x0:1/handles.SR:xEnd;
x(end)=[];
cX=[(handles.currentWin-1)*handles.timeWin handles.currentWin*handles.timeWin];

%% load data
Data = MyLoadDatSeg(handles.FileName, handles.nCh, handles.StartPoint, handles.FinishPoint)/handles.sampsPerVolt;
handles.eeg=Data(handles.EEch,:); 
handles.emg1=Data(handles.EMch1,:);
handles.emg2=Data(handles.EMch2,:);

if handles.Ref==1
    handles.emg=handles.emg1-handles.emg2;
else
    handles.emg=handles.emg1;
end

%% display 
EEGx=str2double(get(handles.EEGx,'String'));
EMGx=str2double(get(handles.EMGx,'String'));
Offset=str2double(get(handles.offset,'String'));

%% compute range
handles.Yrange=[-1*EMGx-Offset EEGx*1];
Yr=handles.Yrange;
%fill([cX(1) cX(2) cX(2) cX(1)],[Yr(2) Yr(2) Yr(1) Yr(1)],'FaceColor',[255,153,153],'Parent',handles.EEMGplot);
plot(handles.EEMGplot,[cX(1) cX(1)],[Yr(2) Yr(1)],'r:','LineWidth',5);
hold(handles.EEMGplot,'on');
plot(handles.EEMGplot,[cX(2) cX(2)],[Yr(2) Yr(1)],'r:','LineWidth',5);
plot(handles.EEMGplot,x,EEGx*handles.eeg,'k');
plot(handles.EEMGplot,x,EMGx*handles.emg-1*Offset,'k');
hold(handles.EEMGplot,'off');
axis(handles.EEMGplot,[x(1) x(end) Yr]);
set(handles.EEMGplot,'XTick',x(1):handles.timeWin:x(end));

%% score (local)
MyScore=handles.Score(handles.StartWin:handles.FinishWin);
MyCurrentScore=handles.Score(handles.currentWin);

bar(handles.score,[handles.StartWin:handles.FinishWin],MyScore,'k');
hold(handles.score,'on');
plot(handles.score,[handles.currentWin handles.currentWin],[-0.1 3.1],'r:','LineWidth',5);
%bar(handles.score,handles.currentWin,MyCurrentScore,'r');
hold(handles.score,'off');
axis(handles.score,[handles.StartWin-0.5 handles.FinishWin+0.5 -0.1 3.1]);
set(handles.score,'XTick',handles.StartWin:handles.FinishWin);
    
%% score (overall)
bar(handles.TotalScore,handles.Score,'k');
hold(handles.TotalScore,'on');
plot(handles.TotalScore,[handles.currentWin handles.currentWin],[-0.1 3.1],'r:','LineWidth',5);
hold(handles.TotalScore,'off');
axis(handles.TotalScore,[0 handles.TotalNumWin+1 -0.1 3.1]);

%% FFT (overall in the window)
L=length(handles.eeg);
nFFT=2^nextpow2(L);
Y=fft(handles.eeg,nFFT)/L;
f=handles.SR/2*linspace(0,1,nFFT/2+1);
Power=2*abs(Y(1:nFFT/2+1));
plot(handles.FFT,f,Power/max(Power),'k');
axis(handles.FFT,[0 20 0 1]);
xlabel(handles.FFT,'Freq (Hz)');
ylabel(handles.FFT,'Normalized Power');

%% FFT (current segment)
Data = MyLoadDatSeg(handles.FileName, handles.nCh, handles.currentWin*handles.timeWin*handles.SR+1,(handles.currentWin+1)*handles.timeWin*handles.SR)/handles.sampsPerVolt;
CurrentEEG=Data(1,:); %% ASSUMPTION

L=length(CurrentEEG);
nFFT=2^nextpow2(L);
Y=fft(CurrentEEG,nFFT)/L;
f=handles.SR/2*linspace(0,1,nFFT/2+1);
Power=2*abs(Y(1:nFFT/2+1));
hold(handles.FFT,'on');
plot(handles.FFT,f,Power/max(Power),'r','LineWidth',2);
hold(handles.FFT,'off');
legend(handles.FFT,'entire signals','current target'); 

%% update values
MyState=handles.Score(handles.currentWin);
if MyState==0
    MyState=4;
end
set(handles.State,'Value',MyState);
set(handles.Slide,'Value',(handles.currentWin-1)/handles.TotalNumWin);

set(handles.Where,'String',[num2str(handles.currentWin),' / ',num2str(handles.TotalNumWin)]);

%%%%%%%%%%%%
function SegDat = MyLoadDatSeg(FileName, nCh, StartPoint, FinishPoint)

% SegDat (ch#, datapoint)

StartByte=2*nCh*StartPoint;	% 16bits --> 2bytes
SampleToLoad = FinishPoint-StartPoint + 1;	%caution: It does include finishPoint data!

if StartPoint == 1
	SegDat = bload(FileName, [nCh, SampleToLoad]);
else
	SegDat = bload(FileName, [nCh, SampleToLoad], StartByte);
end
