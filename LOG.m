function LOG(SavePath,Date_Start,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%                             LOG SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% WRITE INSTRUCTIONS HERE !!

% Process Secondary Arguments
if nargin > 2
    for i = 1:2:length(varargin)
        Param = varargin{i};
        Value = varargin{i+1};
        if ~ischar(Param)
          error('Flag arguments must be strings')
        end
        Param = lower(Param);
    
        switch Param
            case 'logmst'
                LogMST  = Value;
            case 'logic'
                LogIC  = Value;
            case 'errorstudy'
                ErrorSTUDY = Value;
            case 'errorspreproc'
                ErrorArtifacts = Value;
            otherwise
                display (['Unknown parameter setting: ' Param])
        end
    end
end

%-------------------------------------------------------------------------%
% LOADING GUI MATRICES
%-------------------------------------------------------------------------%
% Retrieving the list of all matrices
ParametersPath=dir([SavePath '\Parameters\' Date_Start '\']);
ParametersPath=ParametersPath(~cell2mat({ParametersPath.isdir}));

% Loading the matrices stored with the GUI
for k=1:length(ParametersPath)
    load([ParametersPath(k).folder '\' ParametersPath(k).name]);
end

% If creating new parameters and only running GroupStudy, the
% MainWorkspace.mat might be missing. 
if nnz(contains({ParametersPath.name},'MainWorkspace.mat')) < 1
    % Get time and date
    CurrentDate = datenum(clock);
    
    % List of possible matrices
    TempList = dir([SavePath '\Parameters\' '**/*' 'MainWorkspace.mat']);
    
    if length(TempList) > 1
        for m=1:length(TempList)
            % Find the file whose date is closest to today
            SplitDateFolder = strsplit(TempList(m).folder,'\');
            SplitDateFolder = strrep(SplitDateFolder{end},'_','-');
            SplitDateFolder = strsplit(SplitDateFolder,'-');
            SplitDateFolder = [str2double(SplitDateFolder{3}) str2double(SplitDateFolder{2}) ...
                str2double(SplitDateFolder{1}) str2double(SplitDateFolder{4 }) ...
                str2double(SplitDateFolder{5}) 00];
            Datefolders(m) = datenum(SplitDateFolder);
        end
        
        % Compute distance from dates to now and take index of smallest
        [~, SmallDistDateIdx] = min(CurrentDate - Datefolders);
        load([TempList(SmallDistDateIdx).folder '\' TempList(SmallDistDateIdx).name]);
    else
        load([TempList.folder '\' TempList.name]);
    end
end

% Initialize variables
time_end = datestr(now);
username=getenv('USERNAME');
% Decision if doing basic analyses
if exist('AnalysesSwitch','var')
    Analysis = 1;
end

% Group analyses log data
if exist('LogMST','var')
    MSTVersions = fieldnames(LogMST);
end

% Converting Date_Start to comparable format
SplitDateFolder = strrep(Date_Start,'_','-');
SplitDateFolder = strsplit(SplitDateFolder,'-');
SplitDateFolder = [str2double(SplitDateFolder{3}) str2double(SplitDateFolder{2}) ...
    str2double(SplitDateFolder{1}) str2double(SplitDateFolder{4 }) ...
    str2double(SplitDateFolder{5}) 00];
time_start = datestr(SplitDateFolder);

% Creating the log file
% date_name = datestr(now,'dd-mm-yy_HHMM');
if sum(~SavePath)<1
    fid = fopen([SavePath '\RESTINGLablog_' Date_Start '.txt'],'w');    
else
    fid = fopen([CurrentPWD '\RESTINGLablog_' Date_Start '.txt'],'w');
end
% date, starting time, finished time, number of analyzed files
fprintf(fid,'%s\t%s\r\n',['Start : ',time_start],['End: ',time_end]);
fprintf(fid,'\r\n%s\r\n',['Windows username : ' username]);

% Summary of design ----------
BetweenFactors = BetweenFactors(~cellfun('isempty',BetweenFactors));
WithinFactors = WithinFactors(~cellfun('isempty',WithinFactors));

% Between subject
fprintf(fid,'\r\n\r\n%s\r\n','------ DESIGN SUMMARY ------');
if ~isempty(BetweenFactors)
    BetweenLevels = BetweenLevels(~cellfun('isempty',BetweenLevels));
    for k=1:length(BetweenFactors)
        LevelsTemp = '';
        for t=1:length(BetweenLevels(:,k))
            if t~=length(BetweenLevels(:,k))
                LevelsTemp = strcat(LevelsTemp, BetweenLevels{t,k}, ' Vs. ');
            else
                LevelsTemp = strcat(LevelsTemp, BetweenLevels{t,k});
            end
        end
        fprintf(fid,'\r\n%s',sprintf('Between-subject factor %d : %s, with levels : %s',...
            k,upper(BetweenFactors{k}),LevelsTemp)); 
    end
end

% Within subject
if ~isempty(WithinFactors)
    WithinLevels = WithinLevels(~cellfun('isempty',WithinLevels));
    for k=1:length(WithinFactors)
        fprintf(fid,'\r\n\r\n%s',sprintf('Within-subject factor %d : %s, with levels : %s',...
            k,upper(WithinFactors{k}),strjoin(WithinLevels(:,k))));
    end
end

% Summary of parameters ---------
fprintf(fid,'\r\n\r\n%s\r\n','------ PARAMETERS SUMMARY ------');
fprintf(fid,'\r\n%s',sprintf('Common files name: %s%s',FileNames, Extension)); 
fprintf(fid,'\r\n%s',sprintf('Sampling rate : %d Hz',SamplingRate));
fprintf(fid,'\r\n%s',sprintf('Channels list : %s',num2str(Channels))); 
if ~isempty(Beginning) || ~isempty(End)
   fprintf(fid,'\r\n%s',sprintf('Data length restriction : %d - %d (seconds)',Beginning, End));
else
    fprintf(fid,'\r\n%s',sprintf('Data length restriction : NO'));
end

% Summary of preprocessing -------
fprintf(fid,'\r\n\r\n%s\r\n','------ PREPROCESSING SUMMARY ------');
fprintf(fid,'\r\n%s',sprintf('Filtering frequencies : %d-%d Hz',HighPass,LowPass));
fprintf(fid,'\r\n%s',sprintf('Subjects included : %s',num2str(Subjectslist)));
fprintf(fid,'\r\n%s',sprintf('1) Reject asleep epochs : %s',upper(Sleep)));
if strcmpi(Sleep,'Yes')
    fprintf(fid,'\r\n%s',sprintf('Seconds/epoch : %d',SecondsEpoch));
end
fprintf(fid,'\r\n%s',sprintf('2) Make processing semi-automatic : %s',upper(Automaticity)));
fprintf(fid,'\r\n%s',sprintf('3) Compute ICA : %s',upper(ICA)));
fprintf(fid,'\r\n%s',sprintf('4) Use ASR : %s',upper(ASR)));
fprintf(fid,'\r\n%s',sprintf('4) Use BLINKER : %s',upper(BLINKER)));

if strcmpi(ICA,'Yes')
   % ICA PARAMETERS 
   fprintf(fid,'\r\n\r\n%s','ICA parameters ----');
   fprintf(fid,'\r\n%s',sprintf('- Number of iterations : %s',AnswerICA{:}));
end

if strcmpi(ASR,'Yes')
   % ASR PARAMETERS 
   fprintf(fid,'\r\n\r\n%s','ASR parameters ----');
   fprintf(fid,'\r\n%s',sprintf('- Repair bursts using ASR (std) : %s',AnswerASR{1}));
   fprintf(fid,'\r\n%s',sprintf('- Remove time windows (0-1) : %s',AnswerASR{2}));
end

if strcmpi(BLINKER,'Yes')
   % BLINKER PARAMETERS 
   BlinkerParam = {'Reject portions of data containing blinks','Interpolate portions containing blinks using ASR'};
   fprintf(fid,'\r\n\r\n%s','BLINKER parameters ----');
   fprintf(fid,'\r\n%s',sprintf('- BLINKER applied as following : %s',BlinkerParam{str2double(BLINKERParam{1})}));
end

% CLEANLINE PARAMETERS 
fprintf(fid,'\r\n\r\n%s','CLEANLINE parameters ----');
fprintf(fid,'\r\n%s',sprintf('- Line noise frequencies to remove : %s',AnswerCleanLine{1}));
fprintf(fid,'\r\n%s',sprintf('- P-value for detection of significant sinusoid : %s',AnswerCleanLine{2}));
fprintf(fid,'\r\n%s',sprintf('- Bandwidth [Hz] : %s',AnswerCleanLine{3}));
fprintf(fid,'\r\n%s',sprintf('- Sliding window length (sec) : %s',AnswerCleanLine{4}));
fprintf(fid,'\r\n%s',sprintf('- Sliding window step size (sec) : %s',AnswerCleanLine{5}));
fprintf(fid,'\r\n%s',sprintf('- Window overlap smoothing factor : %s',AnswerCleanLine{6}));
fprintf(fid,'\r\n%s',sprintf('- FFT padding factor : %s',AnswerCleanLine{7}));
            

% Summary of analyses ----------
if Analysis
    fprintf(fid,'\r\n\r\n%s\r\n','------ ANALYSIS SUMMARY ------');
    fprintf(fid,'\r\n\r\n%s\r\n','1) Single-subject level ------');
    fprintf(fid,'\r\n%s',sprintf('1) Compute Global Power Spectra (GPS) : %s',AnalysesSwitch{1,end}));
    fprintf(fid,'\r\n%s',sprintf('2) Compute GPS on scalp areas : %s',AnalysesSwitch{2,end}));
    if strcmpi(AnalysesSwitch{1,end},'Yes')
        fprintf(fid,'\r\n%s',sprintf('Plotted frequencies : %sHz',num2str(FreqsToPlot)));
        fprintf(fid,'\r\n\r\n%s','- Frequency Ranges -');
        for k=1:size(FreqData,1)
            fprintf(fid,'\r\n%s',sprintf('%d) %s : %sHz',k,FreqData{k,1},FreqData{k,2}));
        end
    end
    
    if strcmpi(Steps,'Both') || strcmpi(Steps,'Group analyses') 
        % STUDY analyses
        fprintf(fid,'\r\n\r\n%s\r\n','2) Groups/Conditions levels ------');
        fprintf(fid,'\r\n%s\r\n','A) MicroStates analyses ------');
        for p=1:length(MSTVersions)
            CurrentV = regexp(MSTVersions{p},'\d*','Match');
            CurrentStruct = LogMST.(MSTVersions{p});
            InternalFields = fieldnames(CurrentStruct);
            fprintf(fid,'\r\n%s',sprintf('VERSION %s ------',CurrentV{:}));
            for m=1:length(InternalFields)
                % Changing data to string
                if ~isstring(CurrentStruct.(InternalFields{m}))
                    TempData = num2str(CurrentStruct.(InternalFields{m}));
                else
                    TempData =CurrentStruct.(InternalFields{m});
                end
                fprintf(fid,'\r\n%s',sprintf('%d) %s : %s',...
                    m,InternalFields{m},TempData));
            end
        end
        
        % WRITE ALSO FOR IC CLUSTERING !! 
        LogIC
    end

end

% FILES PROCESSED ----------
fprintf(fid,'\r\n\r\n%s\r\n','------ PROCESSED FILES SUMMARY ------');
GroupsTemp=fieldnames(Participant_load);
% For each BS levels
% WHAT ABOUT MULTIPLE FACTORS ???? 
for k=1:numel(GroupsTemp)
    SubjTemp = fieldnames(Participant_load.(GroupsTemp{k}));
    fprintf(fid,'\r\n\r\n%s',sprintf('- %s -',upper(GroupsTemp{k})));
    % For each subjects
    for l=1:length(SubjTemp)
        % For each processed file
        FilesTemp = '';
        for t=1:length(Participant_load.(GroupsTemp{k}).(SubjTemp{l}).FileList)
            if t~=length(Participant_load.(GroupsTemp{k}).(SubjTemp{l}).FileList)
                FilesTemp = strcat(FilesTemp,Participant_load.(GroupsTemp{k}).(SubjTemp{l}).FileList{t},', ');
            else
                FilesTemp = strcat(FilesTemp,Participant_load.(GroupsTemp{k}).(SubjTemp{l}).FileList{t});
            end
        end
        fprintf(fid,'\r\n%s',sprintf('%s) : %s',upper(SubjTemp{l}),FilesTemp));
    end
end

% ERRORS ----------
% IMPLEMENT ERRORS LOG IN SCRIPT! e.g. files that cannot be processed for
% unknown reasons, maybe list potential reasons (e.g. not enough epochs,
% etc)
fprintf(fid,'\r\n\r\n%s\r\n','------ ERRORS ------');
fprintf(fid,'\r\n%s\r\n','------ 1) Preprocessing: ------');
if exist('ErrorArtifacts','var') && ~isempty(ErrorArtifacts)
    for k=1:size(ErrorArtifacts,1)
        fprintf(fid,'\r\n\r\n%s',sprintf(['File %s could not be processed due to massive artifacts'...
            newline 'We recommend that you manually remove the time period containing the artifacts'],...
            ErrorArtifacts{k}));
    end
else
    fprintf(fid,'\r\n No errors reported, congrats!');
end
fprintf(fid,'\r\n%s\r\n','------ 2) Group Analyses: ------');
if ~isempty('ErrorSTUDY')
     fprintf(fid,'\r\n%s',sprintf(['The following files in the STUDY displayed a different' ...
         'sampling rate than others, hence we downsampled all of them to %dHz:'],ErrorSTUDY{1}));
    for k=2:size(ErrorSTUDY,1) 
        fprintf(fid,'\r\n\r\n%s',sprintf('%d) %s',k-1,ErrorSTUDY{k}));
    end
else
    fprintf(fid,'\r\n No errors reported, congrats!');
end

% Closing the file
fclose(fid);

end