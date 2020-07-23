%%------------------------------RESTINGLAB-------------------------------%%

% Version 0.62
% Developped by <Corentin Wicht>
% 23.10.2019
% Author: Corentin Wicht (corentin.wicht@unifr.ch)
% Contributor: Christian Mancini (christian.mancini@unifr.ch)
%-------------------------------------------------------------------------%
% The script allows the user to perform semi-automatic artifact rejection 
% including ICA on BioSemi resting-state EEG recordings. If ICA is computed
% the script runs a second time through all files to enable the user to
% reject components which will be pre-selected with dedicated algorithms.

                              % ENJOY %

% DEPENDENCIES:
% EEGLAB V. XXX (REF)
% MORE


% LICENSE !!!
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%      MAIN SCRIPT FOR PREPROCESSING AND ANALYSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function MAIN(App,CurrentPWD,Date_Start,SavePath)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       GLOBAL VARIABLES DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get time
time_start = datestr(now);

% Hide warning messages
warning off MATLAB:subscripting:noSubscriptsSpecified

% Adding path to dependencies
addpath([CurrentPWD '\Functions\eeglab14_1_2b']);
addpath([CurrentPWD '\Functions\']);
addpath(genpath([CurrentPWD '\Functions\Dependencies']));
% addpath(genpath([CurrentPWD '\Functions\eeglab14_1_2b\plugins\blinkerv1.1.2']));
addpath([CurrentPWD '\Functions\EEGInterp']); 

%-------------------------------------------------------------------------%
% LOADING GUI MATRICES
%-------------------------------------------------------------------------%
% Retrieving the list of all matrices
ParametersPath=dir([SavePath '\Parameters\' Date_Start '\']);
ParametersPath=ParametersPath(~cell2mat({ParametersPath.isdir}));

% Loading the matrices stored with the GUI
for k=1:length(ParametersPath)
    load([ParametersPath(k).folder '\' ParametersPath(k).name])
end

%-------------------------------------------------------------------------%
% INSTALL MPICH2 for AMICA plugin
%-------------------------------------------------------------------------% 

% Installing additionnal softwares for AMICA
ComputerType = computer;

% Determining the number and letters of harddrives
F = getdrives('-nofloppy');

% Finding out whether it's already installed or not
Result = 0;
for k=1:length(F)
    if exist([F{k} 'Program Files\MPICH2'],'dir') || ...
        exist([F{k} 'Program Files (x86)\MPICH2'],'dir') 
        Result = 1;
    end
end

if ispc && ~Result && strcmpi(ICA,'yes')
    
    % Temporary folder to delete
    mkdir([CurrentPWD '\Mpich2_1.4']);
    
    % AMICA
    % https://sccn.ucsd.edu/~jason/amica_web.html
    % Saving data from web & installing
    if str2double(cell2mat(regexp(ComputerType,'\d*','Match')))==64 
        if ~exist([CurrentPWD '\Mpich2_1.4\mpich2-1.4-win-x86-64.msi'],'file')       
            urlwrite('http://www.mpich.org/static/downloads/1.4/mpich2-1.4-win-x86-64.msi',...
                [CurrentPWD '\Mpich2_1.4\mpich2-1.4-win-x86-64.msi']);
        end
        system([CurrentPWD '"\Mpich2_1.4\mpich2-1.4-win-x86-64.msi"']);
    elseif str2double(cell2mat(regexp(ComputerType,'\d*','Match')))==32 
        if ~exist([CurrentPWD '\Mpich2_1.4\mpich2-1.4-win-ia32.msi'],'file')   
            urlwrite('http://www.mpich.org/static/downloads/1.4/mpich2-1.4-win-ia32.msi',...
                [CurrentPWD '\Mpich2_1.4\mpich2-1.4-win-ia32.msi'])
        end
        system([CurrentPWD '\Mpich2_1.4\mpich2-1.4-win-ia32.msi'])
    end
 
    % Delete temporary directory
    cd(CurrentPWD)
    % ERROR HERE !! 
    % rmdir('Mpich2_1.4','s')
    
end

%% DESIGN / EEG PARAMETERS
% Between-subject factor(s)
BetweenFactors=BetweenFactors(~cellfun('isempty',BetweenFactors));
if ~isempty(BetweenFactors)
    for k=1:length(BetweenFactors)
        Groups_Names(:,k) = strrep(BetweenLevels(k,~cellfun('isempty',BetweenLevels(k,:))),' ','');
    end
else
    Groups_Names = strsplit(FilesPath{1},'\') ;
    Groups_Names = Groups_Names(end);
end

% Within-subject factor(s)
WithinFactors=WithinFactors(~cellfun('isempty',WithinFactors));
if ~isempty(WithinFactors)
    for k=1:length(WithinFactors)
        Conditions_Names(:,k) = WithinLevels(k,~cellfun('isempty',WithinLevels(k,:)));
        Conditions=1:length(Conditions_Names);
    end
else
    Conditions = 1;
    % Conditions_Names =Session1; % MIGHT CAUSE AN ERROR !!!!
end

% EEG data
Channels=str2num(Channels);
FreqNames = FreqData(~cellfun('isempty',FreqData(:,1)),1);
% FreqRanges = FreqData(~cellfun('isempty',FreqData(:,2)),2);
% FreqRanges = cellfun(@str2num,FreqRanges,'UniformOutput',false);
% FreqRanges = cell2mat(FreqRanges);
FreqRanges = cell2mat(cellfun(@(x) str2num(x), FreqData(:,2),'UniformOutput',0));
Frequency_export = cellfun(@(x) ['S%d_' x '_'],FreqNames,'UniformOutput',false);
if max(Channels)<128
    DirectoryTemp = dir(['**/*' '.locs']);
    FolderTemp = DirectoryTemp(contains({DirectoryTemp.folder},'ChanLocs')).folder;
    Channel_load = [FolderTemp '\biosemi64.locs'];
else
    DirectoryTemp = dir(['**/*' '.xyz']);
    FolderTemp = DirectoryTemp(contains({DirectoryTemp.folder},'ChanLocs')).folder;
    Channel_load=[FolderTemp '\biosemi128.xyz'];
end

% Plugins settings (if user pressed cancel the settings will be empty)
if isempty(AnswerCleanLine)
   CleanLineParam = {'50 100','0.01','2','4','4','100','2'}; % Default 
end
if isempty(AnswerASR)
   AnswerASR = {'10',''}; % Default 
else
    if isempty(AnswerASR{1})
        AnswerASR{1} = '-1';
    end
    if isempty(AnswerASR{2})
        AnswerASR{2} = '-1';
    end
end
if strcmpi(ICA,'Yes') && isempty(AnswerICA)
   AnswerICA = {'2000'}; % Default 
end
if str2double(AnswerBLINKER{1}) == 1
    AnswerBLINKER = 'reject';
else
    AnswerBLINKER = 'interp';
end

% Decision if doing basic analyses
if exist('AnalysesSwitch','var')
    Analysis = 1;
end

% List of electrode groups
if strcmpi(AnalysesSwitch{2,end},'Yes')
    AreasList = AreasList(cellfun(@(x) ~isempty(x),AreasList));
    SplitHeaders = strsplit(AreasList{1},' ');
    for k = 2:length(AreasList)
        if ~isempty(str2num(AreasList{k})) % If user entered numbers
            SplitChans=str2num(AreasList{k});
            NewAreasList.(SplitHeaders{k-1})=SplitChans';
        else
            Delimiter = '\t';FormatSpec = '%q%q%q%[^\n\r]';
            FileID = fopen(Channel_load,'r');
            ChanList = textscan(FileID, FormatSpec, 'Delimiter', Delimiter,'EndOfLine', '\r\n');
            ChanList = ChanList{end};
            fclose(FileID);
            SplitChans = strsplit(AreasList{k},' ');
            for t=1:length(SplitChans)
                Idx = find(ismember(lower(ChanList),lower(SplitChans{t})));
                NewAreasList.(SplitHeaders{k-1})(t)=Idx;
            end
        end
    end
end
%% DIRECTORIES 
% If save path different than CurrentPWD
if sum(~SavePath)<1
    
    % Excel directory
%     mkdir([SavePath '\Excel\' Date_Start])
    ExcelDirectory=[SavePath '\Excel\' Date_Start '\'];
    
    % Exports directory
    StatsDirectory=[SavePath '\Exports\' Date_Start '\SourceLocalisation\'];
    TopoplotsDirectory=[SavePath '\Exports\' Date_Start '\Topoplots\'];
    
    % Finds out if all folder for source localisation
    StatsFolders=dir(StatsDirectory);
    StatsFolders=StatsFolders([StatsFolders.isdir]);
    for k=1:length(FreqNames)
       if ~ismember(FreqNames{k},{StatsFolders.name}) 
           % Create missing directories
           mkdir([StatsDirectory FreqNames{k}])
       end
    end
    
    % Additionnal directories
    TopoDipFitDirectory=[SavePath '\Exports\'  Date_Start '\TopoplotsDipFit\'];
    
else
    % If save path same as CurrentPWD
    % Excel directory
    ExcelDirectory=[pwd '\Excel\' Date_Start '\'];

    % Exports directory
    StatsDirectory=[pwd '\Exports\SourceLocalisation\' Date_Start '\'];
    TopoplotsDirectory=[pwd '\Exports\Topoplots\' Date_Start '\'];

    % Finds out if all folder for source localisation
    StatsFolders=dir(StatsDirectory);
    StatsFolders=StatsFolders([StatsFolders.isdir]);
    for k=1:length(FreqNames)
       if ~ismember(FreqNames{k},{StatsFolders.name}) 
           % Create missing directories
           mkdir([StatsDirectory FreqNames{k}])
       end
    end

    % Additionnal directories
    TopoDipFitDirectory=[pwd '\TopoplotsDipFit\' Date_Start '\'];
end

% Datasets templates (Adding the subject number before)
Dataset_filtCleaned = ['%s' FileNames '%d_filtered'];
Dataset_interp = ['%s' FileNames  '%d_filtered_cleaned'];
Dataset_filtCleaned_ICAed = ['%s' FileNames '%d_filtered_cleaned_ICAed'];
Dataset_filtCleaned_ICAedRejected = ['%s' FileNames '%d_filtered_cleaned_ICAedRejected'];
PreprocessedEEG=['%s' FileNames '%d_Preprocessed.bdf'];
PSDEEG=['%s' FileNames '%d_PowerSpectDensity.mat'];

%% SUBJECTS TEMPLATES
Conditions_Order=readtable([DirectoryCond FileCond]);
Conditions_OrderCell=table2cell(Conditions_Order(:,2:end));
Conditions_Labels = Conditions_Order.Properties.VariableNames(2:end-1);
TempSubjectslist = table2cell(Conditions_Order(:,1));

% If only 1 condition
if Conditions<2
    Conditions_Names = {FileNames};
end

% Only keeping subjects selected for analysis
j=1;t=1;
% TempSubjectslist=cellfun(@(x) str2double(x),TempSubjectslist);
TempSubjectslist=cell2mat(TempSubjectslist); % Error since first line is integer and not string anymore ?!
for k=1:length(TempSubjectslist)
    for l=1:size(ProcessData,2) % Number of sessions  
        if ProcessData{k,l}
            Subjectslist(j) = TempSubjectslist(k);
            j=j+1;
            break
        end
        
        % Scan the ProcessData line to remplace empty fields by 0
        ProcessData(k,cellfun(@(x) isempty(x),ProcessData(k,:))) = {false};
        
        % This will detect files for which there is no data to analyze
        if sum(cell2mat(ProcessData(k,:)))<1
%             SubjectsRemoved(t) = TempSubjectslist(k);
            SubjectsRemovedIdx(t) = k;
            t=t+1;
            break
        end
    end
end

% Removing lines from Conditions_Order (if completely excluded)
if exist('SubjectsRemovedIdx','var')
    Conditions_OrderCell(SubjectsRemovedIdx,:)=[];
end

% Building the subjects selection list based on GUI data
SubjectsIncluded = [num2cell(TempSubjectslist) ProcessData];
FilesToProcess = 0;
% SubjectsIncluded = cellfun(@(x) num2str(x), SubjectsIncluded,'UniformOutput',false);

% No idea where it comes from but sometimes it is loaded
if exist('Participant_load','var')
   clear Participant_load 
end

% Participants group directory
for k=1:length(FilesPath)
    
    % Path of the folders containing the datasets for each group
    FileList = dir([FilesPath{k} '\**/*' lower(Extension)]); 
    
    % Removing unnecessary files:
    % 1) Removing Preprocessed.bdf files if preprocessing selected
    FileList = FileList(~contains({FileList.name},'Preprocessed'));
    
    % 2) If loading .set, remove unused file versions
    if strcmpi(Extension,'.set')
        % non-ICA files
        IdxNonICA = contains({FileList.name}','filtered_cleaned');
        
        % ICA files
        IdxICA = contains({FileList.name}','filtered_cleaned_ICAedRejected');
        
        % Grouping indexes and removing unused files
        IdxRem = or(IdxNonICA,IdxICA);
        FileList = FileList(IdxRem);
    end
    
    % If save path different than CurrentPWD
    if sum(~SavePath)<1
        % Creating folder templates for each subject
        SplitTemp =  strsplit(FilesPath{k},'\');
        FoldPath = SplitTemp{end};
        mkdir(SavePath,FoldPath);
        for p=1:length(FileList)
            SplitTemp=strsplit(FileList(p).folder,'\');
            mkdir([SavePath '\' FoldPath],SplitTemp{end})
        end
    end

    % Only keeping the ones containg the common name
%     if strcmpi(Extension,'.bdf')
        % If file extension is .bdf 
        
    % Works for .bdf! needs to be tested for .set !!! Might be a problem 
    % due to more than 1 file per data    
    FileList = FileList(contains({FileList.name},FileNames));
    
%         p=1;
%         for m=1:length(FileList)
%             if sum(ismember(FileList(m).name,FileNames))<=length(FileNames)
%                 NewFileList(p,1) = FileList(m);
%                 p=p+1;
%             end
%         end
%         FileList = NewFileList;
%     else
%         % If file extension is .set
%         FilesLength = cellfun(@(x) length(x), {FileList.name});
%         [~,SmallPos] = mink(FilesLength,length(Conditions)*length(Groups_Names)*length(Subjectslist));
%         FileList = FileList(SmallPos);
%     end
    
    % List of all folders in the group
    GroupFoldersTemp = {FileList(contains({FileList.folder},FilesPath{k})).folder};
    UniqueFoldersTemp = sort_nat(unique(GroupFoldersTemp));
    
    for j=1:length(UniqueFoldersTemp)
        SplitFileTemp = strsplit(UniqueFoldersTemp{j},'\');
        % UniqueNamesTemp = unique({FileList.name});
        
        % Position of files
        Pos = 1;
        
        % Current subject
        CurrentSubj = str2num(cell2mat(regexp(SplitFileTemp{end},'\d*','Match')));
        
        % Number of files for the current subject
        NumofFiles = find(strcmpi({FileList.folder},UniqueFoldersTemp{j}));
        
        % For each unique file (for each participant/folder) 
        for l=1:length(NumofFiles)
            CurrentFile = {FileList(NumofFiles(l)).name};
            CurrentFileSplit = strsplit(CurrentFile{:},'.');
            
            % Position in excel files 
            TempFilePos = cellfun(@(x) x==CurrentSubj, SubjectsIncluded(:,1)); % Subject
%             if length(Conditions) > 1 % Session
%                 CondPos = find(contains(upper(Conditions_Labels),upper(CurrentFileSplit{1})));
%             else
%                 CondPos = 1;
%             end
            
            % Preallocating array
            Condname_i = zeros([size(Conditions_Labels,2) 1]);

            % For each condition, see if we find it in the name or subpath
            if length(Conditions) > 1
                for f = 1:size(Conditions_Labels,2)
                    Condname_i(f) = contains(upper(CurrentFileSplit{1}),upper(Conditions_Labels{f}));
                end
                
                % In case of multiple positives, take the lengthier condition name
                if sum(Condname_i) > 1
                    [~,CondPos] = max(cellfun('length',Conditions_Labels) .* Condname_i);
                else 
                    CondPos = find(Condname_i == 1);
                end
            else
                CondPos = 1;
            end
            
            % Only including data that was selected in the GUI
            if cell2mat(SubjectsIncluded(TempFilePos,CondPos+1)) == 1 % initially str2double 
                
                % HERE WRITE FOR THE LOG ! 
                Participant_load.(Groups_Names{k}).(SplitFileTemp{end}).Path = UniqueFoldersTemp(j);
                
                % Create the folder list content structure called FileList
                Participant_load.(Groups_Names{k}).(SplitFileTemp{end}).FileList(Pos) = CurrentFile;
                
                % Retrieving the conditions assignement values
                if length(Conditions_Names)>1
                    Participant_load.(Groups_Names{k}).(SplitFileTemp{end}).CondAssign(Pos) = ...
                        Conditions_OrderCell(TempFilePos,l); 
                else
                    Participant_load.(Groups_Names{k}).(SplitFileTemp{end}).CondAssign(Pos) = ...
                        Conditions_Names(l); 
                end
                
                % Writting the export path
                if sum(~SavePath)<1
                    SplitFile = strsplit(UniqueFoldersTemp{j},'\');
                    splitFolder = strsplit(FilesPath{k},'\');
                    Participant_load.(Groups_Names{k}).(SplitFileTemp{end}).ExportPath =...
                        [SavePath '\' splitFolder{end} '\' SplitFile{end}]; % UNSURE WILL ALWAYS WORK!!!
                else
                    Participant_load.(Groups_Names{k}).(SplitFileTemp{end}).ExportPath =  UniqueFoldersTemp(j);
                end
                
                % Count the number of files to process below
                FilesToProcess = FilesToProcess + 1; Pos = Pos + 1;
            end
        end
    end
end

%% Excel templates

try 
    % If save path different than CurrentPWD
    if sum(~SavePath)<1
        ExcelPath = [SavePath '\Excel\' Date_Start];
        ExcelFiles=dir([ExcelPath  '\**/*' '.xlsx']);

        % Creating excel templates if they do not exist
        if isempty(ExcelFiles) % sum(contains({ExcelFiles.name}, 'AreaAmplitude'))==0 && sum(contains({ExcelFiles.name}, 'GPS'))==0
            CreateTemplates(num2cell(TempSubjectslist),Conditions_Names,FreqNames,Channels,[SavePath '\Excel\' Date_Start '\'],ExcelFiles);
        end
    else
        Temp=what('Excel');
        ExcelPath = [Temp.path '\' Date_Start];
        ExcelFiles=dir([ExcelPath  '\**/*' '.xlsx']);

        % Creating excel templates if they do not exist
        if isempty(ExcelFiles) % sum(contains({ExcelFiles.name}, 'AreaAmplitude'))==0 && sum(contains({ExcelFiles.name}, 'GPS'))==0
            CreateTemplates(num2cell(TempSubjectslist),Conditions_Names,FreqNames,Channels,[Temp.path '\' Date_Start '\'],ExcelFiles);
        end
    end
catch
    warning(['Since excel templates already exist and you just changed the FileName, the script crashed.',...
        newline 'To avoid this, please remove all previously generated templates from the following folder and re-run the script:', ...
        newline ExcelPath])
    
end

% LOADING THE TEMPLATES
SleepNoSleepTable=readtable(backslash([ExcelDirectory 'AsleepAwakeTrials.xlsx']));

for m=1:size(FreqRanges,1)

    % Importing areas amplitude templates
    AreaAmpTable.(FreqNames{m})=readtable(backslash([ExcelDirectory ['AreaAmplitude' FreqNames{m} '.xlsx']]));

    for n=1:length(Conditions)
      % Temporary sheet name
        TempSheetsName = [FreqNames{m} '_' lower(Conditions_Names{n})];

        % sheets name longer than or equal to 31 characters will crash
        if length(TempSheetsName) >= 31
            TempSheetsName = TempSheetsName(1:30);
        end
        
       % Importing Global Power Spectra templates
        GPS.(FreqNames{m}).(Conditions_Names{n})=readtable(backslash([ExcelDirectory ['GPS_' FreqNames{m} '.xlsx']]),...
            'Sheet',TempSheetsName);  
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epitome of UI
WaitBarApp = uiprogressdlg(App.MainGUIUIFigure,'Title','Progress Bar',...
    'Message','','Cancelable','on');
File = 0;

% If the extension is .set, only the analyses are performed
if ~strcmpi(Extension,'.set') && strcmpi(Steps,'Preprocessing') || strcmpi(Steps,'Both')
    
    % Participants loop
    for g=1:length(Subjectslist)

        % Current Subject
        ParticipantNumber = Subjectslist(g);

        % Position of the current subject in the subject list
        Pos = find(Subjectslist == ParticipantNumber);

        % Retrieving current subject information
        % WILL NOT WORK WITH 2 BETWEEN-SUBJ FACTORS !! 
        
        FolderTemplate = regexp(SplitFileTemp{end},'\D*','Match');
        if length(FilesPath)>1
            CurrentFile = Participant_load.(Conditions_OrderCell{Pos,end}). ...  
            ([FolderTemplate{:} num2str(ParticipantNumber)]);
        else
            CurrentFile = Participant_load.(Groups_Names{1}). ...  
            ([FolderTemplate{:} num2str(ParticipantNumber)]);
        end

        % Conditions loop
        for h=1:length(CurrentFile.FileList)
            %-----------------------------------------------------------------%    
            % Loading directory and templates
            %-----------------------------------------------------------------% 

            % Finding current condition number
            WhichCond = find(contains(lower(Conditions_Names),...
                  lower(CurrentFile.CondAssign{h})));

            % Selecting appropriate groups folder     
            Dir_load = backslash(CurrentFile.Path{:});
            Subj_load = [Dir_load backslash(['\' CurrentFile.FileList{h}])];
            
            % If save path different than CurrentPWD
            if ~isempty(SavePath)
                Dir_save = backslash(CurrentFile.ExportPath);
            else
                Dir_save = Dir_load;
            end     
                    
            % Update progress, report current estimate
            File = File + 1;
            WaitBarApp.Value = File/FilesToProcess;
            WaitBarApp.Title = '1. PREPROCESSING: Loading data';
            WaitBarApp.Message = sprintf('Sbj %d/%d : %s',g,...
                length(Subjectslist),CurrentFile.CondAssign{WhichCond});

            % Check for Cancel button press
            if WaitBarApp.CancelRequested
                return
            end
            
            
            %% IMPORT
            % Creating first dataset
            [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
            close gcf

            % set double-precision parameter
            pop_editoptions('option_single', 0);

            % Import the .bdf file
            % SHOULD NEVER SPECIFY A REF CHAN FOR BIOSEMI (since BioSemi uses CMS-DRL which cannot be imported)
            % https://sccn.ucsd.edu/pipermail/eeglablist/2014/007822.html
            % https://blricrex.hypotheses.org/files/2015/03/Pr%C3%A9sentationCREx_EEGLAB-corrig%C3%A9.pdf
            EEG = pop_biosig(Subj_load, 'channels', Channels); 

            % Loading BioSemi channel location
            EEG=pop_chanedit(EEG, 'load',{Channel_load 'filetype' 'autodetect'});

            %% RESTRICTING DATA LENGTH (optional)
            if End~=0
                EEG = pop_select( EEG,'time',[Beginning ...
                    End]);  
            end
            
            %% Changing events structure (double to strings)
            EventFields={'event','urevent'};
            for t = 1:length(EventFields)
                if ~isempty(EEG.(EventFields{t}))
                    TempEvents = {EEG.(EventFields{t}).type};
                    for f=1:length(TempEvents)
                        if isa(TempEvents{f},'double') 
                            EEG.(EventFields{t})(f).type = num2str(EEG.(EventFields{t})(f).type);
                        end
                    end
                end
            end
            %% FILTERING
            
            % Waitbar updating 
            WaitBarApp.Title = '1. PREPROCESSING: Filtering';

            % Sampling rate
            if EEG.srate>SamplingRate
                EEG = pop_resample(EEG, SamplingRate);
            end
            
            % High-pass/Low-pass filter
            EEG = pop_eegfiltnew(EEG,HighPass,LowPass); 
           
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %                    ARTIFACTS REJECTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
            % Waitbar updating 
            WaitBarApp.Title = '1. PREPROCESSING: CleanLine';

            % CLEANLINE
            % CleanLine sinusoidal stationary noise removal
            % (http://www.antillipsi.net/research/software#TOC-Cleanline)
            PromptCleanLine=cellfun(@(x) str2num(x),AnswerCleanLine,'UniformOutput',false);
            EEG = pop_cleanline(EEG,'bandwidth',PromptCleanLine{3},'chanlist',1:EEG.nbchan ,...
                'computepower',0,'linefreqs',PromptCleanLine{1},'normSpectrum',0,'p',PromptCleanLine{2},...
                'pad',PromptCleanLine{7},'plotfigures',0,'scanforlines',1,'sigtype','Channels',...
                'tau',PromptCleanLine{6},'verb',1,'winsize',PromptCleanLine{5},'winstep',PromptCleanLine{4});
            close gcf;
                 
            % Finding current subject export name
%             CurrentSession = CurrentFile.FileList{h};
%             CurrentSession = strsplit(CurrentSession,'.');
%             SubjName = strtok(CurrentSession{1},FileNames);
%             % If there is no file specific name
%             if sum(isstrprop(SubjName,'digit'))>0
                Temp = strsplit(CurrentFile.Path{:},'\');
                SubjName = [Temp{end} '_'];
%             end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%                    INTERPOLATION/CHANNELS REJECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%n%%%%%%%%%%%%             

            % Temporary save the original EEG (non-robust averaged)
            OriginalEEG = EEG;
            
            % Waitbar updating 
            WaitBarApp.Title = '1. PREPROCESSING: Channels interpolation';
            
            % Detection of bad channels 
            ErrorArtifacts = {};
            r = 1;
            try
                % Function from the PrepPipeline to perform 1) Bad Channels
                % Detection and rejection
%                 [~, InterpChanStruct] =  performReference(EEG);
                     
                % Parameters
                InterpChanStruct = getReferenceStructure();
                defaults = getPrepDefaults(EEG, 'reference');
                InterpChanStruct = checkDefaults(struct(), InterpChanStruct, defaults);
                defaults = getPrepDefaults(EEG, 'detrend');
                InterpChanStruct = checkDefaults(struct(), InterpChanStruct, defaults);
                InterpChanStruct.rereferencedChannels = sort(InterpChanStruct.rereferencedChannels);
                InterpChanStruct.referenceChannels = sort(InterpChanStruct.referenceChannels);
                InterpChanStruct.evaluationChannels = sort(InterpChanStruct.evaluationChannels);
                
                % Bad channels identification by robust reference
                InterpChanStruct.referenceSignalOriginal = ...
                    nanmean(EEG.data(InterpChanStruct.referenceChannels, :), 1);
                InterpChanStruct = robustReference(EEG, InterpChanStruct);
                
                % Saving the bad channels data to reintroduce them later
                % Saving the bad channels data to reintroduce them later
                EEG.BadChans.chanlocs = EEG.chanlocs; 
                EEG.BadChans.nbchan = EEG.nbchan;
                EEG.BadChans.data = EEG.data(InterpChanStruct.badChannels.all,:);
                EEG.BadChans.InterpChans = InterpChanStruct.badChannels.all;
                
                % Removing the bad channels
                EEG = pop_select(EEG,'nochannel',InterpChanStruct.badChannels.all);
            catch
                % Write error to the LOG
                ErrorArtifacts{r} = [SubjName FileNames num2str(WhichCond)];
                r = r+1;
                break
            end
            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%                    NON-SINUSOIDAL NOISE REMOVAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%n%%%%%%%%%%%%  

            % Saving data for visual comparision below
            OriginalEEG = EEG;
           
            % Waitbar updating 
            WaitBarApp.Title = '1. PREPROCESSING: ASR';
            
            % ARTIFACT SUBSPACE RECONSTRUCTION (ASR)
            % Automated bad channels detection and non-stationary noise removal
            % Christian's method (http://sccn.ucsd.edu/eeglab/plugins/ASR.pdf)
            % Non-stationary artifacts removal
            
            % Fixing the maximum available memory enhances reproducibility
            MaxMemory = round(hlp_memfree()/2000000,-3);
 
            if strcmpi(ASR,'Yes') 
                
                % ASR 
                EEG = clean_rawdata(EEG, -1, -1, -1, -1, str2double(AnswerASR{1}),...
                    str2double(AnswerASR{2}));          

                % Removing the events that were created by ASR 
                EEG = RemovEvents(EEG,'EventTypes',{'X','boundary'});
                
                % Allows for visual inspection of old/new EEG (Optional)
                if strcmpi(Automaticity,'Yes')
                    vis_artifacts(EEG,OriginalEEG);
                    disp('Ignore errors above !!');
                    % Holds the figure until inspection is over
                    Fig=msgbox('THE CODE WILL CONTINUE ONCE YOU PRESS OK','WAIT','warn'); 
                    uiwait(Fig);
                    close gcf
                end 
            end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%             EYE BLINKS ARTIFACTS REJECTION/INTERPOLATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%n%%%%%%%%%%%% 
            
            % BLINKER: EYE BLINKS REMOVAL
            if strcmpi(BLINKER,'Yes')
                
                % Waitbar updating 
                WaitBarApp.Title = '1. PREPROCESSING: BLINKER';

                % Setting parameters
                Params = checkBlinkerDefaults(struct(), getBlinkerDefaults(EEG));
                Params.fileName = [Dir_save '\\' sprintf(Dataset_filtCleaned, SubjName, WhichCond)];
                Params.blinkerSaveFile = [Dir_save '\\' sprintf(Dataset_filtCleaned, SubjName, WhichCond) '_blinks.mat'];
                Params.showMaxDistribution = true;
                Params.verbose = false;
                Params.fieldList = {'leftBase','rightBase'};
                
                % Run BLINKER algorithm
                try
                    [EEG, ~, blinks, blinkFits,blinkProperties] = pop_blinker(EEG, Params);
                    close gcf;fprintf('%d blinks identified by BLINKER\n',length(blinkFits));
                catch
                    disp('0 blink identified by BLINKER\n');
                end
                
                % Add the blinks to EEG.event
%                 [EEG, ~] = addBlinkEvents(EEG, blinks, blinkFits,...
%                     blinkProperties, Params.fieldList);
%                 eegplot(EEG.data,'events',EEG.event,'winlength',30)

                % Retrieve blinks latency
                for m=1:length(blinkFits)
                    AbsLatency(m,:) = [blinkFits(m).leftBase blinkFits(m).rightBase];
                end
                
                if strcmpi(AnswerBLINKER,'reject') && ~isempty(blinkFits)
                    
                    % Removing the data containing blinks
                    EEG = eeg_eegrej(EEG, AbsLatency);
                end
                          
                % Removing the events that were created by eeg_eegrej()
                EEG = RemovEvents(EEG,'EventTypes',{'boundary'});
            end

            if strcmpi(BLINKER,'Yes') && strcmpi(AnswerBLINKER,'interp') ...
                    && ~isempty(blinkFits)
                
                % THIS IS CURRENTLY NOT WORKING!!! I CANNOT FORCE TO
                % INTERPOLATE ALL THE SIGNAL I AM PROVIDING IT!!!!
                % SHIT....
                
                % 1) ASR 1st run
                % ASR settings
                asr_windowlen = max(0.5,1.5*EEG.nbchan/EEG.srate);
                BurstCriterion = str2double(AnswerASR{1});asr_stepsize = 4;
                maxdims = 1;usegpu = false;

                % Creating a clean reference section
                EEGCleanRef = clean_windows(EEG,0.075,[-3.5 5.5],1); 

               %  Building the Blink EEG dataset
                BlinkEEG = EEG; TempEEG = EEG;
                BlinkEEG.data = []; BlinkEEG.event = []; Pos = 1;
                for p=1:size(AbsLatency,1)
%                     BlinkEEG.data(:,AbsLatency(p,1):AbsLatency(p,2)) = ...
%                         EEG.data(:,AbsLatency(p,1):AbsLatency(p,2));
                    BlinkEEG.data = [BlinkEEG.data EEG.data(:,AbsLatency(p,1):AbsLatency(p,2))];
                    BlinkEEG.event(Pos).type = 'leftBase';
                    if p>1
                        BlinkEEG.event(Pos).latency = BlinkEEG.event(Pos-1).latency+1;
                    else
                        BlinkEEG.event(Pos).latency = 1;
                    end
                    Pos = Pos + 1;
                    BlinkEEG.event(Pos).type = 'rightBase';
                    BlinkEEG.event(Pos).latency = BlinkEEG.event(Pos-1).latency + (AbsLatency(p,2)-AbsLatency(p,1));
                    Pos = Pos + 1;
                end
            
                % Calibrate on the reference data
                state = asr_calibrate_r(EEGCleanRef.data, EEGCleanRef.srate,...
                    BurstCriterion, [], [], [], [], [], [], []);

                % Extrapolate last few samples of the signal
                sig = [BlinkEEG.data bsxfun(@minus,2*BlinkEEG.data(:,end),...
                    BlinkEEG.data(:,(end-1):-1:end-round(asr_windowlen/2*BlinkEEG.srate)))];
                
                % Process signal using ASR
                [BlinkEEG.data,state] = asr_process_r(sig,BlinkEEG.srate,state,...
                    asr_windowlen,asr_windowlen/2,asr_stepsize,maxdims,MaxMemory,usegpu);
                
                % Shift signal content back (to compensate for processing delay)
                BlinkEEG.data(:,1:size(state.carry,2)) = [];

                % Replace original data with the Blink corrected data
                NEWEEG = EEG; Pos = 1:2:size(AbsLatency,1)*2;
                for p=1:size(AbsLatency,1)
                    NEWEEG.data(:,AbsLatency(p,1):AbsLatency(p,2)) = ...
                        BlinkEEG.data(:,BlinkEEG.event(Pos(p)).latency:...
                        BlinkEEG.event(Pos(p)+1).latency);
                end

                % Allows for visual inspection of old/new EEG (Optional)
                if strcmpi(Automaticity,'Yes')
                    vis_artifacts(NEWEEG,EEG);
                    disp('Ignore errors above !!');
                    % Holds the figure until inspection is over
                    Fig=msgbox('THE CODE WILL CONTINUE ONCE YOU PRESS OK','WAIT','warn'); 
                    uiwait(Fig);
                    close gcf
                end 
                
                % Replace artifact data
                EEG.data = NEWEEG.data;
            end
            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%                    SLEEP / NO SLEEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Here we implemented a method to determine segments of
            % recording where the participant was awake or asleep.
            % The analysis is based on Diaz et al. (2016).
            % Overall higher amplitude of alpha band compared to theta band
            % would mean that the participant is still awake. If the opposite
            % is true, then the epoch will be classified as sleep.

           if strcmpi(Sleep, 'Yes') 
                
                % TO CONTINUE:                
                
                % TESTING THE NEW ALGORITHM:
                % https://github.com/alexander-malafeev/feature-based-sleep-scoring  
                
%                 % Channel to use
%                 SleepChan = ismember({EEG.chanlocs.labels},'C1'); % GUI!!
%                 LOCChan = ismember({EEG.chanlocs.labels},'Fp1');
%                 ROCChan = ismember({EEG.chanlocs.labels},'Fp2');
%                 
%                 % Main electrode data
%                 EEGSleep = squeeze(EEG.data(SleepChan,:));
%                 
%                 % truncate signals, i.e. it should contain integer number of epochs
%                 EpochsTF = SecondsEpoch*EEG.srate;
%                 MaxExp=floor(size(EEG.data,2)/EpochsTF); % how many X second epochs are in our file
%                 EEGSleep = EEGSleep(1,1:MaxExp*EpochsTF);
%                 
%                 % Additional electrodes
%                 EOG = zeros(1,size(EEGSleep,2));
%                 EMG = zeros(1,size(EEGSleep,2));
%                 LOC = squeeze(EEG.data(LOCChan,:));
%                 LOC = LOC(1,1:MaxExp*EpochsTF);
%                 ROC = squeeze(EEG.data(ROCChan,:));
%                 ROC = LOC(1,1:MaxExp*EpochsTF);
%                 
%                 % Channel power spectra
%                 FreqResol = 5; % GUI!!
%                 SamplingRate = EEG.srate;
%                 Spect = spectopo(EEGSleep,0,EEG.srate,'freqfac',FreqResol,'plot','off');
%                 
%                 % Extracting sleep-related features on 1 EEG channel
%                 [ X, FeaturesNames] = extractFeatures(Spect', EEGSleep,  EOG,  EMG,...
%                     SecondsEpoch, SamplingRate, FreqResol, LOC, ROC);

                %% 
                % Temporary epoching the file 
                EEG = eeg_regepochs(EEG, SecondsEpoch);

                % WAVELET MORLET DECOMPOSITION
                % composed or real and imaginery numbers : Gaussian and complex sine wave
                % see : https://jallen.faculty.arizona.edu/sites/jallen.faculty.arizona.edu/files/Chapter_13_Complex_Morlet_Wavelets_Power_Phase.pdf
                % and : http://www.mikexcohen.com/left_toc.html
                
                % Waitbar updating 
                WaitBarApp.Title = '1. PREPROCESSING: Sleep rejection';
                
                % Preallocation
                WTData.Morlet=[];
                SleepNoSleep=zeros(size(EEG.data,3),1); 

                % Looping over epochs
                for j=1:size(EEG.data,3)
                   % Computing Wavelet Morlet Decomposition
                   [WTData.Morlet(:,:,j), Freq]=wt(EEG.data(:,:,j),EEG.srate,'fmin',2.6,'fmax',...
                        48,'Preprocess','on','Wavelet','Morlet','Padding',0,'plot','off','Display','off','f0',0.1);
                end

                % Calculating Power of complex number (i.e. absolute value)
                WTData.Absolute.All=abs(WTData.Morlet);

                % Calculating EEG bandpass filtered (i.e. only taking real number)
                % WTData.Absolute.All=real(WTData.Morlet);
                % Calculating Phase of complex number (arctan(imag/real))
                % WTData.Absolute.All=arctan(imag(WTData.Morlet)/real(WTData.Morlet));

                % Looking for indexes for Theta and Alpha bands
                Index.Theta=find(Freq>=5 & Freq<=7);
                Index.Alpha=find(Freq>=8 & Freq<=12);

                % Computing average WTdata for Theta and Alpha indexes
                WTData.Absolute.Theta=mean(WTData.Absolute.All(Index.Theta,:,:),1); 
                WTData.Absolute.Alpha=mean(WTData.Absolute.All(Index.Alpha,:,:),1);
                WTData.Absolute.Theta=squeeze(mean(WTData.Absolute.Theta,2));
                WTData.Absolute.Alpha=squeeze(mean(WTData.Absolute.Alpha,2));

                % Plotting the amplitude values (absolute)
                figure
                plot(1:length(WTData.Absolute.Theta),WTData.Absolute.Theta)
                hold on 
                plot(1:length(WTData.Absolute.Alpha),WTData.Absolute.Alpha)
                title('Morlet Wavelet Decomposition')
                xlabel('Epochs')
                ylabel('Amplitude (absolute values)')
                legend('Theta','Alpha')

                % Exporting SleepNoSleep figures
                if sum(~SavePath)<1
                    SaveFigures(gcf,[SavePath '\Exports\' Date_Start '\SleepNoSleep\'...
                        sprintf('WaveletMorletSleep_%d_%s',ParticipantNumber,...
                        Conditions_OrderCell{Pos,WhichCond})],'w','bmp');
                else
                    SaveFigures(gcf,[CurrentPWD '\Exports\' Date_Start '\SleepNoSleep\'...
                        sprintf('WaveletMorletSleep_%d_%s',ParticipantNumber,...
                        Conditions_OrderCell{Pos,WhichCond})],'w','bmp');
                end

                % Statistical decision if asleep or not
                % If = 1 means asleep / = 0 means awake
                for i=1:size(EEG.data,3)
                    if WTData.Absolute.Alpha(i)<WTData.Absolute.Theta(i) 
                        SleepNoSleep(i)=1;
                    else
                        SleepNoSleep(i)=0;
                    end
                end

                % Preparing list of trials in Awake/Asleep states
                TrialsAwake=[];
                TrialsAsleep=[];
                i=1;
                j=1;
                for n=1:length(SleepNoSleep)
                    if SleepNoSleep(n)==0
                        TrialsAwake(i)=n;
                        i=i+1;
                    else
                        TrialsAsleep(j)=n;
                        j=j+1;
                    end
                end

                % EEG will only contain the trials corresponding to AWAKE state
                if ~isempty(TrialsAsleep)
                    EEG = pop_select(EEG, 'trial', TrialsAwake);
                end

                % Indexes to store data in correct columns
                TempIndex = find(contains(lower(SleepNoSleepTable.Properties.VariableNames),...
                    lower(Conditions_Names{WhichCond})));

                % Filling the SleepNoSleep table
                SleepNoSleepTable(Pos,TempIndex(1))={length(TrialsAsleep)};
                SleepNoSleepTable(Pos,TempIndex(2))={length(TrialsAwake)};

                % Unepoching the file 
                % Epoching was only temporary to perform Wavelet Morlet Conv
                EEG = eeg_epoch2continuous(EEG);

                % Removing the events that were created by Epoching
                EventFields={'event','urevent'};
                for t=1:length(EventFields)
                    IdxX = contains({EEG.(EventFields{t}).type},'X'); 
                    IdxBound = contains({EEG.(EventFields{t}).type},'boundary');
                    EEG.(EventFields{t})(or(IdxX,IdxBound))=[];
                end
            end
            
            % Save DATASET 1 - FILTERED (MANDATORY)
            if strcmpi(ExpFilt,'Yes')
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',...
                    sprintf(Dataset_filtCleaned, SubjName,WhichCond)...
                ,'gui','off', 'savenew', [Dir_save '\\' sprintf(Dataset_filtCleaned, SubjName, WhichCond)]);  
            end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  INDEPENDENT COMPONENT ANALYSIS (ICA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ONLY IF ICA SELECTED IN GUI
            if strcmpi(ICA,'Yes')        
                
                % Waitbar updating 
                WaitBarApp.Title = '1. PREPROCESSING: ICA computation';

                % Running ICA decomposition with best algorithm so far (AMICA)
                [W,S,mods] = runamica15(EEG.data,'outdir',...
                    [Dir_load '\\AmicaResults\\' sprintf('Session%d',WhichCond)],...
                    'max_iter',str2double(AnswerICA{:})); 

                % Storing amica results in EEG structure
                EEG.icaweights = W;
                EEG.icasphere = S(1:size(W,1),:);
                EEG.icawinv = mods.A(:,:,1);
                EEG.mods = mods;

                % DIPOLE FITTING OF ICA DERIVED COMPONENTS (SOURCE LOCALIZATION)
                % Tutorials : 
                % https://github.com/kevmtan/EEGpipeline/wiki/Stage-2-%E2%80%93-ICA-&-source-localization#ICA_procedure
                % https://sccn.ucsd.edu/wiki/A08:_DIPFIT#Setting_up_DIPFIT_model_and_preferences
                % Calls the function Automated_Dipfit
                EEG= Automated_Dipfit(EEG);

                % Exporting Dipole Fitting Figures
                ExportString={'LabelDipoles\\LabelsFit%d_%s','DipFit%d_%s',...
                    'LabelDipoles\\LabelsFitBar%d_%s'};
                Color={'w','k','w'};
                for k=length(findobj('type','figure'))-2:length(findobj('type','figure'))
                    SaveFigures(figure(k),[TopoDipFitDirectory...
                        sprintf(ExportString{k},Subjectslist(g),...
                        Conditions_OrderCell{Pos,WhichCond})],Color{k},'bmp');
                end

                % Save DATASET 2 - ICA COMPUTED
                if strcmpi(ExpICAED,'Yes')
                    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',...
                        sprintf(Dataset_filtCleaned_ICAed, SubjName, WhichCond),...
                    'gui','off', 'savenew',[Dir_save '\\' sprintf(Dataset_filtCleaned_ICAed, SubjName, WhichCond)]);  
                end
            end
        end
    end           

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%                            IC REJECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Epitome of UI
    WaitBarApp = uiprogressdlg(App.MainGUIUIFigure,'Title','Progress Bar',...
    'Message','','Cancelable','on');
    File = 0;
    
    % Reloading each file
    for g=1:length(Subjectslist)

        %Current Subject
        ParticipantNumber = Subjectslist(g);

        % Position of the current subject in the subject list
        Pos = find(Subjectslist == ParticipantNumber);

        % Retrieving current subject information 
        % WILL NOT WORK WITH 2 BETWEEN-SUBJ FACTORS !! 
        % IT SHOULD DEPEND ON THE NUMBER OF DIFFERENT FOLDER NAMES !
        FolderTemplate = regexp(SplitFileTemp{end},'\D*','Match');
        
        if length(FilesPath)>1
            CurrentFile = Participant_load.(Conditions_OrderCell{Pos,end}). ...  
            ([FolderTemplate{:} num2str(ParticipantNumber)]);
        else
            CurrentFile = Participant_load.(Groups_Names{1}). ...  
            ([FolderTemplate{:} num2str(ParticipantNumber)]);
        end

        for h=1:length(CurrentFile.FileList)

            % Restart EEGLAB
            eeglab

            % Closing EEGLAB GUI
            close gcf

            %-----------------------------------------------------------------%    
            % Loading directory and templates
            %-----------------------------------------------------------------% 
            % Finding current condition number
            WhichCond = find(contains(lower(Conditions_Names),...
                  lower(CurrentFile.CondAssign{h})));

            % Selecting appropriate groups folder     
            Dir_load = backslash(CurrentFile.Path{:});
            
            % Finding current subject export name
%             CurrentSession = CurrentFile.FileList{h};
%             CurrentSession = strsplit(CurrentSession,'.');
%             SubjName = strtok(CurrentSession{1},FileNames);
%             % If there is no file specific name
%             if sum(isstrprop(SubjName,'digit'))>0
                Temp = strsplit(CurrentFile.Path{:},'\');
                SubjName = [Temp{end} '_'];
%             end
                  
            % If save path different than CurrentPWD
            if ~isempty(SavePath)
                Dir_save = backslash(CurrentFile.ExportPath);
            else
                Dir_save = Dir_load;
            end

            % RE-IMPORT LAST DATASET 
            if strcmpi(ICA,'Yes')
                % IF ICA computed : Dataset_filtered_cleaned_ICAed
                EEG = pop_loadset('filename',[sprintf(Dataset_filtCleaned_ICAed,SubjName,WhichCond) '.set'],...
                    'filepath',Dir_save);
            else
                % IF no ICA: Dataset_filtered_cleaned
                EEG = pop_loadset('filename',[sprintf(Dataset_filtCleaned,SubjName,WhichCond) '.set'],...
                    'filepath',Dir_save);
            end
            
            % Waitbar updating 
            File = File + 1;
            WaitBarApp.Value = File/FilesToProcess;
            WaitBarApp.Title = '1. PREPROCESSING: ICA rejection';
            WaitBarApp.Message = sprintf('Sbj %d/%d : %s',g,...
                length(Subjectslist),CurrentFile.CondAssign{WhichCond});
            
            % Check for Cancel button press
            if WaitBarApp.CancelRequested
                return
            end

            % ONLY IF ICA COMPUTED
            if strcmpi(ICA,'Yes')  
                % Initialize the analysis for each data
                Restart=0;
                
                % ICLabel plugin
                EEG=iclabel(EEG);
                
                % Retrieve results (pre-select all non-brain components)
                Idx=1;
                for l=1:size(EEG.icaact,1)
                    [~,CompType] = max(EEG.etc.ic_classification.ICLabel.classifications(l,:));
                    if CompType ~= 1 % Brain component, see EEG.etc.ic_classification.ICLabel.classes
                        CompsToRej(Idx) = l;
                        Idx = Idx + 1;
                    end
                end
                
                % Variance accounted for by each components
                % Taken from pop_prop_extended.m
                for m=1:size(EEG.icaact,1)
                    maxsamp = 1e5;
                    n_samp = min(maxsamp, EEG.pnts*EEG.trials);
                    try
                        samp_ind = randperm(EEG.pnts*EEG.trials, n_samp);
                    catch
                        samp_ind = randperm(EEG.pnts*EEG.trials);
                        samp_ind = samp_ind(1:n_samp);
                    end
                    if ~isempty(EEG.icachansind)
                        icachansind = EEG.icachansind;
                    else
                        icachansind = 1:EEG.nbchan;
                    end
                    icaacttmp = EEG.icaact(m, :, :);
                    datavar = mean(var(EEG.data(icachansind, samp_ind), [], 2));
                    projvar = mean(var(EEG.data(icachansind, samp_ind) - ...
                    EEG.icawinv(:, m) * icaacttmp(1, samp_ind), [], 2));
                    PVaf(m,1) = 100 *(1 - projvar/ datavar);
                end

                while Restart<1


                    % Visualize the results
                    pop_viewprops( EEG, 0, 1:size(EEG.icaact,1), {'freqrange', [1 60]}, {}, 1, 'ICLabel' )

                    % Move figures all over the screen
                    ScreenPos={'northwest','northeast','southeast','southwest'};
                    for k=1:length(findobj('type','figure'))        
                        movegui(figure(k),ScreenPos{k})
                    end

                    % Matrix to integrate in the following uitable
                    Response = repmat({false},[size(EEG.icaact,1) 1]);
                    Response(CompsToRej)={true};
                    to_display = [num2cell(1:size(EEG.icaact,1))', Response];
                    CompList = {};
                    CompsToRej = [];
                    
                    % Select components to reject
                    % PROBLEM HERE IF YOU DO NOT MODIFY THE SUGGESTION! THE
                    % UITABLE DOESN'T RECORD THE RESPONSES !!! 
                    Screensize = get( groot, 'Screensize' );
                    f = figure('Position', [Screensize(3)/2-200 Screensize(4)/2-200 400 500]);
                    p=uitable('Parent', f,'Data',to_display,'ColumnEdit',[false true],'ColumnName',...
                        {'Components', 'REJECTION?'},'CellEditCallBack','CompList = get(gco,''Data'');');
                    uicontrol('Style', 'text', 'Position', [0 400 400 80], 'String',...
                            {'SELECTION OF COMPONENTS TO REJECT',...
                            'Click on the box corresponding to the component(s) you want to reject.',...
                            'Components already selected correspond to all non-brain components detected by the algorithm.'});
                    % Adding the Clear button
                    p=ClearButton(p);
                    
                    % Wait for t to close until running the rest of the script
                    waitfor(p)
                    
                    % Saving the changes
                    if ~isempty(CompList)
                        CompsToRej = find(cell2mat(CompList(:,2))~=0)';
                        RemainPfav = sum(PVaf)-sum(PVaf(CompsToRej));
                    else
                        CompsToRej = find(cell2mat(Response)==1)';
                        RemainPfav = sum(PVaf)-sum(PVaf(CompsToRej));
                    end

                    % close all figures
                    close all

                    % Reject the marked components and save the new dataset
                    if ~isempty(CompsToRej)

                        % Removing components
                        TEMPEEG = pop_subcomp( EEG, CompsToRej, 0);

                        % Visual check before/after interpolation
                        vis_artifacts(TEMPEEG,EEG);    

                        % Wait Bar 
                        Fig=msgbox(['Take time to visualize the difference before and after IC rejection!'... 
                            newline 'THE CODE WILL CONTINUE ONCE YOU PRESS OK' ...
                            newline newline sprintf('! Currently you are keeping %s%% of the data !',...
                            num2str(round(RemainPfav,2)))],'WAIT','warn'); 
                        uiwait(Fig);
                        close all
                    end

                    % Check to be sure of rejection
                    PromptSureRejICA = inputdlg(['Are you sure you want to reject to following components? '...
                        newline 'If yes, press OK'....
                        newline 'If not, press CANCEL and the code will restart!'],...
                        'Components to reject',5,{num2str(CompsToRej)});

                    % Restarts the loop if user pressed "Cancel"
                    if ~isempty(PromptSureRejICA)
                        % Ends the loop 
                        Restart = 1;

                        % Reject the selected IC
                        EEG=pop_subcomp( EEG, CompsToRej, 0);
                    end

                    % Close all figures
                    close all
                end

                % Save DATASET 3 - PRUNED WITH ICA
                if strcmpi(ExpICARej,'Yes')
                    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',...
                        sprintf(Dataset_filtCleaned_ICAedRejected, SubjName, WhichCond),...
                    'gui','off', 'savenew',[Dir_save '\\' sprintf(Dataset_filtCleaned_ICAedRejected,SubjName,WhichCond)]);
                end
            end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%             BAD CHANNELS INTERPOLATION AND AVERAGE REFERENCING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             ErrorArtifacts = {};
%             r = 1;
%             try
                % Reintroduce the bad channels data 
                Temp = zeros(length(Channels),size(EEG.data,2)); PosGood = 1; PosBad = 1;
                for m=Channels
                   if ~ismember(m,EEG.BadChans.InterpChans) 
                       Temp(m,:) = EEG.data(PosGood,:); PosGood = PosGood + 1;
                   else
                       % Restricting channel data length since EEG.data
                       % size might have changed with artifact rejection
                       Temp(m,:) = EEG.BadChans.data(PosBad,1:size(EEG.data,2));PosBad = PosBad + 1;
                   end
                end
                
                % Adjust the EEG structure
                EEG.data = Temp; EEG.chanlocs = EEG.BadChans.chanlocs;
                EEG.nbchan = EEG.BadChans.nbchan;
                EEG = eeg_checkset(EEG);
                
                % Multiquadratics bad channels interpolation
                EEG.data = EEGinterp('MQ',0.05,EEG,EEG.BadChans.InterpChans);
                
                % Average referencing
                EEG = average_ref(EEG);
                
                % Use this function from the PrepPipeline to perform "Robust
                % Referencing" and Detection of Bad Channels Rejection/Interpolation
%                 [EEG, InterpChanStruct] =  performReference(EEG);
%                 EEG.etc.noiseDetection.reference = InterpChanStruct;
%                 EEG.etc.noiseDetection.fullReferenceInfo = true;
%                 EEG.etc.noiseDetection.interpolatedChannelNumbers = ...
%                     InterpChanStruct.interpolatedChannels.all;
%                 EEG.etc.noiseDetection.stillNoisyChannelNumbers = ...
%                     InterpChanStruct.noisyStatistics.noisyChannels.all;
%             catch
                % Write error to the LOG
%                 ErrorArtifacts{r} = [SubjName FileNames num2str(WhichCond)];
%                 r = r+1;
%                 break
%             end
            
            % Visual check before/after interpolation
            if strcmpi(Automaticity,'Yes')==1
                vis_artifacts(EEG,OriginalEEG);
            % Holds the figure until inspection is over
                Fig=msgbox('THE CODE WILL CONTINUE ONCE YOU PRESS OK','WAIT','warn'); 
                uiwait(Fig);
                close gcf
            end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%             EXPORTING PREPROCESSED RESULTS In .BDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmpi(ExpBDF,'Yes')
                
                % Waitbar updating 
                WaitBarApp.Title = '1. PREPROCESSING: BDF export';

                if End==0
                    EXPORTEEG = EEG;
                else
                    % Epoching files (temporary, to avoid writeeg to fill in empty
                    % blocks with NaNs)
                    EXPORTEEG = eeg_regepochs(EEG, 1); % 1sec epochs (arbitrary)
                    % Removing the events that were created by ASR (may be useless)
                    Eventfields = {'event','urevent'};
                    for p=1:2
                        if ~isempty(EXPORTEEG.(Eventfields{p}))
                            if nnz(contains({EXPORTEEG.(Eventfields{p}).type},'X'))>0
                                IdxX = contains({EXPORTEEG.(Eventfields{p}).type},'X'); 
                                EXPORTEEG.(Eventfields{p})(IdxX)=[];
                            end
                        end
                    end
                end

                % Export as .BDF
                % Running the script below should fix unkown error 
                x = fileparts( which('sopen') );
                rmpath(x);
                addpath(x,'-begin');
                if ~isempty(EXPORTEEG.event)
                    writeeeg([Dir_save '\\' sprintf(PreprocessedEEG,SubjName,WhichCond)],EXPORTEEG.data,...
                        EXPORTEEG.srate,'TYPE','BDF','Label',{EXPORTEEG.chanlocs.labels},...
                        'Patient.id',strrep(SubjName,'_',''),'EVENT',EXPORTEEG.event);
                else
                    writeeeg([Dir_save '\\' sprintf(PreprocessedEEG,SubjName,WhichCond)],EXPORTEEG.data,...
                    	EXPORTEEG.srate,'TYPE','BDF','Label',{EXPORTEEG.chanlocs.labels},...
                        'Patient.id',strrep(SubjName,'_',''));
                end
            end
            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%             EXPORTING PREPROCESSED DATASET 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%             % Finding current subject export name
%             CurrentSession = CurrentFile.FileList{h};
%             CurrentSession = strsplit(CurrentSession,'.');
%             SubjName = strtok(CurrentSession{1},FileNames);
%             % If there is no file specific name
%             if sum(isstrprop(SubjName,'digit'))>0
                Temp = strsplit(CurrentFile.Path{:},'\');
                SubjName = [Temp{end} '_'];
%             end

            % Save DATASET 4 - ClEANED/CHANNELED
            if strcmpi(ExpClean,'Yes')
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',...
                    sprintf(Dataset_interp, SubjName, WhichCond),...
                'gui','off', 'savenew', [Dir_save '\\' sprintf(Dataset_interp, SubjName, WhichCond)]); 
            end
            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         EXCEL FILES EXPORT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
            % Waitbar updating 
            WaitBarApp.Title = '1. PREPROCESSING: Excel templates';
    %-------------------------------------------------------------------------%
            % INTERPOLATED CHANNELS 
    %-------------------------------------------------------------------------%        
            % Stores interp channels and in text file in each participant's
            % folder

            % List of interpolated channels
            if ~isempty(EEG.BadChans.InterpChans)
                if size(EEG.BadChans.InterpChans,1) == 1
                    InterpChans = EEG.BadChans.InterpChans;
                else
                    InterpChans = EEG.BadChans.InterpChans';
                end
            else
                InterpChans = [];
            end

            % This function is used to determine the excel column names 
            % corresponding to the number of interpolated channels
            XlsCompCol=xlsColNum2Str(length(InterpChans)+1); 

            % Adding names of electrodes
            InterpChansStr = {EEG.chanlocs(InterpChans).labels};
            InterpChans = cellfun(@(x) ['(' num2str(x) ')'], num2cell(InterpChans),'UniformOutput',false);
            InterpChans = strcat(InterpChansStr,InterpChans);

            % Exporting the list of rejected components in specific excel files
            if ~isempty(InterpChans)
                xlswrite([ExcelDirectory 'InterpChannels' Conditions_Names{WhichCond} '.xlsx'],...
                    InterpChans,sprintf('B%d:%s%d',Pos+1,...
                    XlsCompCol{1},Pos+1));
            end

%-------------------------------------------------------------------------%
                      % ICA REJECTED COMPONENTS
%-------------------------------------------------------------------------%         
            if strcmpi(ICA,'Yes')         
                % This function is used to determine the excel column names 
                % corresponding to the number of rejected ICA components
                XlsCompCol=xlsColNum2Str(length(CompsToRej)+1); 
                
                % List of rejected components type
                ICALabels = cell(1,length(CompsToRej));
                for t = 1:length(CompsToRej)
                    [~,TypePos] = max(EEG.etc.ic_classification.ICLabel.classifications(CompsToRej(t),:));
                    ICALabels(t) = strcat(EEG.etc.ic_classification.ICLabel.classes(TypePos),sprintf('(%s)',num2str(CompsToRej(t))));
                end

                % Exporting the list of rejected components in specific excel files
                if ~isempty(CompsToRej)
                    xlswrite([ExcelDirectory 'RejectedComponents' Conditions_Names{WhichCond} '.xlsx'],... % Conditions_OrderCell{Pos,WhichCond}
                        ICALabels,sprintf('B%d:%s%d',Pos+1,...
                        XlsCompCol{1},Pos+1));
                end
            end
%-------------------------------------------------------------------------%
                      % SLEEP/NO SLEEP EXPORTS
%-------------------------------------------------------------------------%
            % Exporting to Excel file
            if strcmpi(Sleep, 'Yes') 
                writetable(SleepNoSleepTable,backslash(strcat(ExcelDirectory,'AsleepAwakeTrials.xlsx')));
            end
        end
    end
else
    ErrorArtifacts = [];
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       FREQUENCY BANDS ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Epitome of UI
WaitBarApp = uiprogressdlg(App.MainGUIUIFigure,'Title','Progress Bar',...
    'Message','','Cancelable','on');
File = 0;

if Analysis && strcmpi(Steps,'Preprocessing') || strcmpi(Steps,'Both')
        % Finally, here we perform the power spectra analysis on specific
        % frequency bands. At the same time, topoplots are exported in
        % designed folder. 
    for g=1:length(Subjectslist)

        % Current Subject
        ParticipantNumber = Subjectslist(g);

        % Position of the current subject in the subject list
        Pos = find(Subjectslist == ParticipantNumber);

        % Retrieving current subject information 
        % WILL NOT WORK WITH 2 BETWEEN-SUBJ FACTORS !! 
        FolderTemplate = regexp(SplitFileTemp{end},'\D*','Match');
        if length(FilesPath)>1
            CurrentFile = Participant_load.(Conditions_OrderCell{Pos,end}). ...  
            ([FolderTemplate{:} num2str(ParticipantNumber)]);
        else
            CurrentFile = Participant_load.(Groups_Names{1}). ...  
            ([FolderTemplate{:} num2str(ParticipantNumber)]);
        end

        for h=1:length(CurrentFile.FileList)    

            % Restart EEGLAB
            eeglab
            % Close it automatically (too many figures)
            close gcf

            %-------------------------------------------------------------%    
            % Loading directory and templates
            %-------------------------------------------------------------% 
            % Finding current condition number
            WhichCond = find(contains(lower(Conditions_Names),...
                  lower(CurrentFile.CondAssign{h})));

            % Selecting appropriate groups folder     
            Dir_load = backslash(CurrentFile.Path{:});
            
            % Finding current subject export name
%             CurrentSession = CurrentFile.FileList{h};
%             CurrentSession = strsplit(CurrentSession,'.');
%             SubjName = strtok(CurrentSession{1},FileNames);
%             % If there is no file specific name
%             if sum(isstrprop(SubjName,'digit'))>0
                Temp = strsplit(CurrentFile.Path{:},'\');
                SubjName = [Temp{end} '_'];
%             end
            
            % If save path different than CurrentPWD
            if ~isempty(SavePath)
                Dir_save = backslash(CurrentFile.ExportPath);
            else
                Dir_save = Dir_load;
            end

            % RE-IMPORT LAST DATASET 
            if strcmpi(ICA,'Yes')
                % IF ICA computed : Dataset_filtered_cleaned_ICAedRejected
                EEG = pop_loadset('filename',...
                    [sprintf(Dataset_filtCleaned_ICAedRejected,SubjName,WhichCond) '.set'],...
                    'filepath',Dir_save);
            else
                % IF no ICA: Dataset_filtered_cleaned
                EEG = pop_loadset('filename',...
                    [sprintf(Dataset_interp,SubjName,WhichCond) '.set'],...
                    'filepath',Dir_save);
            end
            
            % Waitbar updating 
            File = File + 1;
            WaitBarApp.Value = File/FilesToProcess;
            WaitBarApp.Title = '2. ANALYSES: Power spectral density (PSD)';
            WaitBarApp.Message = sprintf('Sbj %d/%d : %s',g,...
                length(Subjectslist),CurrentFile.CondAssign{WhichCond});
            
            % Check for Cancel button press
            if WaitBarApp.CancelRequested
                return
            end
            
            %-------------------------------------------------------------%    
            % POWER SPECTRA ANALYSIS
            %-------------------------------------------------------------% 
            
            % Detecting the frequency resolution based on user input
            if nnz(mod(FreqRanges,1) == 0) < numel(FreqRanges)
                FreqToRes = find(mod(FreqRanges,1)~=0); FreqRes = 0;
                for f=1:length(FreqToRes) % For each frequency with resolution lower than integer
                    if FreqRes < 1/mod(FreqRanges(FreqToRes(1)),1)
                        FreqRes = 1/mod(FreqRanges(FreqToRes(1)),1);
                    end
                end
            else
                FreqRes = 1; 
            end
            
            % Computing PSD for all channels (= GPS)
            figure;[SpectOutputs,FreqsSpec]=pop_spectopo(EEG, 1, [], 'EEG','percent', ...
            100, 'freq', FreqsToPlot,'freqfac',FreqRes, 'freqrange',...
            [min(FreqRanges(1,1)) max(FreqRanges(end,2))],'electrodes','on'); 
        
            % Removing the first column of data (0 Freq = DC offset)
            SpectOutputs = SpectOutputs(:,2:end);
            FreqsSpec = FreqsSpec(2:end);

            % Saving topoplots
            SaveFigures(gcf,[TopoplotsDirectory sprintf('PowerSpectrum%d_%s',ParticipantNumber,...
                Conditions_OrderCell{Pos,WhichCond})],'w','bmp');
            
            % Saving matrix with PSD for each subject
            save([Dir_save '\\' sprintf(PSDEEG,SubjName,WhichCond)],'SpectOutputs');

            %-------------------------------------------------------------%
            % FREQUENCIES EXPORT
            %-------------------------------------------------------------%
            % If range has been defined for a specific frequency at the beginning of
            % the script, the frequency power will be exported. Otherwise, it will be
            % omitted. Files are exported in .EP format, which is compatible with
            % STEN and CARTOOL softwares.
                    
            for m=1:length(FreqRanges)

                % LOOP
                LowIdx(m) = find(FreqsSpec==FreqRanges(m,1));
                HighIdx(m) = find(FreqsSpec==FreqRanges(m,2));
                FreqRangesTemp=[LowIdx(m) HighIdx(m)]; % Index of frequency bounds
                FreqTemp=mean(SpectOutputs(:,FreqRangesTemp(1):FreqRangesTemp(2)),2); 
                FreqTempEP=FreqTemp';

                % Save to specific STATS folder 
                DirectoryFreq=[StatsDirectory FreqNames{m} '\',...
                sprintf(Frequency_export{m}, Subjectslist(g)) Conditions_OrderCell{Pos,WhichCond} '.ep']; 
                save(DirectoryFreq, 'FreqTempEP', '-ascii');
                
                for r=1:length(Conditions_Names)
                    % Fill in the GPS Matrix
                    GPS.(FreqNames{m}).(Conditions_Names{WhichCond})(Pos,2:end)=array2table(FreqTempEP);
                    % Rename the channels according to defined standards
                    GPS.(FreqNames{m}).(Conditions_Names{WhichCond}).Properties.VariableNames=[{'Participants'} {EEG.chanlocs.labels}];
                    % Export to GPS excel file
                    writetable(GPS.(FreqNames{m}).(Conditions_Names{WhichCond}),...
                        [ExcelDirectory 'GPS_' FreqNames{m} '.xlsx'],...
                        'Sheet',[FreqNames{m} '_' lower(Conditions_Names{WhichCond})]);
                end
           
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           POWER SPECTRAL AMPLITUDE LOCALIZATION ANALYSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This complementary analysis measures for each frequency band
            % what are the highest and lowest spectral amplitude (peak) 
            % and where they are recorded on the EEG scalp at the level of:
            % Areas (user defined) 
            % Then the amplitude of the contralateral electrode/area is
            % computed in order to compare the results with the ipsilateral
            % maximum amplitude side.

                if strcmpi(AnalysesSwitch{2,2},'Yes')

                    % Creating a structure from the initial AreasList matrix
                    FieldsAreas=fieldnames(NewAreasList);

                    % Waitbar updating 
                    WaitBarApp.Title = '2. ANALYSES: PSD Electrodes clusters';

            %---------------------------% AREA LEVEL %------------------------%
                    % Spectra max amplitude measurement
%                     for m=1:size(FreqNames,1)
                    SpectOutputsStruc.Data.(FreqNames{m})=mean(SpectOutputs(:,[LowIdx(m) HighIdx(m)]),2); 
                        % Maximum amplitude for all freq ranges peak
    %                     SpectOutputsStruc.MaxAmp.(FreqNames{m})=...
    %                         max(SpectOutputsStruc.Data(:).(FreqNames{m})); 
%                     end

                    % Creating a structure from the initial AreasList matrix
                   %  FieldsAreas=fieldnames(NewAreasList);
    %                 AllElectrodes = [];
    %                 for t=1:length(FieldsAreas)
    %                     AllElectrodes = [AllElectrodes NewAreasList.(FieldsAreas{t})];
    %                 end

                    % Finding the electrodes-cluster with max/min amplitude
%                     for m=1:size(FreqNames,1)
                    for f=1:length(FieldsAreas)
                        % [MeanClust(m,f)] = mean(SpectOutputsStruc.Data.(FreqNames{m})(AllElectrodes(:,f)));
                        [MeanClust(m,f)] = mean(SpectOutputsStruc.Data.(FreqNames{m})(NewAreasList.(FieldsAreas{f})));
                    end
%                     end
            
                
        %----------------------% AMPLITUDE TABLES %-----------------------%            
                % Amplitude tables for frequency ranges consisting
                % of min/max amplitude values and corresponding electrodes cluster index

                    % Find positions of columns to store values according to
                    % conditions
                    Max = ismember(lower(AreaAmpTable.(FreqNames{m}).Properties.VariableNames),...
                        lower([Conditions_Names{WhichCond} '_Max']));
                    Min =ismember(lower(AreaAmpTable.(FreqNames{m}).Properties.VariableNames),...
                        lower([Conditions_Names{WhichCond} '_Min']));
                    MaxAmpClust = ismember(lower(AreaAmpTable.(FreqNames{m}).Properties.VariableNames),...
                        lower([Conditions_Names{WhichCond} '_MaxClust']));
                    MinAmpClust = ismember(lower(AreaAmpTable.(FreqNames{m}).Properties.VariableNames),...
                        lower([Conditions_Names{WhichCond} '_MinClust']));

                % Store values in corresponding tables
                % for m=1:size(FreqRanges,1)
                    
                    % Finding corresponding electrode-cluster area
                    [MaxClust,FindMax] = max(MeanClust(m,:));
                    [MinClust,FindMin] = min(MeanClust(m,:));
                    TempMaxClustName = FieldsAreas(FindMax);
                    TempMinClustName = FieldsAreas(FindMin);
       
                    % Temp variables 
                    TempCell=table2cell(AreaAmpTable.(FreqNames{m}));
                    TempHeader=AreaAmpTable.(FreqNames{m}).Properties.VariableNames;

                    % Areas level
                    TempCell(Pos,Max)={MaxClust}; % Max
                    TempCell(Pos,Min)={MinClust}; % Min
                    TempCell(Pos,MaxAmpClust)=TempMaxClustName; % Max amplitude area
                    TempCell(Pos,MinAmpClust)=TempMinClustName; % Min amplitude area

                    % Storing the new values in the AreaAmpTable structure
                    AreaAmpTable.(FreqNames{m})=cell2table(TempCell);

                    % Restores headers
                    AreaAmpTable.(FreqNames{m}).Properties.VariableNames=TempHeader;

                    % Exporting the tables in predefined excel files
                    writetable(AreaAmpTable.(FreqNames{m}),...
                        [ExcelDirectory 'AreaAmplitude' FreqNames{m} '.xlsx']);
                % end  
                end
            end
        end
    end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  SAVING VARIABLES FOR FURTHER USE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([SavePath '\Parameters\' Date_Start '\MAINWorkspace.mat'],...
    'Dataset_filtCleaned_ICAedRejected','Dataset_filtCleaned',...
    'Participant_load','Groups_Names','Conditions_Names','Subjectslist',...
    'FreqRanges','ErrorArtifacts','ExcelDirectory');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  RUNNING THE STUDY (GROUP ANALYSES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(Steps,'Both') || strcmpi(Steps,'Group analyses') 
    % Epitome of UI
    WaitBarApp = uiprogressdlg(App.MainGUIUIFigure,'Title','Progress Bar',...
        'Message','','Cancelable','on');
    if sum(~SavePath)<1
        GROUPSTUDY(SavePath,Date_Start,WaitBarApp)
    else
        GROUPSTUDY(CurrentPWD,Date_Start,WaitBarApp)
    end
else
    % Running the log
    if sum(~SavePath)<1
        LOG(SavePath,Date_Start,'ErrorsPreproc',ErrorArtifacts)
    else
        LOG(CurrentPWD,Date_Start,'ErrorsPreproc',ErrorArtifacts)
    end
end
end