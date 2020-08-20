%%------------------------------RESTINGLAB-------------------------------%%

% Version 0.60
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
                                
                                
% TO IMPLEMENT:
% DO THE SOURCE LOCALISATION AND INCLUDE IN PLOT ??? 
% Could do that w/ ICA using MST toolbox ? 
% Or implement Cartool method! 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%                     STUDY CREATION FOR GROUP ANALYSIS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function GROUPSTUDY(App,SavePath,Date_Start)

% Epitome of UI
WaitBarApp = uiprogressdlg(App.MainGUIUIFigure,'Title','Progress Bar',...
    'Message','','Cancelable','on');

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
        load([TempList(SmallDistDateIdx).folder '\' TempList(SmallDistDateIdx).name])
    else
        load([TempList.folder '\' TempList.name])
    end
end

% Adding path to dependencies
addpath([pwd '\Functions\eeglab2020_0']);
addpath([pwd '\Functions\']);
addpath(genpath([pwd '\Functions\Dependencies']));
addpath([pwd '\Functions\EEGInterp']); 

% Need to intialize it otherwise it doesn't recognise std_editset
STUDY = []; CURRENTSTUDY = 0; ALLEEG=[]; EEG=[]; CURRENTSET=[]; 
eeglab
close gcf

% This allows to process more datasets while only keeping 1 in memory! 
pop_editoptions( 'option_storedisk',1);

% Looking for ICA computed file
ICAexist = 'No';
if sum(~SavePath)<1
    TempList =  dir([SavePath '\**/*' '.set']); 
    for l=1:length({TempList.name})

        % Looking for ICA pattern in files
        if nnz(contains({TempList.name},'ICA')) > 0
            ICAexist = 'Yes';
        end
    end
else
    for k=1:length(FilesPath)
        TempList =  dir([FilesPath{k} '\**/*' '.set']); 
        for l=1:length({TempList.name})

            % NEW METHOD
            if nnz(contains({TempList.name},'ICA')) > 0
                ICAexist = 'Yes';
            end
        end
    end
end

% Update progress
WaitBarApp.Title = '3. STUDY ANALYSES: Building structure';

% Check for Cancel button press
if WaitBarApp.CancelRequested
    return
end

% Dataset to load
if strcmpi(ICAexist,'Yes')
    DatasetToLoad = Dataset_filtered_cleaned_ICAedRejected;
else
    DatasetToLoad = Dataset_filtered_cleaned;
end

% Only keeping subjects selected for analysis
Conditions_Order=readtable([DirectoryCond FileCond]);
TempSubjectslist = table2cell(Conditions_Order(:,1));
if isinteger(TempSubjectslist{1})
TempSubjectslist=cell2mat(TempSubjectslist); 
else; TempSubjectslist=cellfun(@(x) str2double(x),TempSubjectslist); 
end
SubjectsRemoved = []; Subjectslist = [];
j=1;
t=1;

for k=1:length(TempSubjectslist)
    for l=1:size(ProcessDataSTUDY,2) % Number of sessions  
        if ProcessDataSTUDY{k,l}
            Subjectslist(j) = TempSubjectslist(k);
            j=j+1;
            break
        end
        
        % Scan the ProcessData line to remplace empty fields by 0
        ProcessDataSTUDY(k,cellfun(@(x) isempty(x),ProcessDataSTUDY(k,:))) = {false};
        
        % This will detect files for which there is no data to analyze
        if sum(cell2mat(ProcessDataSTUDY(k,:)))<1
            SubjectsRemoved(t) = TempSubjectslist(k);
            t=t+1;
            break
        end
    end
end

% Building the STUDY structure
Increment=1;
Sbj = 0;
% Number of between-subject (BS) factor
FieldsBS = fieldnames(Participant_load);

% Update waitbar
WaitBarApp.Title = '3. STUDY ANALYSES: Loading data';

% Loop for BS factor
for i=1:numel(FieldsBS)

    % Number of subjects in each BS levels
    FieldsBSSubj=fieldnames(Participant_load.(Groups_Names{i}));
    
    % Remove subjects not included in analysis (NEW)
    TempFieldsBSSubj=cellfun(@(x) regexp(x,'\d*','Match'),FieldsBSSubj);
    TempFieldsBSSubj=cellfun(@(x) str2num(x),TempFieldsBSSubj);
    FieldsBSSubj(ismember(TempFieldsBSSubj,SubjectsRemoved))=[];
    
    % Sort in natural order
    FieldsBSSubj = sort_nat(FieldsBSSubj);

    for l=1:numel(FieldsBSSubj)
        
        % Subjects counter
        Sbj = Sbj + 1;
        
        % Number of within-subject (WS) levels
        Cond=Participant_load.(Groups_Names{i}).(FieldsBSSubj{l}).CondAssign;

        %Current Subject
        ParticipantNumber = regexp(FieldsBSSubj{l},'\d*','Match');

        for m=1:length(Cond)

            % Finding current condition number
            WhichCond = find(contains(lower(Conditions_Names),...
                  lower(Cond{m})));

            % NEED TO IMPROVE IT IF 2 WS FACTORS ! (1 should be session and
            % the second condition) ?
            CurrentFolder = Participant_load.(Groups_Names{i}).(FieldsBSSubj{l});

            % Finding current subject export name
            CurrentSession = CurrentFolder.FileList{m};
            CurrentSession = strsplit(CurrentSession,'.');
            SubjName = strtok(CurrentSession{1},FileNames);
            % If there is no file specific name
            if sum(isstrprop(SubjName,'digit'))>0
                Temp = strsplit(CurrentFolder.Path{:},'\');
                SubjName = [Temp{end} '_'];
            end
            
             % Update waitbar
             WaitBarApp.Value = Increment/length(TempList);
             WaitBarApp.Message = sprintf('Sbj %d/%d : %s',Sbj,...
                length(Subjectslist),CurrentFolder.CondAssign{WhichCond});

            % Loading each dataset
            [STUDY ALLEEG] = std_editset( STUDY, ALLEEG,'name', StudyName, 'commands',...
            {{'index' Increment 'load' [CurrentFolder.ExportPath '\' sprintf(DatasetToLoad, SubjName, WhichCond) '.set'] ...
            'subject' num2str(ParticipantNumber{:}) 'session' WhichCond 'condition' Conditions_Names{WhichCond} 'group' FieldsBS{i}}},...  
            'updatedat','off','rmclust','on' );           
             Increment=Increment+1;
        end
    end
end    

% Check for downsampling mismatch between the datasets (may cause errors)
UniqueSampRate = unique(cell2mat({ALLEEG.srate}));
ErrorLog = {};
if length(UniqueSampRate) > 1
    ErrorLog = {num2str(min(UniqueSampRate))};
    for k=1:length(ALLEEG)
        if ALLEEG(k).srate > min(UniqueSampRate)
            ALLEEG(k) = pop_resample(ALLEEG(k), SamplingRate);
            ALLEEG(k) = pop_saveset(ALLEEG(k),'filename',ALLEEG(k).filename,...
                'filepath',ALLEEG(k).filepath);
            
            % Write to LOG
            ErrorLog = [ErrorLog;ALLEEG(k).setname];
        end
    end
end

% Load a template EEG
TemplateEEG = pop_loadset('filename',ALLEEG(1).filename,'filepath',ALLEEG(1).filepath);

%% DEFINING DESIGN

% Update progress
WaitBarApp.Title = '3. STUDY ANALYSES';
WaitBarApp.Value = 1/3;
WaitBarApp.Message = 'Create design';

% Retrieve number of factors
WithinFactors = WithinFactors(~cellfun('isempty',WithinFactors));
BetweenFactors = BetweenFactors(~cellfun('isempty',BetweenFactors));
SubjectsDesign = cellfun(@(x) num2str(x),num2cell(Subjectslist),'UniformOutput', 0);

% BUILD DESIGN
if length(STUDY.group) == 2
    Design.Between = [{'B'};BetweenFactors;STUDY.group'];
end
if length(STUDY.condition) == 2
    Design.Within = [{'W'};WithinFactors;STUDY.condition'];
end

% First variable is always between and second within

% Dependent-samples t-tests & RM ANOVA
if isempty(BetweenFactors) && length(WithinFactors) == 1
    % Design definition
    STUDY = std_makedesign(STUDY, ALLEEG, 1, 'variable1','','variable2','condition',...
        'name',StudyName,'pairing1','off','pairing2','on','delfiles','off',... % delfiles = limites crashes !!! 
        'defaultdesign','off','values1',Conditions_Names,'subjselect',SubjectsDesign);

% Independent-samples t-tests & One-Way ANOVA
elseif isempty(WithinFactors) && length(BetweenFactors) == 1
    % Design definition 
    STUDY = std_makedesign(STUDY, ALLEEG, 1, 'variable1','group','variable2','',...
        'name',StudyName,'pairing1','off','pairing2','on','delfiles','off',...
        'defaultdesign','off','values1',Groups_Names,'subjselect',SubjectsDesign);

% Mixed ANOVA
elseif length(BetweenFactors) == 1 && length(WithinFactors) == 1
    STUDY = std_makedesign(STUDY, ALLEEG, 1, 'variable1','group','variable2',...
        'condition','name',StudyName,'pairing1','off','pairing2','on','delfiles','off',...
        'defaultdesign','off','values1',Groups_Names,'values2',Conditions_Names,...
        'subjselect',SubjectsDesign);

% Repeated-measures ANOVA
% elseif length(WithinFactors) == 2 && isempty(BetweenFactors)
%     STUDY = std_makedesign(STUDY, ALLEEG, 1, 'variable1','condition','variable2',...
%         'session','name',StudyName,'pairing1','on','pairing2','on','delfiles','limited',...
%         'defaultdesign','off','values1',Conditions_Names,'values2',Conditions_Names,...
%         'subjselect',SubjectsDesign); % Error with conditions_names !!! 

% ??? ANOVA
% elseif length(BetweenFactors) == 2 && isempty(WithinFactors)
end

% Retrieve Channels Labels
ChannelsLabels={ALLEEG(1).chanlocs.labels};

% % Selecting which design to compute analysis on 
% STUDY = std_selectdesign(STUDY, ALLEEG, x);

% Retrieve frequency bins
FreqRanges = cell2mat(cellfun(@(x) str2num(x), FreqData(:,2),'UniformOutput',0));

% Define the boundaries based on LowPass if defined or freq bins
if nnz(LowPass)>0
    ImportFreqRange = [0 LowPass];
elseif nnz(FreqRanges)>0
    ImportFreqRange = [0 max(FreqRanges(:,2))];
else
    ImportFreqRange = [];
end

% Update progress
WaitBarApp.Value = 2/3;
WaitBarApp.Message = 'Precompute power spectra';

% Precompute Channel Power Spectra
STUDY.SpecMode = 'psd'; % THIS SHOULD BE IN GUI !
% Don't know what to do since for individual subjects we force use of PSD!!
[STUDY ALLEEG] = std_precomp(STUDY, ALLEEG,'channels','spec','on','recompute','on',...
    'specparams',{'specmode',STUDY.SpecMode,'logtrials','off'});

% Save STUDY
[STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename',[StudyName '.study'],...
'filepath',[SavePath '\STUDY\' Date_Start '\']);

%% PLOT POWER SPECTRA 

% Update progress
WaitBarApp.Value = 3/3;
WaitBarApp.Message = 'Plot subjects power spectra';

% Reading spectral data (precomputed)
try % Sometimes files get corrupted for unknown reasons
    [STUDY,SpectData,SpectFreqs] =  std_readdata (STUDY, ALLEEG,'channels',...
        ChannelsLabels,'datatype','spec','freqrange',ImportFreqRange);
    
catch
    % Precompute Channel Power Spectra
    [STUDY ALLEEG] = std_precomp(STUDY, ALLEEG,'channels','spec','on','recompute','on',...
        'specparams',{'specmode',STUDY.SpecMode,'logtrials','off'});
    [STUDY,SpectData,SpectFreqs] =  std_readdata (STUDY, ALLEEG,'channels',...
        ChannelsLabels,'datatype','spec','freqrange',ImportFreqRange);
end

% Saving spectral data and statistics in a .mat file
save([SavePath '\STUDY\' Date_Start '\SpectralData.mat'],'SpectData');

% SpecData structure is based on STUDY.design.cell order !!!!! 
% Thus, changing the order of STUDY.group! 
STUDY.group = unique({STUDY.datasetinfo.group},'stable');

% Plotting subjects boxplots to identify potential outliers
Pos = 1;
for p=1:length(STUDY.group)
    Dat = []; Color = [];GrpLab = {};
    
    % For each condition
    for m=1:length(STUDY.condition)
        
        % Build data for plotting
        Dat = [Dat;reshape(squeeze(mean(SpectData{p,m},2)),[size(SpectData{p,m},1)*size(SpectData{p,m},3),1])];
        Color = [Color;repmat(STUDY.condition(m),size(Dat,1)/m,1)];
        
        % Groups labels
        for n=1:size(SpectData{p,m},3)
           GrpLab = [GrpLab;repmat({sprintf('P%s',STUDY.datasetinfo(Pos).subject)},...
               size(SpectData{p,m},1),1)];
           Pos = Pos + 1;
        end
       
    end
        
    % Plotting with GRAMM
    Graph=gramm('x',GrpLab,'y',Dat,'color',Color);
    Graph.stat_boxplot(); % 'width',0.1,'dodge',0.1
    Graph.set_names('x','','y','10*Log10(\muV^2/Hz)','color','Conditions'); 
    Graph.set_title(STUDY.group{p});
    figure('units','normalized','outerposition',[0 0 1 1]); Graph.draw();
    clear Graph

    % Channels-averaged data
    SaveFigures(gcf,[SavePath '\STUDY\' Date_Start '\' sprintf('ChanAVG_%s_%s',...
        STUDY.group{p},STUDY.name)],'w','bmp');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STEP 1 : STATISTICS on ALL DATA (channels X frequencies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running Statistics   

% THE GENERALIZATION TO WITHIN/BETWEEN DESIGNS (E.G. MIXED AOV) WILL DEPEND
% ON HOW I STRUCTURE THE DATA! THEN WILL NEED TO ADAPT THE SUBPLOTS BASED
% ON COLUMN/LINES STRUCTURE OF DATA.
% --> Should use same structure as eeglab stats function require! 

if strcmpi(FreqBandsAnalyses,'yes')
    
    % Update progress
    WaitBarApp.Value = 1/(2+size(FreqData,1));
    WaitBarApp.Message = 'Statistics on channels X frequencies data';

    % T-TESTs
    if length(STUDY.group)==2 && length(STUDY.condition)<2 || ...
        length(STUDY.group)<2 && length(STUDY.condition)==2
        Test = 't-test';

        % Building AllData structure
         % THIS WON'T WORK WITH dependent-samples t-tests ??? 
        FullSpectData = [];
        GroupTemp = {};
        for t = 1:length(STUDY.group) 
            GroupTemp = [GroupTemp;repmat(STUDY.group(t),size(SpectData{t},3),1)];

            FullSpectData = cat(3,FullSpectData,SpectData{t});
        end
        AllDataFields = strcat(GroupTemp,'_',{STUDY.datasetinfo.subject}');
        for t = 1:length(AllDataFields)
            AllData.(AllDataFields{t}) = FullSpectData(:,:,t);
        end
        AllData = structfun(@(x) permute(x,[2,1]),AllData,'UniformOutput',0);

        % Determine if independent-/dependent-samples analysis
        if isfield(Design, 'Between')
            StatsIdx = 'i';
        elseif isfield(Design, 'Within')
            StatsIdx = 'd';
        end

        % Run the statistical test
        [StatsResults,Cluster_Results]=Perm_Ttest(AllData,Design,TemplateEEG,...
            'N_Permutes',NPermut,'Pval',AlphaThresh,'root_folder',pwd,'TFCE',...
            TFCE,'StatsType',StatsIdx); 

        % Adding information to the Spect structure
        PermResults.Statistics = StatsResults;
        PermResults.Clusters = Cluster_Results;

        % Permutation threshold (e.g. 95% confidence interval)                     
    %     U = round((1-AlphaThresh)*NPermut); 
    %     MaxTFCE=sort(StatsResults.maxTFCE);
    %     PermResults.Cluster_Threshold= MaxTFCE(U);

        %% Plotting
        STUDY_Figures(STUDY,PermResults,SpectData,SpectFreqs,TemplateEEG.chanlocs,...
            'All','exportpath',[SavePath '\STUDY\' Date_Start '\'],'alphathresh',AlphaThresh);   

    % ANOVAs
    elseif length(STUDY.group)>1 && length(STUDY.condition)==1 || ...
            length(STUDY.group)>1 && length(STUDY.condition)>1
        Test = 'anova';

        % Building AllData structure
        FullSpectData = zeros(size(SpectData{1},1),size(SpectData{1},2), length(STUDY.datasetinfo));
        PosAll = zeros(length(STUDY.group),length(STUDY.condition));
        for t=1:length(STUDY.datasetinfo)
            PosGrp = find(ismember(STUDY.group,STUDY.datasetinfo(t).group));
            PosCnd = find(ismember(STUDY.condition,STUDY.datasetinfo(t).condition));
            PosAll(PosGrp,PosCnd) = PosAll(PosGrp,PosCnd) + 1;
            Idx = PosAll(PosGrp,PosCnd);
            FullSpectData(:,:,t) = SpectData{PosGrp,PosCnd}(:,:,Idx);
        end

    %     for t = 1:length(STUDY.group)
    %         GroupTemp = [GroupTemp;repmat(STUDY.group(t),size(SpectData{t,1},3),1)];
    % 
    %         FullSpectData = cat(3,FullSpectData,SpectData{t,1});
    %     end

        AllDataFields = strcat({STUDY.datasetinfo.group}',...
            {STUDY.datasetinfo.subject}','_',{STUDY.datasetinfo.condition}');
        for t = 1:length(AllDataFields)
            AllData.(AllDataFields{t}) = FullSpectData(:,:,t);
        end
        AllData = structfun(@(x) permute(x,[2,1]),AllData,'UniformOutput',0);

        % Run the statistical test
        [StatsResults,Cluster_Results]=Perm_ANOVA(AllData,Design,TemplateEEG,...
            'N_Permutes',NPermut,'Pval',AlphaThresh,'root_folder',pwd,'TFCE',...
            TFCE);

        % Adding information to the Spect structure
        PermResults.Statistics = StatsResults;
        PermResults.Clusters = Cluster_Results;
        PermResults.AlphaThresh = AlphaThresh;

        % Permutation threshold (e.g. 95% confidence interval)                     
    %     U = round((1-AlphaThresh)*NPermut); 
    %     MaxTFCE=sort(StatsResults.maxTFCE);
    %     PermResults.Cluster_Threshold= MaxTFCE(U);

        %% Plotting
        STUDY_Figures(STUDY,PermResults,SpectData,SpectFreqs,TemplateEEG.chanlocs,...
            'All','exportpath',[SavePath '\STUDY\' Date_Start '\'],'alphathresh',AlphaThresh);   

    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2.1 : STATISTICS on FREQUENCY BANDS (average of frequency bins)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Clear data
    clear PermResults

    % Update progress
    WaitBarApp.Message = 'Statistics on averaged frequency range';

    % For each frequency bands of interest
    for k=1:size(FreqData,1)

        % Update progress
        WaitBarApp.Value = (1+k)/(2+size(FreqData,1));

        % Finding frequency band specific ranges
        [~, LowBound] = min(abs(SpectFreqs-FreqRanges(k,1)));
        [~, HighBound] = min(abs(SpectFreqs-FreqRanges(k,2)));

        %% T-tests
        if strcmpi(Test,'t-test')

            % Averaging over the frequencies
            SpectDataFreq = cellfun(@(x) squeeze(mean(x(LowBound:HighBound,:,:),1)),SpectData,'UniformOutput',0);

            % Computing neighbouring channels 
            ChN = ept_ChN2(TemplateEEG.chanlocs, 0);

            % Virtually duplicating the channel dimensions 
            TestData = cell(1,length(SpectDataFreq));
            for t = 1:length(SpectDataFreq)
                TestData{t} = reshape(repmat(permute(SpectDataFreq{t},[2 1]),...
                    [size(SpectDataFreq{t},1),1]),[size(SpectDataFreq{t},2),...
                    size(SpectDataFreq{t},1),size(SpectDataFreq{t},1)]);
            end

            % Permutation test
            Results = ept_TFCE(TestData{1},TestData{2},TemplateEEG.chanlocs,'nPerm',...
                NPermut,'rSample', TemplateEEG.srate,'ChN', ChN,'flag_tfce',TFCE,...
                'type',StatsIdx);

            % Permutation threshold (e.g. 95% confidence interval)                     
    %         U = round((1-AlphaThresh)*NPermut); 
    %         MaxTFCE=sort(Results.maxTFCE);

            % Retrieving significant results
            PermResults(k).TFCE = Results;
    %         PermResults(k).TFCE.Threshold = MaxTFCE(U);
            PermResults(k).Cluster_Results = ept_calculateClusters(Results, ChN, AlphaThresh);

            %% Plotting the results   
            STUDY_Figures(STUDY,PermResults(k),SpectDataFreq,SpectFreqs,TemplateEEG.chanlocs,...
                'Bands','exportpath',[SavePath '\STUDY\' Date_Start '\'],'freqdata',...
                FreqData(k,:),'alphathresh',AlphaThresh);   

        %% ANOVAs
        elseif strcmpi(Test,'anova')

            % Datasets
            Temp = AllData.(AllDataFields{1});
            AllDataFreq = structfun(@(x) repmat(x,[size(x,1),1]),AllData,'UniformOutput',0);
            AllDataFreq = structfun(@(x) reshape(x,[size(Temp,1) size(Temp,1) size(Temp,2)]),AllDataFreq,'UniformOutput',0);
            % Averaging over the frequencies
            AllDataFreq = structfun(@(x) squeeze(mean(x(:,:,LowBound:HighBound),3)),AllDataFreq,'UniformOutput',0); 
            AllDataFreq = structfun(@(x) x',AllDataFreq,'UniformOutput',0); 
            % Averaging over the frequencies
            SpectDataFreq = cellfun(@(x) squeeze(mean(x(LowBound:HighBound,:,:),1)),SpectData,'UniformOutput',0);

            % Permutation test
            [Results,Cluster_Results]=Perm_ANOVA(AllDataFreq,Design,TemplateEEG,...
                'N_Permutes',NPermut,'Pval',AlphaThresh,'root_folder',pwd,'TFCE',...
                TFCE);

            % Retrieving significant results
            PermResults(k).TFCE = Results;
            PermResults(k).Cluster_Results = Cluster_Results;

            %% Plotting the results   
            STUDY_Figures(STUDY,PermResults(k),SpectDataFreq,SpectFreqs,TemplateEEG.chanlocs,...
                'Bands','exportpath',[SavePath '\STUDY\' Date_Start '\'],'freqdata',...
                FreqData(k,:),'alphathresh',AlphaThresh);   
        end
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 3 : STATISTICS over ALL FREQUENCIES (average of electrodes)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % There is a conceptual problem here with the neighbouring matrix. Since we
    % do not have a channel/source surface anymore, i don't know what to do.
    % Should all frequency bins be considered independent of each other?

    % Clear data
    clear PermResults

    % Update progress
    WaitBarApp.Message = 'Statistics on averaged electrodes';
    WaitBarApp.Value = (2+size(FreqData,1))/(2+size(FreqData,1));

    % Averaging over the channels
    SpectDataChan = cellfun(@(x) squeeze(mean(x,2)),SpectData,'UniformOutput',0);

    %% T-tests
    if strcmpi(Test,'t-test')

         % Virtually duplicating the channel dimensions    
        TestData{1} = reshape(repmat(permute(SpectDataChan{1},[2 1]),...
            [size(SpectDataChan{1},1),1]),[size(SpectDataChan{1},2),size(SpectDataChan{1},1),size(SpectDataChan{1},1)]);
        TestData{2} = reshape(repmat(permute(SpectDataChan{2},[2 1]),...
            [size(SpectDataChan{2},1),1]),[size(SpectDataChan{2},2),size(SpectDataChan{2},1),size(SpectDataChan{2},1)]);

        % Permutation test
        PermResults.TFCE = ept_TFCE(TestData{1},TestData{2},TemplateEEG.chanlocs,'nPerm',...
            NPermut,'rSample', TemplateEEG.srate,'flag_tfce',TFCE,'flag_ft',1,'type',StatsIdx);
        
        % Using EEGLAB functions + FDR correction for multiple comparisons
        [stats, df, pvals] = statcond(SpectDataChan','paired',fastif(strcmpi(StatsIdx,'d'),'on','off'),...
             'method','perm','naccu',NPermut,'verbose','off','alpha',AlphaThresh);
        [p_masked, ~, ~, pvals_FDR]=fdr_bh(pvals,AlphaThresh);
        
        % Using Fieldtrip statistics + max cluster correction (ERRORS)
        % See: https://github.com/sccn/eeglab/issues/184
%         
%         if strcmpi(StatsIdx,'d')
%             [pcond, ~, ~, statscond] = ...
%             std_stat(SpectDataChan, 'condstats','on', 'fieldtripnaccu',NPermut,'fieldtripmethod',...
%             'montecarlo','fieldtripmcorrect','max','fieldtripalpha',AlphaThresh,'mode','fieldtrip');
%         else
%             [~, pgroup, ~, ~, statsgroup] = ...
%             std_stat(SpectDataChan', 'groupstats','on','fieldtripnaccu',NPermut,'fieldtripmethod',...
%             'montecarlo','fieldtripmcorrect','max','fieldtripalpha',AlphaThresh,'mode','fieldtrip');
%         end
%         [stats, df, pvals] = statcondfieldtrip(SpectDataChan','paired',fastif(strcmpi(StatsIdx,'d'),'on','off'),...
%         'method','permutation','naccu',NPermut,'alpha',AlphaThresh,'fieldtripmcorrect','max',...
%         'avgoverchan','yes','avgovertime','yes');

        % Permutation threshold (e.g. 95% confidence interval)                     
    %     U = round((1-AlphaThresh)*NPermut); 
    %     MaxTFCE=sort(PermResults.TFCE.maxTFCE);

        % Finding clusters
        PermResults.AlphaThresh = AlphaThresh;
        PermResults.Cluster_Results = ept_calculateClusters(PermResults.TFCE, ChN, AlphaThresh);

        %% Plotting the results
        STUDY_Figures(STUDY,PermResults,SpectDataChan,SpectFreqs,TemplateEEG.chanlocs,...
            'AvgFreqs','exportpath',[SavePath '\STUDY\' Date_Start '\'],'freqdata',...
            FreqData(k,:),'alphathresh',AlphaThresh);   

    %% ANOVAs
    elseif strcmpi(Test,'anova')

        % Datasets
        AllDataChan = structfun(@(x) permute(x,[2,1]),AllData,'UniformOutput',0);
        AllDataChan = structfun(@(x) repmat(x,[size(x,1),1]),AllDataChan,'UniformOutput',0);
        AllDataChan = structfun(@(x) reshape(x,[size(Temp,2) size(Temp,2) size(Temp,1)]),AllDataChan,'UniformOutput',0);
        % Averaging over the channels
        AllDataChan = structfun(@(x) squeeze(mean(x,3)),AllDataChan,'UniformOutput',0);
        AllDataChan = structfun(@(x) x',AllDataChan,'UniformOutput',0); % UNSURE OF THIS !!! 

        % Permutation test
        [Results,Cluster_Results]=Perm_ANOVA(AllDataChan,Design,TemplateEEG,...
            'N_Permutes',NPermut,'Pval',AlphaThresh,'root_folder',pwd,'TFCE',...
            TFCE);

        % Retrieving significant results
        PermResults.TFCE = Results;
        PermResults.Cluster_Results = Cluster_Results;

        % COMPARING WITH MATLAB AOV FUNCTION
    %     for m=1:size(SpectDataChan{1},1)
    %         AllPlotData = [];
    %         GrpLab = {};
    %         for p=1:length(STUDY.group)
    %             Dat = SpectDataChan{p};
    %             GrpLab = [GrpLab; repmat(STUDY.group(p),[size(Dat,2),1])];
    %             AllPlotData = [AllPlotData;Dat(m,:)'];
    %         end
    %         Pvalues(m) = anova1(AllPlotData,GrpLab);
    %         close all
    %     end

        %% Plotting the results
        STUDY_Figures(STUDY,PermResults,SpectDataChan,SpectFreqs,TemplateEEG.chanlocs,...
            'AvgFreqs','exportpath',[SavePath '\STUDY\' Date_Start '\'],'freqdata',...
            FreqData(k,:),'alphathresh',AlphaThresh);   
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MICROSTATES ANALYSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I STOPPED HERE FOR THE 2x2 mixed ANOVAs tests !!!! 

% Calling the function
if strcmpi(MicroStatesSwitch,'Yes')
    MicroStatesParam = [MicroStatesParam;{Algo};{MicroStatesComputGroups}];
    LogMST = MicroStates(STUDY,'ExportPath',[SavePath '\Exports\' Date_Start '\MicroStatesSegment\'],...
       'ExcelDirectory',[ExcelDirectory Date_Start],'MicroStatesParam',MicroStatesParam,...
       'SamplingRate',MSTSamplingRate,'LowPass',MSTLowPass);  % ,'GUI',WaitBarApp
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      IC CLUSTERS SOURCE LOCALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calling the function
if strcmpi(ICclusteringSwitch,'Yes') && strcmpi(ICAexist,'Yes')
    % Precompute components measures
    [STUDY ALLEEG] = std_precomp(STUDY, ALLEEG,'components','spec','on','scalp','on',...
    'recompute','on','specparams',{'specmode' 'psd','logtrials','off'});

    % Create preclustering array
    [STUDY,ALLEEG] = std_preclust(STUDY,ALLEEG,[],{ 'spec'  'npca' 10 'norm' 1 ...
        'weight' 1},{ 'scalp' 'npca' 10 'norm' 1 'weight' 1 'abso' 1 },...
        { 'dipoles' 'norm' 1 'weight' 10 });
    
    LogICCLust = ICClustLocalisation(STUDY,ALLEEG,'ExportPath',[SavePath '\Exports\' Date_Start '\ICClust\'],...
        'ExcelDirectory',[ExcelDirectory Date_Start],'AllParameters',ICClustParam); % ,'GUI',WaitBarApp
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            RUNNING THE LOG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sum(~SavePath)<1
    LOG(SavePath,Date_Start,'LogMST',LogMST,'LogIC',LogICCLust,'ErrorStudy',ErrorLog)
else
    LOG(CurrentPWD,Date_Start,'LogMST',LogMST,'LogIC',LogICCLust,'ErrorStudy',ErrorLog)
end

end