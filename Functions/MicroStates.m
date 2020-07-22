function Log = MicroStates(STUDY,varargin)

% This function computes semi-automated microstates analyses based on
% data from the GROUPSTUDY script and the MST1.0 EEGLab toolbox:
% https://github.com/atpoulsen/Microstate-EEGlab-toolbox

% Usage:
%    >> Log = MicroStates(STUDY,varargin);
%
% Inputs (mandatory):
%   STUDY            = EEGlab STUDY structure

% Inputs (Optional):
%   ExportPath       = Path where to export results.
%   ExcelDirectory   = Path containing the excel template.
%   LowPass          = Applied a low pass filter (in Hz).
%   Sampling Rate    = Resamples the data for faster computation.
%   MicroStatesParam = Parameters imported from StudyGUI.mlapp. 

% Outputs:
%   Log              = Structure containing parameters information for each
%                      version that the user ran.

% Author: Corentin Wicht, LCNS, 2019
% - corentin.wicht@unifr.ch
% - https://github.com/CorentinWicht

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

%% Set Defaults
ExcelDirectory = pwd;
ExportPath = pwd;

% Parameters
FitMeasList = {'GEV', 'CV', 'W', 'KL', 'KL_nrm'};
AlgoList = {'kmeans','modkmeans','aahc','taahc','varmicro'};
NormIndiv = 0;
MinPeakDist = 10;
Npeaks = 1000;
GFPthresh = 1;
Algo = 'modkmeans';
Nmicrostates = 2:8;
NormAll = 0;
Nrepet = 50;
MaxIter = 1000;
Thresh = 1e-06;
FitMeas = 'CV';
TimeRange = [];
Determin = 0;
Polarity = 0;
Sig2_0 = [];
P0 = 0;
ComputGroups = 'off';
MSTLowPass = [];
MSTSamplingRate = [];

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
            case 'exportpath'
                ExportPath  = Value;
            case 'exceldirectory'
                ExcelDirectory  = Value;
            case 'lowpass'
                MSTLowPass = Value;
            case 'samplingrate'
                MSTSamplingRate = Value;
            case 'microstatesparam'
                MicroStatesParam = Value;
                
                % Parameters
                Algo = MicroStatesParam{end-1};
                ComputGroups = MicroStatesParam{end};
                NormIndiv = strcmpi(MicroStatesParam{1},'Y');
                MinPeakDist = str2double(MicroStatesParam{2});
                Npeaks = str2double(MicroStatesParam{3});
                GFPthresh = str2double(MicroStatesParam{4});
                Nmicrostates = str2num(MicroStatesParam{5});
                NormAll = strcmpi(MicroStatesParam{6},'Y');
                TimeRange = str2num(MicroStatesParam{7});
                
                % Depending on selected algorithm 
                if strcmpi(Algo,'kmeans')
                    Nrepet = str2double(MicroStatesParam{8});
                    MaxIter = str2double(MicroStatesParam{9});
                    
                elseif strcmpi(Algo,'modkmeans')
                    Nrepet = str2double(MicroStatesParam{8});
                    MaxIter = str2double(MicroStatesParam{9});
                    Thresh = str2double(MicroStatesParam{10});
                    FitMeas = FitMeasList{str2double(MicroStatesParam{11})};
                    
                elseif strcmpi(Algo,'aahc') || strcmpi(Algo,'taahc') 
                    Determin = strcmpi(MicroStatesParam{8},'Y');
                    Polarity = strcmpi(MicroStatesParam{9},'Y');
                    
                elseif strcmpi(Algo,'varmicro')
                    Nrepet = str2double(MicroStatesParam{8});
                    MaxIter = str2double(MicroStatesParam{9});
                    Thresh = str2double(MicroStatesParam{10});
                    Sig2_0 = str2num(MicroStatesParam{11});
                    P0 = str2double(MicroStatesParam{12}); 
                end
                
            otherwise
                display (['Unknown parameter setting: ' Param])
        end
    end
end

%% 
    % TO DO:
    % 1) Should include possibility to re-run w/ different parameters (like
    % ICA).
    % 4) Implement the synthax output (if doesn't exist yet) ?
    
    % Initialize the while loop enabling to re-run analysis w/ other
    % parameters
    Restart=0;
    Increment=1;

    % Initialize EEGLab
    eeglab
    close gcf
    
    % If MATLAB is out of memory, use single precision instead
    try
        
        % EEGLAB options
        pop_editoptions('option_storedisk', 0, 'option_single', 0);

        % Load STUDY
        [STUDY ALLEEG] = pop_loadstudy('filename', STUDY.filename,'filepath', STUDY.filepath);
    catch
        
        warning(['Warning: Matlab is out of memory. Trying now to load the STUDY data' ...
            'in single-precision. If it fails, increase your system swap space, ' ...
            'i.e. use harddrive to compensate for lack of physical memory space (RAM).'])
        
        % EEGLAB options (set to single precision)
        pop_editoptions('option_storedisk', 0, 'option_single', 1);
        
        % Load STUDY
        [STUDY ALLEEG] = pop_loadstudy('filename', STUDY.filename,'filepath', STUDY.filepath);

    end  
    
    % Apply LowPass and/or resampling(optional)
    if  ~isempty(MSTLowPass) || ~isempty(MSTSamplingRate)
        for k=1:length(ALLEEG) 
            % Load each file
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, ALLEEG(1), k,'retrieve',k,'study',0,...
                'setname',[ALLEEG(1).setname '_Processed']); 

            % High-pass/Low-pass filter   
            if ~isempty(MSTLowPass) 
                EEG = pop_eegfiltnew(EEG,[],MSTLowPass); 
                [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            end
            
            % Resampling
            if ~isempty(MSTSamplingRate)

                EEG = pop_resample(EEG, MSTSamplingRate);
            end
            
            % Store the new data
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        end
    end
    
    % If the user requires the analysis to be performed separately for each
    % level of the between-subject factor
    % SHOULD ALSO ADAPT THE EXCEL FILE OUTPUT!!!! 
    AllGroups = {STUDY.datasetinfo.group};
    if strcmpi(ComputGroups,'on')
        
        % Number of runs
        Runs = length(STUDY.group);    
        
        % Find subject list corresponding to each groups
        for k = 1:length(STUDY.group)
            GroupSubjects{k} = ismember(AllGroups,STUDY.group{k});
        end
    else
        Runs = 1;
    end
    
    % Until the user is satisfied with the results
    while Restart<1
        %% CLUSTERING  
        % Clustering or each group (optional)
        % THE IDEA IS TO CLUSTER A THIRD TIME THE TWO PREVIOUS CLUSTERING.
        % HENCE NEED TO FIND HOW TO COMBINE BOTH CLUSTERING IN
        % ALLEEG(end)!!
        % TO DO:
        % 1) TO HAVE BOTH FIGURES AT THE SAME TIME ON THE SCREEN !!!
        % 2) COMPARE THE FACT OF MERGING DATASETS TO CARTOOL OUTPUTS, IT
        % SEEMS TO WORK BUT I AM UNSURE IT IS THE CORRECT WAY !

        for f=1:Runs

            % Only if separated by groups
            if Runs > 1
                if f == 1
                    ALLEEGTEMP = ALLEEG;
                    ALLEEG = ALLEEGTEMP(GroupSubjects{f});
                else
                    ALLEEG = ALLEEGTEMP(GroupSubjects{f});
                end

                % Saving subjects/conditions list
                SubjList = cell(length(ALLEEG),1);
                CondList = cell(length(ALLEEG),1);
                GrpList = cell(length(ALLEEG),1);
                Pos = 1;
                for p=find(GroupSubjects{f})
                    TempPath = strsplit(STUDY.datasetinfo(p).filepath,'\');
                    SubjList(Pos) = TempPath(end); 
                    CondList{Pos} = STUDY.datasetinfo(p).condition;
                    GrpList{Pos} = STUDY.datasetinfo(p).group;
                    Pos = Pos + 1;
                end
            else

                % Saving subjects/conditions list
                SubjList = cell(length(STUDY.datasetinfo),1);
                CondList = cell(length(STUDY.datasetinfo),1);
                GrpList = cell(length(STUDY.datasetinfo),1);
                for p=1:length(STUDY.datasetinfo)
                    TempPath = strsplit(STUDY.datasetinfo(p).filepath,'\');
                    SubjList(p) = TempPath(end); 
                    CondList{p} = STUDY.datasetinfo(p).condition;
                    GrpList{p} = STUDY.datasetinfo(p).group;
                end
            end

            EEG = [];
            % Data selection
            [~, ALLEEG] = pop_micro_selectdata(EEG, ALLEEG, 'datatype', 'spontaneous',...
                'avgref', 1,'normalise', NormIndiv, 'MinPeakDist', MinPeakDist, 'Npeaks', Npeaks, ...
                'GFPthresh', GFPthresh, 'dataset_idx', 1:length(ALLEEG));

            % Save for further use in loop
            SavedExportPath = ExportPath;

            % select the "GFPpeak" dataset and make it the active set
            % Free RAM memory space by removing unused datasets
            ALLEEG = ALLEEG(end);

            % Creating export folder according to version of analysis
            ExportPath = [SavedExportPath sprintf('Version%d',Increment) '\'];
            mkdir (ExportPath)
            mkdir ([ExportPath 'FIG'])

            % Microstate segmentation (heavy computation)
            if strcmpi(Algo,'modkmeans')
                ALLEEG = pop_micro_segment(ALLEEG, 'algorithm', 'modkmeans', 'sorting',...
                    'Global explained variance', 'Nmicrostates', Nmicrostates, 'verbose', 1, ...
                    'normalise', NormAll, 'Nrepetitions', Nrepet, 'max_iterations', MaxIter, ...
                    'threshold', Thresh, 'fitmeas', FitMeas,'optimised',1);

            elseif strcmpi(Algo,'kmeans')
                  ALLEEG = pop_micro_segment(ALLEEG, 'algorithm', 'kmeans', 'sorting',...
                    'Global explained variance', 'Nmicrostates', Nmicrostates, 'verbose', 1, ...
                    'normalise', NormAll, 'Nrepetitions', Nrepet, 'max_iterations', MaxIter);

            elseif strcmpi (Algo,'aahc') || strcmpi (Algo,'taahc')
                  ALLEEG = pop_micro_segment(ALLEEG, 'algorithm', Algo, 'sorting',...
                    'Global explained variance', 'Nmicrostates', Nmicrostates, 'verbose', 1, ...
                    'normalise', NormAll, 'determinism', Determin, 'polarity', Polarity);

            elseif strcmpi(Algo,'varmicro')
                ALLEEG = pop_micro_segment(ALLEEG, 'algorithm', 'varmicro', 'sorting',...
                    'Global explained variance', 'Nmicrostates', Nmicrostates, 'verbose', 1, ...
                    'normalise', NormAll, 'Nrepetitions', Nrepet, 'max_iterations', MaxIter, ...
                    'threshold', Thresh, 'sig2_0', Sig2_0,'p0',P0);
            end

            % Plot microstate prototype topographies
            figure;MicroPlotTopo(ALLEEG, 'plot_range', []);
            % Holds the figure until inspection is over
%             Fig=msgbox('THE CODE WILL CONTINUE ONCE YOU PRESS OK','WAIT','warn'); 
%             uiwait(Fig);
            movegui(gcf,'northwest')

            % microstates fit measures
            Measures = {'CV', 'GEV','W', 'KL', 'KL_nrm'};
            ALLEEG = pop_micro_selectNmicro(ALLEEG,'Measures',Measures,'do_subplots',1);
            
            % Saves the topographies in export folder
            if Runs > 1
                SaveFigures(gcf,[STUDY.filepath '\MicroStates\' sprintf('ClusterTopo_%s_%s_v%d',STUDY.group{f},...
                    STUDY.name,Increment)],'w','bmp',1);
            else
                SaveFigures(gcf,[STUDY.filepath '\MicroStates\' sprintf('ClusterTopo_%s_v%d',STUDY.name,Increment)],'w','bmp',1);
            end

            % Replotting the fit measures to export them
            [~, plot_idx] = ismember(Nmicrostates, ALLEEG.microstate.algorithm_settings.Nmicrostates);
            figure('Units', 'normalized','position',[.2 .2 .6 .6]);
            for m=1:length(Measures)
                subplot(length(Measures),1,m)
                MicroPlotFitmeas(ALLEEG.microstate.Res,Measures(m),Nmicrostates,plot_idx)
            end

            % Saves the topographies in export folder
            if Runs > 1
                SaveFigures(gcf,[STUDY.filepath '\MicroStates\' sprintf('FitMeasures_%s_%s_v%d',STUDY.group{f},...
                    STUDY.name,Increment)],'w','bmp');
            else
                SaveFigures(gcf,[STUDY.filepath '\MicroStates\' sprintf('FitMeasures_%s_v%d',STUDY.name,Increment)],'w','bmp');
            end

            % Saving the clustering dataset
            TempEEG(f) = ALLEEG;
        end
        
        % Saving subjects/conditions list
        SubjList = cell(length(STUDY.datasetinfo),1);
        CondList = cell(length(STUDY.datasetinfo),1);
        GrpList = cell(length(STUDY.datasetinfo),1);
        for p=1:length(STUDY.datasetinfo)
            TempPath = strsplit(STUDY.datasetinfo(p).filepath,'\');
            SubjList(p) = TempPath(end); 
            CondList{p} = STUDY.datasetinfo(p).condition;
            GrpList{p} = STUDY.datasetinfo(p).group;
        end

        % If separated by groups, concatenate the clustering datasets
        if Runs > 1
            
            % Merge datasets
            TempEEG = pop_mergeset(TempEEG,1:Runs);
            
            % Re-run the microstate segmentation on the merged sets
            if strcmpi(Algo,'modkmeans')
                TempEEG = pop_micro_segment(TempEEG, 'algorithm', 'modkmeans', 'sorting',...
                    'Global explained variance', 'Nmicrostates', Nmicrostates, 'verbose', 1, ...
                    'normalise', NormAll, 'Nrepetitions', Nrepet, 'max_iterations', MaxIter, ...
                    'threshold', Thresh, 'fitmeas', FitMeas,'optimised',1);

            elseif strcmpi(Algo,'kmeans')
                  TempEEG = pop_micro_segment(TempEEG, 'algorithm', 'kmeans', 'sorting',...
                    'Global explained variance', 'Nmicrostates', Nmicrostates, 'verbose', 1, ...
                    'normalise', NormAll, 'Nrepetitions', Nrepet, 'max_iterations', MaxIter);

            elseif strcmpi (Algo,'aahc') || strcmpi (Algo,'taahc')
                  TempEEG = pop_micro_segment(TempEEG, 'algorithm', Algo, 'sorting',...
                    'Global explained variance', 'Nmicrostates', Nmicrostates, 'verbose', 1, ...
                    'normalise', NormAll, 'determinism', Determin, 'polarity', Polarity);

            elseif strcmpi(Algo,'varmicro')
                TempEEG = pop_micro_segment(TempEEG, 'algorithm', 'varmicro', 'sorting',...
                    'Global explained variance', 'Nmicrostates', Nmicrostates, 'verbose', 1, ...
                    'normalise', NormAll, 'Nrepetitions', Nrepet, 'max_iterations', MaxIter, ...
                    'threshold', Thresh, 'sig2_0', Sig2_0,'p0',P0);
            end
            
            % Store the final segmentation
            ALLEEG = TempEEG;
        end
        
        % I STOPPED HERE (14.11.2019)
              
        % Creating the excel template
        NumMST = ALLEEG(end).microstate.Res.K_act;
        HeadersMicroStates = sprintf('TempVar%d_',1:1+6*NumMST);
        HeadersMicroStates = [HeadersMicroStates sprintf('TempVar%d_',2+6*NumMST:(2+6*NumMST)+NumMST*(NumMST-1)-1)]; % for TP
        HeadersMicroStates = [{'Participants';'Conditions';'Groups'};strsplit(HeadersMicroStates,'_')];
        HeadersMicroStates(end) = [];
        TemplateTable = cell2table([SubjList CondList GrpList repmat({[]},[length(STUDY.datasetinfo),length(HeadersMicroStates)-3])]); 
        TemplateTable.Properties.VariableNames = HeadersMicroStates;

        % Save outputs for the LOG
        Log.(sprintf('V%d',Increment)).AlgoList = NormIndiv;
        Log.(sprintf('V%d',Increment)).MinPeakDist = MinPeakDist;
        Log.(sprintf('V%d',Increment)).Npeaks = Npeaks;
        Log.(sprintf('V%d',Increment)).GFPthresh = GFPthresh;
        Log.(sprintf('V%d',Increment)).Algo = Algo;
        Log.(sprintf('V%d',Increment)).Nmicrostates = Nmicrostates;
        Log.(sprintf('V%d',Increment)).NormAll = NormAll;
        Log.(sprintf('V%d',Increment)).Nrepet = Nrepet;
        Log.(sprintf('V%d',Increment)).MaxIter = MaxIter;
        Log.(sprintf('V%d',Increment)).Thresh = Thresh;
        Log.(sprintf('V%d',Increment)).FitMeas = FitMeas;
        Log.(sprintf('V%d',Increment)).TimeRange = TimeRange;
        Log.(sprintf('V%d',Increment)).FinalNumMST = NumMST;
        Log.(sprintf('V%d',Increment)).ComputGroups = ComputGroups;
        
        % Reload the study (to integrate all datasets again)
        ALLEEG = [];
        pop_editoptions( 'option_storedisk', 1);
        [STUDY,ALLEEG] = pop_loadstudy('filename', STUDY.filename,'filepath', STUDY.filepath);
        CURRENTSET = 1;
        EEG = ALLEEG(1);

        % Add new ALLEEG empty field (microstates)
        for k=1:length(ALLEEG)
           ALLEEG(k).microstate = []; 
        end

        % Add the avgEEG to the ALLEEG structure (at the end)
        ALLEEG = [ALLEEG,TempEEG];

        %% BACK-FITTING

        % Import microstate prototypes from other dataset to the datasets that 
        % should be back-fitted
        if exist('StatsMat','var')
            clear StatsMat
        end
        for i = 1:length(ALLEEG)-1

            % Import data
            fprintf('Importing prototypes and backfitting for dataset %i\n',i)
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',i,'study',0,...
                'setname',[ALLEEG(1).setname '_Processed']);
            EEG = pop_micro_import_proto(EEG, ALLEEG, length(ALLEEG));

            % Back-fit microstates on EEG
            EEG = pop_micro_fit( EEG, 'polarity', 0); % ISNT POLARITY SUPPOSED TO BE IN THE PROMPT?????

            % Temporally smooth microstates labels
            EEG = pop_micro_smooth( EEG, 'label_type', 'backfit', 'smooth_type',...
                'reject segments', 'minTime', 30, 'polarity', 0);

            % Calculate microstate statistics
            EEG = pop_micro_stats( EEG, 'label_type', 'backfit', 'polarity', 0);
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                     
            % Storing the transition probabilitie
            Pos = 1;TP=zeros(1,nnz(EEG.microstate.stats.TP));TPLabel=cell(1,nnz(EEG.microstate.stats.TP));
            for m=1:size(EEG.microstate.stats.TP,1)
                for p=find(EEG.microstate.stats.TP(m,:)~=0)
                   TP(Pos) = EEG.microstate.stats.TP(m,p);
                   TPLabel{Pos} = sprintf('M%dtoM%d',m,p);
                   Pos = Pos + 1;
                end
            end

            % Creating the table contaning stats results
            StatsField = fieldnames(EEG.microstate.stats);
            TemplateTable(i,4:end) = array2table(num2cell([EEG.microstate.stats.(StatsField{1}) EEG.microstate.stats.(StatsField{2}) ...
                EEG.microstate.stats.(StatsField{3}) EEG.microstate.stats.(StatsField{4}) ...
                EEG.microstate.stats.(StatsField{5}) EEG.microstate.stats.(StatsField{6}) ...
                EEG.microstate.stats.(StatsField{7}) TP]));

            % Current subject
            CurrSubj = table2array(TemplateTable(i,1));

            % Current condition
            CurrCond = table2array(TemplateTable(i,2));

            % Saving each subject result as a .mat file in Export folder
            TempMST = EEG.microstate;
            save([ExportPath sprintf('MicroStateRes_S%s_%s_v%d.mat',...
                CurrSubj{:},CurrCond{:},Increment)],'TempMST');

            % Plotting GFP for each file
            figure;MicroPlotSegments( EEG, 'label_type', 'backfit', ...
                'plotsegnos', 'first', 'plottopos', 1, 'plot_time',TimeRange);

            % Saving as .fig for further use in Matlab
            savefig(gcf,[ExportPath 'FIG\' sprintf('GFP_S%d_%s_v%d',...
                CurrSubj{:},CurrCond{:},Increment)]);

            % Saving each subject's GFP figure
            SaveFigures(gcf,[ExportPath sprintf('GFP_S%d_%s_v%d',...
                CurrSubj{:},CurrCond{:},Increment)],'w','bmp');
        end

        % Exporting stats results in excel template
        % Header
        StatsField = fieldnames(EEG.microstate.stats);
        StatsHeader = {}; 
        for i=1:length(StatsField)-5 % Export only half in excel, the rest in .mat
            if size(EEG.microstate.stats.(StatsField{i}),2) == size(EEG.microstate.stats.TP,1)
                for k = 1:size(EEG.microstate.stats.TP,1)
                    StatsHeader = [StatsHeader strcat(StatsField(i),sprintf('_Map%d',k))];
                end
            else
                StatsHeader = [StatsHeader StatsField(i)];
            end
        end
        TemplateTable.Properties.VariableNames(4:end) = [StatsHeader TPLabel];

        % Write to excel template
        writetable(TemplateTable,[ExcelDirectory '\' sprintf('MicroStates_v%d.xlsx',Increment)],'Sheet','Data');    

        % Satisfied with current results`?
        ConfirmAnalyses = questdlg(['Are you happy with the MicroStates results?' newline ...
            'If no, the code will start again'], 'Restart analyses ?', 'No','Yes','Yes');
        if strcmpi(ConfirmAnalyses,'Yes')

            % End the loop
            Restart = 1;

        else
            Increment = Increment + 1;
            disp('Restarting the analyses with different parameters')

            % Select algorithm
            [~,Algo] = listdlg('ListString',AlgoList);

            % Prompting for new parameters
            % Initialize the prompt 
            PromptMicroStates={'Normalize each dataset by average channel SD within each datatset [Y/N] ?',...
            'Minimal distance between GFP peaks',...
            'Number of GFP peaks maximally extracted per dataset',...
            'Discards maps with GFP exceeding [X] times the SD of the GFPs of all maps',...
            'Number of microstates to cluster (range)',...
            'Normalize across all datasets by average channel SD [Y/N] ?',...
            'Time range (is ms, separated by a space) to plot individual microstate segments over the GFP'};

            % Default parameters
            MicroStatesParam = {'N','10','1000','1','2:8','N', ''};

            % Prompt changes depending on selected algorithm
            if strcmpi(Algo,'kmeans')
                 % Default parameters
                MicroStatesParam = [MicroStatesParam {'50','1000'}];

                % Prompt
                PromptMicroStates = [PromptMicroStates ...
                {'Number of random initializations (set it high!)',...
                'Maximum number of iterations for all restarts'}];

            elseif strcmpi(Algo,'modkmeans')
                % Default parameters
                MicroStatesParam = [MicroStatesParam {'50','1000','1e-06','2'}];

                % Prompt
                PromptMicroStates = [PromptMicroStates ...
                {'Number of random initializations (set it high!)',...
                'Maximum number of iterations for all restarts',...
                'Threshold of convergence based on relative change in noise variance',...
                ['Measure of fit, criterion for selecting best segmentation:' newline ...
                '1) Global explained variance (GEV)' newline ...
                '2) Cross-validation (CV)' newline ...
                '3) Dispersion (W)' newline ...
                '4) Krzanowski-Lai (KL)' newline ...
                '5) Normalised Krzanowski-Lai (KLnrm)']}];

            elseif strcmpi(Algo,'aahc') || strcmpi(Algo,'taahc')
                % Default parameters
                MicroStatesParam = [MicroStatesParam {'N','N'}];

                % Prompt
                PromptMicroStates = [PromptMicroStates ...
                {['TAAHC initialisation scheme for making the clustering determinate.' ...
                'Initialises by so every cluster consists of two samples, by agglomerating the most correlated samples [Y/N]'],...
                'Account for polarity? If set to N, the sign of correlation is ignored [Y/N]'}];

            elseif strcmpi(Algo,'varmicro')
                % Default parameters
                MicroStatesParam = [MicroStatesParam {'50','1000','1e-06','','0'}];

                % Prompt
                PromptMicroStates = [PromptMicroStates ...
                {'Number of random initializations (set it high!)',...
                'Maximum number of iterations for all restarts',...
                'Threshold of convergence based on relative change in noise variance',...
                'Prior variance of activations (default: average EEG channel variance)'...
                'Probability for having same microstate as last timepoint. Setting to zero turns smoothing off'}];

            end

            % Running prompt
            MicroStatesParam=inputdlg(PromptMicroStates,'MicroStates (MST1.0 toolbox)',1,MicroStatesParam);  

            % Storing new parameters
            NormIndiv = strcmpi(MicroStatesParam{1},'Y');
            MinPeakDist = str2double(MicroStatesParam{2});
            Npeaks = str2double(MicroStatesParam{3});
            GFPthresh = str2double(MicroStatesParam{4});
            Nmicrostates = str2num(MicroStatesParam{5});
            NormAll = strcmpi(MicroStatesParam{6},'Y');
            TimeRange = str2num(MicroStatesParam{7});

            % Depending on selected algorithm 
            if strcmpi(Algo,'kmeans')
                Nrepet = str2double(MicroStatesParam{8});
                MaxIter = str2double(MicroStatesParam{9});

            elseif strcmpi(Algo,'modkmeans')
                Nrepet = str2double(MicroStatesParam{8});
                MaxIter = str2double(MicroStatesParam{9});
                Thresh = str2double(MicroStatesParam{10});
                FitMeas = FitMeasList{str2double(MicroStatesParam{11})};

            elseif strcmpi(Algo,'aahc') || strcmpi(Algo,'taahc') 
                Determin = strcmpi(MicroStatesParam{8},'Y');
                Polarity = strcmpi(MicroStatesParam{9},'Y');

            elseif strcmpi(Algo,'varmicro')
                Nrepet = str2double(MicroStatesParam{8});
                MaxIter = str2double(MicroStatesParam{9});
                Thresh = str2double(MicroStatesParam{10});
                Sig2_0 = str2num(MicroStatesParam{11});
                P0 = str2double(MicroStatesParam{12}); 
            end
        end
    end
end