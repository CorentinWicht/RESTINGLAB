function STUDY_Figures(STUDY,PermResults,SpectData,SpectFreqs,ChanLocs,PlotType,varargin)
% This function generates graphical outputs for the GROUPSTUDY script. The
% data should come from EEGLAB study structure.

% Usage:
%    >> STUDY_Figures(STUDY,PermResults,SpectData,SpectFreqs,ChanLocs,PlotType,varargin);
%
% Inputs (mandatory):
%   STUDY      = EEGLab STUDY (group analysis) data structure 
%   PermResults = Structure containing the results of the permutation
%                 statistics and the identification of clusters using TFCE.
%                 Results are extracted from the ept_TFCE-MATLAB toolbox.
%   SpectData  = Cell matrix containing for each group (columns), condition
%               (rows) PSD computed using the function "std_precomp" and
%                retrieved as the 2nd output of "std_readspec"
%   SpectFreqs = Array of PSD frequencies retrieved as 3rd output of 
%               "std_readspec"
%   ChanLocs   = EEGLAB structure of channel locations (EEG.chanlocs)
%   PlotType   = String indicating type of plot to compute ('all', 'bands'
%                or 'avgfreqs'). 

% Inputs (Optional):
%   FreqData   = Cell array containing in the first column a string with
%                the frequency band name (e.g. Delta) and second column 
%                a string with the frequency range (e.g. '2 4'). Default is
%                []. Careful, this is mandatory for PlotType 'bands' and
%                'avgfreqs'.
%   ExportPath  = String indicating the path where to save figures. Default
%                 is current folder [pwd].
%   AlphaThresh = Alpha significance threshold. Default is [0.05].


% Author: Corentin Wicht, LCNS, 2019
% - corentin.wicht@unifr.ch
% - https://github.com/CorentinWicht

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

%% Set Defaults
FreqData    = []; % default for frequency bins/ranges
ExportPath  = pwd; % default export path
AlphaThresh = 0.05; % default threshold for significance

% Initialize variables
ChannelsLabels = {ChanLocs.labels};
XTicks = 5:5:length(round(SpectFreqs));
if strcmpi(STUDY.SpecMode,'fft')
    Units = '10*Log10(\muV^2/Hz)';
else
    Units = '10*Log10(\muV^2)'; 
end
 
% Process Secondary Arguments
if nargin > 7
    for i = 1:2:length(varargin)
        Param = varargin{i};
        Value = varargin{i+1};
        if ~ischar(Param)
          error('Flag arguments must be strings')
        end
        Param = lower(Param);
    
        switch Param
            case 'freqdata'
                FreqData  = Value;
            case 'exportpath'
                ExportPath  = Value;
            case 'alphathresh'
                AlphaThresh = Value;
            otherwise
                display (['Unknown parameter setting: ' Param])
        end
    end
end

if strcmpi(PlotType,'all')
    %% STEP 1 : ALL DATA (channels X frequencies)
    
    % Settings
    PlotData = cellfun(@(x) squeeze(mean(x,3))',SpectData,'UniformOutput',0);
    MinColor = min(min(cellfun(@(x) min(min(x)),PlotData)));
    MaxColor = max(max(cellfun(@(x) max(max(x)),PlotData)));
    NPlots = size(PlotData,1)*size(PlotData,2);
    YTicks = 1:2:length(ChannelsLabels);
    TypeStats = {['Main-' STUDY.design.variable(2).label],...
        ['Main-' STUDY.design.variable(1).label],'Interaction'};
    
    % PARENT FIGURE
    figure; 
    annotation('textbox',[.42 .5 .5 .5],'String','All frequencies/channels',...
        'FitBoxToText','on','FontSize',20,'EdgeColor','w');

    % Group/Condition-specific data
    Pos = 0;
    for k=1:size(PlotData,1) % Groups
        for p = 1:size(PlotData,2) % Conditions
            Pos = Pos + 1;
            subplot(2,NPlots,Pos);
            imagesc(PlotData{k,p});caxis([MinColor MaxColor]); colormap(gca,'jet(50)');
            set(gca,'ydir','norm');xlabel('Frequencies'), ylabel('Channels')
            title([strrep(STUDY.group{k},'_','-') '-' STUDY.condition{p}]);
            set(gca,'XTick',XTicks,'XTickLabel',round(SpectFreqs(5:5:end)));
            set(gca,'YTick',YTicks,'YTickLabel',ChannelsLabels(1:2:end));
                    
            % contours
            if isfield(PermResults.Statistics,'Obs')
                for m=1:length(PermResults.Clusters)
                   hold on; contour(PermResults.Clusters(m).cluster_locations,...
                       1,'linecolor','k'); hold off; 
                end
            else % mixed/rm-ANOVAs
                FieldsTemp = fieldnames(PermResults.Statistics);
                for t = 1:numel(FieldsTemp)
                    for m=1:length(PermResults.Clusters.(FieldsTemp{t}))
                        hold on; contour(PermResults.Clusters(m).(FieldsTemp{t}).cluster_locations,...
                           1,'linecolor','k'); hold off; 
                    end
                end
            end
        end
    end

    % Colorbar
    hcb1=colorbar;set(get(hcb1,'XLabel'),'String',Units);
    CoLMap=flipud(hot(50));

    % Statistics
    if isfield(PermResults.Statistics,'Obs') % t-tests + one-way ANOVA
        PlotStats = PermResults.Statistics.Obs;
        PlotStats(PermResults.Statistics.P_Values>AlphaThresh)=0;
        subplot(2,NPlots,[NPlots+1,NPlots*2]);
        h=imagesc(PlotStats);hcb2=colorbar; colormap(gca,CoLMap);
        set(gca,'ydir','norm');set(h,'alphadata',PlotStats~=0);xlabel('Frequencies'), ylabel('Channels')
        title(sprintf('%d sig. clusters of electrodes',length(PermResults.Clusters)))
        set(gca,'XTick',XTicks,'XTickLabel',round(SpectFreqs(5:5:end)));
        set(gca,'YTick',YTicks,'YTickLabel',ChannelsLabels(1:2:end));
        set(get(hcb2,'XLabel'),'String','t-values');
    else % mixed/rm-ANOVAs
        FieldsTemp = fieldnames(PermResults.Statistics);
        for t = 1:numel(FieldsTemp)
            PlotStats = PermResults.Statistics.(FieldsTemp{t}).Obs;
            PlotStats(PermResults.Statistics.(FieldsTemp{t}).P_Values>AlphaThresh)=0;
            subplot(2,NPlots,NPlots+t);
            h=imagesc(PlotStats);hcb2=colorbar; colormap(gca,CoLMap);
            set(gca,'ydir','norm');set(h,'alphadata',PlotStats~=0);xlabel('Frequencies'), ylabel('Channels')
            title(sprintf('%d sig. clusters of electrodes: %s',length(PermResults.Clusters.(FieldsTemp{t})),TypeStats{t}))
            set(gca,'XTick',XTicks,'XTickLabel',round(SpectFreqs(5:5:end)));
            set(gca,'YTick',YTicks,'YTickLabel',ChannelsLabels(1:2:end));
            set(get(hcb2,'XLabel'),'String','t-values');
        end
    end

    % Saves the averaged topoplots in export folder
    SaveFigures(gcf,[ExportPath sprintf('FreqByChan_%s',STUDY.name)],'w','bmp');  
    
elseif strcmpi(PlotType,'bands')
    %% STEP 2 : FREQUENCY BANDS (average of frequency bins)
    
    % Initialize variables
    FreqRanges = str2num(FreqData{2});
    FreqName = FreqData{1};
%     PlotData = cellfun(@(x) squeeze(mean(x,2)),SpectData,'UniformOutput',0); 
    PlotData = SpectData;
    MinColor = min(min(cellfun(@(x) min(min(x)),PlotData)));
    MaxColor = max(min(cellfun(@(x) max(max(x)),PlotData)));
    NPlots = size(PlotData,1)*size(PlotData,2);
    SigClust=[];
    ColorString = repmat({'b','r','y','m','c','g'},[1,20]);
    ColorMap = hsv(10);
    TypeStats = {['Main-' STUDY.design.variable(2).label],...
        ['Main-' STUDY.design.variable(1).label],'Interaction'};
    
    % PARENT FIGURE
    figure;
    annotation('textbox',[.4 .5 .5 .5],'String',sprintf('Frequency band: %s [%d - %dHz] ',...
        FreqName,FreqRanges),'FitBoxToText','on','FontSize',20,'EdgeColor','w');
    
    % Group/Condition-specific data
    Pos = 0;
    for k=1:size(PlotData,1) % Groups
        for p = 1:size(PlotData,2) % Conditions
            Pos = Pos + 1;
            subplot(2,NPlots+1,Pos);
            title([strrep(STUDY.group{k},'_','-') '-' STUDY.condition{p}]);
            topoplotIndie(PlotData{k,p},ChanLocs,'sigelect',[]);
            caxis([MinColor MaxColor]); colormap(gca,'jet(50)');
            hcb1=colorbar;set(get(hcb1,'XLabel'),'String',Units);
        end
    end
    
    % Significant clusters
    if isfield(PermResults.TFCE,'Obs') % t-tests + one-way ANOVA
        subplot(2,NPlots+1,Pos+1);
        if ~isempty(PermResults.Cluster_Results)
            SigClust = PermResults.Cluster_Results;
            title(sprintf('%d sig. clusters of electrodes',length(SigClust)));
            ClustRange = {};
            for m=1:length(SigClust)
                TempRange = str2num(strrep(SigClust(m).sample_range,' -',''));
                ClustRange(m,1) = {TempRange(1):TempRange(end)};
                if m==length(SigClust)
                    topoplotIndie(zeros(size(PlotData{1})),ChanLocs,'sigelect',...
                        ClustRange,'colsigelect',ColorString(1:length(SigClust)),'shading','interp');
                end
            end
        else 
            title('0 sig. clusters of electrodes')
            topoplotIndie(zeros(size(PlotData{1})),ChanLocs,'sigelect',[],'shading','interp');
        end

        % Since we duplicate the channel dimension, we only keep the dimension
        % including p-values which are not all identical
        PlotStats = squeeze(PermResults.TFCE.P_Values(1,:)); 

        % Statistics        
        subplot(2,NPlots+1,[NPlots+2,(NPlots+1)*2]);
        plot(PlotStats,'-.o','Color','k');xlabel('Channels'), ylabel('P-values');
        hold on
        ThreshPlot = NaN(size(PlotStats));
        for m=1:length(SigClust)
            ThreshPlot(ClustRange{m}) = PlotStats(ClustRange{m});
            plot(ThreshPlot,'-.o','Color','r');
            ThreshPlot = NaN(size(PlotStats));
        end
        hold off; axis tight; title('Permutation-corrected P-Values')
        set(gca,'XTick',1:length(ChannelsLabels),'XTickLabel',ChannelsLabels);
        line([0 size(PlotStats,2)], [AlphaThresh AlphaThresh],'Color','k','LineStyle','--');
        
    else % mixed/rm-ANOVAs
        FieldsTemp = fieldnames(PermResults.Cluster_Results);
        for t = 1:numel(FieldsTemp)
            subplot(2,NPlots+1,Pos+t);
            if ~isempty(PermResults.Cluster_Results.(FieldsTemp{t}))
                SigClust = PermResults.Cluster_Results.(FieldsTemp{t});
                title(sprintf('%d sig. clusters of electrodes: %s',length(SigClust),TypeStats{t}));
                ClustRange = {};
                for m=1:length(SigClust)
                    TempRange = str2num(strrep(SigClust(m).sample_range,' -',''));
                    ClustRange(m,1) = {TempRange(1):TempRange(end)};
                    if m==length(SigClust)
                        topoplotIndie(zeros(size(PlotData{1})),ChanLocs,'sigelect',...
                            ClustRange,'colsigelect',ColorString(1:length(SigClust)),'shading','interp');
                    end
                end
            else 
                title(sprintf('0 sig. clusters of electrodes: %s',TypeStats{t}))
                topoplotIndie(zeros(size(PlotData{1})),ChanLocs,'sigelect',[],'shading','interp');
            end
        end
        
        subplot(2,NPlots+1,[NPlots+4,(NPlots+1)*2]);
        for t = 1:numel(FieldsTemp)
            % Since we duplicate the channel dimension, we only keep the dimension
            % including p-values which are not all identical
            PlotStats = squeeze(PermResults.TFCE.(FieldsTemp{t}).P_Values(1,:)); 

            % Statistics        
            Fig(t) = plot(PlotStats,'-.o','Color',ColorMap(t,:));xlabel('Channels'), ylabel('P-values');
            hold on
        end
        axis tight; title('Permutation-corrected P-Values')
        set(gca,'XTick',1:2:length(ChannelsLabels),'XTickLabel',ChannelsLabels(1:2:end));
        line([0 size(PlotStats,2)], [AlphaThresh AlphaThresh],'Color','k','LineStyle','--');
        hold off; legend(Fig,TypeStats,'location','best');legend('boxoff');
    end

    % Saves the averaged topoplots in export folder
    SaveFigures(gcf,[ExportPath sprintf('Topoplots_%s_%s',STUDY.name,FreqName)],'w','bmp');
    
    % ------------------------------------------------------------------- %
    % STEP 2.2 : For each sig. result, extract data of the electrode
    
    if isfield(PermResults.TFCE,'Obs') % t-tests + one-way ANOVA
        if ~isempty(PermResults.Cluster_Results)

            for t = 1:size({PermResults.Cluster_Results.cluster_locations},2)
                % BOXPLOTS for each channels
                Temp{t} = PermResults.Cluster_Results(t).cluster_locations;
                SigChans{t} = find(Temp{t}(1,:)==1);
            end
            AllSigChans = sort([SigChans{:}]);

            % Min/Max values for graph
    %         Min = cellfun(@(x) min(min(x(AllSigChans,:))),SpectData,'UniformOutput',0);
    %         Max = cellfun(@(x) max(max(x(AllSigChans,:))),SpectData,'UniformOutput',0);

            % NEW GRAMM TOOLBOX (GGPLOT2-LIKE) FOR BOXPLOTS
            clear Graph;
            Columns = 1; Lines = 1; PlotNum = 1;
            for m=1:length(AllSigChans)

                % Initialize variables
                PlotData = cellfun(@(x) x(AllSigChans(m),:)',SpectData,'UniformOutput',0); 
%                 x1 = SpectData{1}(AllSigChans(m),:)';
%                 x2 = SpectData{2}(AllSigChans(m),:)';
                GrpNum = []; % cellfun(@(x) length(x),PlotData)'; 
                GrpLab = {};
                AllPlotData = [];
                for p=1:length(STUDY.group)
                    GrpLab = [GrpLab; repmat(STUDY.group(p),[length(PlotData{p}),1])];
                    AllPlotData = [AllPlotData; PlotData{p}];
                    GrpNum = [GrpNum; repmat(p-1,[length(PlotData{p}),1])];
                end
    %             grp = [zeros(1,length(x1)),ones(1,length(x2))];

    %             % Setup the graph properties/data
    %             Graph(Lines,Columns)=gramm('x',[repmat(STUDY.group(1),...
    %                 [length(x1),1]);repmat(STUDY.group(2),[length(x2),1])],...
    %             'y',[x1;x2],'color',GrpNum);
                Graph(Lines,Columns)=gramm('x',GrpLab,'y',AllPlotData,'color',GrpNum);
                Graph(Lines,Columns).stat_boxplot('width',0.1,'dodge',0.1);
                Graph(Lines,Columns).stat_violin('fill','transparent');
                Graph(Lines,Columns).set_title(ChannelsLabels{AllSigChans(m)});     

                % Removing labels
                Graph(Lines,Columns).set_names('x','','y','','color',''); 
                if Columns == 1
                    Graph(Lines,Columns).set_names('x','','y',Units,'color','');
    %             else
    %                 Graph(Lines,Columns).set_names('x','','y','','color','');
                end
                Graph(Lines,Columns).no_legend();
    %             Graph(Lines,Columns).axe_property('YLim',[min([Min{:}]) max([Max{:}])]);
                % Graph(Lines,Columns).geom_point();

                % Create the structure of subplots
                if  mod(Columns,6) == 0              
                    Lines = Lines + 1;
                    Columns = 1;
                else
                    Columns = Columns + 1 ;
                end

                % If number of graphs > 30, plot separately
                if m == 30 || m == length(AllSigChans) 

                    % Plotting
                    figure('units','normalized','outerposition',[0 0 1 1]);
    %                 title(sprintf('Frequency band: %s [%d - %dHz]',...
    %                         FreqData{k,1},FreqRanges(k,1),FreqRanges(k,2)));
                    Graph.draw();

                    % Aestethics modifications
                    Columns = 1; Lines = 1;   
                    for r=(PlotNum-1)*30+1:m

                        set([Graph(Lines,Columns).results.stat_boxplot.outliers_handle],'visible','off')
                        % set([Graph(Lines,Columns).results.geom_point_handle],'MarkerEdgeColor','k','MarkerSize',2)

                        % Create the structure of subplots
                        if mod(Columns,6) == 0               
                            Lines = Lines + 1;
                            Columns = 1;
                        else
                            Columns = Columns + 1 ;
                        end
                    end

                    % Saves the boxplot in export folder
                    SaveFigures(gcf,[ExportPath sprintf('Boxplot_%s_%s_%d',FreqName,...
                        STUDY.name, PlotNum)],'w','bmp');  

                    % New graph
                    clear Graph;
                    close all
                    Columns = 1; Lines = 1;   
                    PlotNum = PlotNum + 1;
                end
            end
        end
    else % mixed/rm-ANOVAs
        FieldsTemp = fieldnames(PermResults.Cluster_Results);
        AllSigChans = [];
        for t = 1:numel(FieldsTemp)
            if ~isempty(PermResults.Cluster_Results.(FieldsTemp{t}))
                ClustRes = PermResults.Cluster_Results.(FieldsTemp{t});
                for g = 1:size({ClustRes.cluster_locations},2)
                    % BOXPLOTS for each channels
                    Temp{g} = ClustRes(g).cluster_locations;
                    SigChans{g} = find(Temp{g}(1,:)==1);
                end
                AllSigChans = sort([SigChans{:}]);
            end
        end
        if ~isempty(AllSigChans)
            
            % NEW GRAMM TOOLBOX (GGPLOT2-LIKE) FOR BOXPLOTS
            clear Graph;
            Columns = 1; Lines = 1; PlotNum = 1;
            for m=1:length(AllSigChans)

                % Initialize variables
                PlotData = cellfun(@(x) x(AllSigChans(m),:)',SpectData,'UniformOutput',0); 
                GrpNum = []; 
                GrpLab = {};
                AllPlotData = [];
                for p=1:length(STUDY.group)
                    for f=1:length(STUDY.condition)
                        GrpLab = [GrpLab; repmat(STUDY.group(p),[length(PlotData{p,f}),1])];
                        AllPlotData = [AllPlotData; PlotData{p,f}];
                        GrpNum = [GrpNum; repmat(p-1,[length(PlotData{p}),1])];
                    end
                end
                Graph(Lines,Columns)=gramm('x',GrpLab,'y',AllPlotData,'color',GrpNum);
                Graph(Lines,Columns).stat_boxplot('width',0.1,'dodge',0.1);
                Graph(Lines,Columns).stat_violin('fill','transparent');
                Graph(Lines,Columns).set_title(ChannelsLabels{AllSigChans(m)});     

                % Removing labels
                Graph(Lines,Columns).set_names('x','','y','','color',''); 
                if Columns == 1
                    Graph(Lines,Columns).set_names('x','','y',Units,'color','');
                end
                Graph(Lines,Columns).no_legend();

                % Create the structure of subplots
                if  mod(Columns,6) == 0              
                    Lines = Lines + 1;
                    Columns = 1;
                else
                    Columns = Columns + 1 ;
                end

                % If number of graphs > 30, plot separately
                if m == 30 || m == length(AllSigChans) 

                    % Plotting
                    figure('units','normalized','outerposition',[0 0 1 1]);
                    Graph.draw();

                    % Aestethics modifications
                    Columns = 1; Lines = 1;   
                    for r=(PlotNum-1)*30+1:m

                        set([Graph(Lines,Columns).results.stat_boxplot.outliers_handle],'visible','off')
                        % Create the structure of subplots
                        if mod(Columns,6) == 0               
                            Lines = Lines + 1;
                            Columns = 1;
                        else
                            Columns = Columns + 1 ;
                        end
                    end

                    % Saves the boxplot in export folder
                    SaveFigures(gcf,[ExportPath sprintf('Boxplot_%s_%s_%d',FreqName,...
                        STUDY.name, PlotNum)],'w','bmp');  

                    % New graph
                    clear Graph;
                    close all
                    Columns = 1; Lines = 1;   
                    PlotNum = PlotNum + 1;
                end
            end
        end
    end
    
elseif strcmpi(PlotType,'avgfreqs')
   
    % Initialize variables
    if length(STUDY.design.variable) > 1
        TypeStats = {['Main-' STUDY.design.variable(2).label],...
            ['Main-' STUDY.design.variable(1).label],'Interaction'};
    end
    ColorMap = hsv(10);
    
    % PARENT FIGURE
    figure;
    
    % Group/Condition-specific data
    subplot(1,3,[1,2]); LegendLabel = [];
    set(gca,'XTick',XTicks,'XTickLabel',round(SpectFreqs(5:5:end)));
    Pos = 1;
    for k=1:length(STUDY.group)
        for f=1:length(STUDY.condition)
            Serror=std(SpectData{k,f},0,2) / sqrt( length(SpectData{k,f}));
            YVal=mean(SpectData{k,f},2); 
            XVal=1:size(SpectData{k,f},1);
            shadedErrorBar(XVal,YVal,Serror,'lineprops',{'-bo','MarkerFaceColor',ColorMap(Pos,:),...
                'color',ColorMap(Pos,:)},'patchSaturation',0.2,'transparent',1);
            hold on
            
            LegendLabel = [LegendLabel;{[strrep(STUDY.group{k},'_','-') '-' strrep(STUDY.condition{f},'_','-')]}];
            Pos = Pos + 1;
        end
    end
    axis tight; xlabel('Frequencies'); ylabel(Units)
    
    if isfield(PermResults,'pinter') % mixed/rm-ANOVAs
%         FieldsTemp = fieldnames(PermResults.Cluster_Results);
%         for t = 1:numel(FieldsTemp)
%             % Finding the significant clusters
%             SigClust={};
%             if ~isempty(PermResults.Cluster_Results.(FieldsTemp{t}))
%                 for m=1:length(PermResults.Cluster_Results.(FieldsTemp{t}))
%                     TempRange = str2num(strrep(PermResults.(FieldsTemp{t}).Cluster_Results(m).sample_range,' -',''));
%                     SigClust(m,1) = {TempRange(1):TempRange(end)};
%                 end
%             end
% 
%             % Adding vertical patches (for each sig. cluster) 
%             if ~isempty(PermResults.Cluster_Results.(FieldsTemp{t}))
%                 YLimits = get(gca,'YLim');
%                 for m=1:length(PermResults.Cluster_Results.(FieldsTemp{t}))
%                     XCoord = [SigClust{m}(1) SigClust{m}(end) SigClust{m}(end) SigClust{m}(1)];
%                     YCoord = [YLimits(1)  YLimits(1) YLimits(2) YLimits(2)];
%                     patch(XCoord,YCoord, [([17 17 17])/255], 'FaceAlpha', 0.2); % [([17 17 17])/255] is gray color
%                 end
%             end
%         end
    else % t-tests + one-way ANOVA
        % Finding the significant clusters
        SigClust={};
        Idx = bwlabel(PermResults.mask{:});
        for m = 1:length(unique(Idx))-1 % -1 is for 0s
            TempRange = find(Idx == m);
            SigClust(m,1) = {TempRange};
        end

        % Adding vertical patches (for each sig. cluster) 
        if ~isempty(SigClust)
            YLimits = get(gca,'YLim');
            for m=1:length(SigClust)
                XCoord = [SigClust{m}(1) SigClust{m}(end) SigClust{m}(end) SigClust{m}(1)];
                YCoord = [YLimits(1)  YLimits(1) YLimits(2) YLimits(2)];
                patch(XCoord,YCoord, [([17 17 17])/255], 'FaceAlpha', 0.2); % [([17 17 17])/255] is gray color
            end
        end
    end
    % Legend/Title
    hold off; legend(LegendLabel,'location','best');legend('boxoff');title('Average over all electrodes');

    if isfield(PermResults.TFCE,'pinter') % mixed/rm-ANOVAs
%         FieldsTemp = fieldnames(PermResults.Cluster_Results);
%         for t = 1:numel(FieldsTemp)
%             % Since we duplicate the channel dimension, we only keep the dimension
%             % including p-values which are not all identical
%             PlotStats = squeeze(PermResults.TFCE.(FieldsTemp{t}).P_Values(1,:)); 
% 
%             % Stats results
%             subplot(133); 
%             title(sprintf('%d sig. clusters of frequencies',length(PermResults.Cluster_Results.(FieldsTemp{t}))))
%             Fig(t) = plot(PlotStats,'-.o','Color',ColorMap(end-t+1,:));xlabel('Frequencies'), ylabel('P-values');
%             hold on
%             ThreshPlot = NaN(size(PlotStats));
%             for m=1:length(SigClust)
%                 ThreshPlot(SigClust{m}) = PlotStats(SigClust{m});
%                 plot(ThreshPlot,'-.o'); % ,'Color',ColorString{m}
%                 ThreshPlot = NaN(size(PlotStats));
%             end
%         end
%         hold off; axis tight; title('Permutation-corrected P-Values');
%         line([0 size(PlotStats,2)], [AlphaThresh AlphaThresh],'Color','k','LineStyle','--');
%         set(gca,'XTick',XTicks,'XTickLabel',round(SpectFreqs(5:5:end))); 
%         legend(Fig,TypeStats','location','best');legend('boxoff');
     else % t-tests + one-way ANOVA
         
        % Since we duplicate the channel dimension, we only keep the dimension
        % including p-values which are not all identical
        PlotStats = PermResults.stats{:}; 

        % Stats results
        subplot(133); 
        title(sprintf('%d significant clusters of frequencies',length(SigClust)))
        plot(PlotStats,'-.o','Color','k');xlabel('Frequencies'), ylabel('T-values');
        hold on
        ThreshPlot = NaN(size(PlotStats));
        for m=1:length(SigClust)
            ThreshPlot(SigClust{m}) = PlotStats(SigClust{m});
            Fig=plot(ThreshPlot,'-.o','Color','r'); 
            ThreshPlot = NaN(size(PlotStats));
        end
        hold off; axis tight; title('Permutation-corrected P-Values');
        line([0 size(PlotStats,2)], [AlphaThresh AlphaThresh],'Color','k','LineStyle','--');
        set(gca,'XTick',XTicks,'XTickLabel',round(SpectFreqs(5:5:end))); 
    end

    % Saves the averaged topoplots in All folder
    SaveFigures(gcf,[ExportPath sprintf('WholeSpect_%s',STUDY.name)],'w','bmp');
end

end