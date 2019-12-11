function Log = ICClustLocalisation(STUDY,ALLEEG,varargin)
% This function is used to cluster independent components across all
% subjects. Data are issued from EEGLab STUDY structure and the Measure
% Projection Toolbox (MST):
% https://sccn.ucsd.edu/wiki/MPT

% Usage:
%    >> ICClustLocalisation(STUDY,varargin);
%
% Inputs (mandatory):
%   STUDY            = EEGlab STUDY structure

% Inputs (Optional):
%   ROI                           = Logical indicating whether to XXX
%   headGridSpacing               =
%   NPermut                       =
%   stdOfDipoleGaussian           =
%   numOfStdsToTruncateGaus       =
%   normalizeInBrainDipoleDensity =

% Outputs:
%   Log              = Structure containing parameters information for each
%                      version that the user ran.

% Author: Corentin Wicht, LCNS, 2019
% - corentin.wicht@unifr.ch
% - https://github.com/CorentinWicht

%% Parameters
ROI = 0;
headGridSpacing = 8; % Default from headGrid.m
NPermut = 5000;
stdOfDipoleGaussian = 12; % Default from meanProjection.m
numOfStdsToTruncateGaus = 3; % Default from meanProjection.m
normalizeInBrainDipoleDensity = 'on'; % Default from meanProjection.m
cutoffRatio = [1 0.05]; % Default from dipole.m
% [1 0.05]) indicates that we want all the ICs that at least has a 0.05 
% chance of being in the domain. You may want to use 0.1 or even 0.5 to get fewer ICs.
FDRCorrect = false; % Default is not to perform FDR correction 
Pval = 0.05; % Default for significance
SpectMaxCorrel = 0.9; % Default for ???? 
FreqRange = [1 40]; % Default for the frequency bands
cutoffDensity = 0.05; % Anatomical areas with less than this dipole denisty ratio will not be reported
ExportPath = [pwd '\TempICClust\']; % THIS IS TEMPORARY!!!!!!!!!!!!!
mkdir ([ExportPath '\TempICClust\']);
ExcelDirectory = [pwd 'TempICClust\']; % THIS IS TEMPORARY!!!!!!!!!!!!!
Measure =  'spec';

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
            case 'roi'
                ROI  = Value;
            case 'headgridspacing'
                headGridSpacing  = Value;
            case 'npermut'
                NPermut  = Value;
            case 'stdofdipolegaussian'
                stdOfDipoleGaussian  = Value;    
            case 'numofstdstotruncategaus'
                numOfStdsToTruncateGaus  = Value;
            case 'normalizeinbraindipoledensity'
                normalizeInBrainDipoleDensity  = Value;
            case 'cutoffratio'
                cutoffRatio  = Value;
            case 'fdrcorrect'
                FDRCorrect  = Value; 
            case 'pval'
                Pval  = Value;
            case 'spectmaxcorrel'
                SpectMaxCorrel  = Value;
            case 'freqrange'
                FreqRange  = Value; 
            case 'cutoffdensity'
                cutoffDensity = Value;
            case 'exportpath'
                ExportPath = Value;
            case 'exceldirectory'
                ExcelDirectory  = Value;
            case 'allparameters'
                AllParameters = Value;
                
                ROI = strcmpi(AllParameters{1},'Y');
                headGridSpacing = str2double(AllParameters{2}); 
                stdOfDipoleGaussian = str2double(AllParameters{3}); 
                numOfStdsToTruncateGaus = str2double(AllParameters{4}); 
                if strcmpi(AllParameters{5},'Y')
                    normalizeInBrainDipoleDensity = 'on'; 
                else
                    normalizeInBrainDipoleDensity = 'off';
                end
                cutoffRatio = str2num(AllParameters{6});
                Pval = str2double(AllParameters{7}); 
                FDRCorrect = strcmpi(AllParameters{8},'Y');
                SpectMaxCorrel = str2double(AllParameters{9});  
                cutoffDensity = str2double(AllParameters{10});
                Measure =  AllParameters{11};
                FreqRange = str2num(AllParameters{12}); 
            otherwise
                display (['Unknown parameter setting: ' Param])
        end
    end
end

%% TEMPORARY FOR TESTS
% Adding path to dependencies
addpath([pwd '\Functions\eeglab14_1_2b']);
addpath([pwd '\Functions\']);
addpath(genpath([pwd '\Functions\Dependencies']));
addpath([pwd '\Functions\EEGInterp']); 
eeglab; close gcf;
STUDYName = 'STUDY.study'; % 'STUDY_Clustered.study'
STUDYPath = 'E:\Dropbox\DOC1\Frequency_Source_Analysis\Sommeil_Joelle\FINALSoftware\DATASETS\FarfallaOH(1BS0WS)\OUTPUTS\STUDY\3_10_2019-8_45\';
% Load STUDY
pop_editoptions( 'option_storedisk', 1);
[STUDY ALLEEG] = pop_loadstudy('filename', STUDYName,'filepath', STUDYPath);

%% STEP 0 : ICA rejection for clustering

% Prompt
ConfirmReject = questdlg(['Do you need to reject artifactual independent components?' newline ...
    'For the IC clustering, you should only keep 5-15 "true" brain IC for each subject.'], ...
     'Reject more components ?', 'No','Yes','Yes');
if strcmpi(ConfirmReject,'Yes')
    ConfirmAuto = questdlg(['Would you like to make the rejection automatic?' newline ...
        'If yes, artifactual ICs and those with dipole RV > 15% and/or outside the brain will be rejected.'], ...
         'Automatic rejection', 'No','Yes','No'); 
end

if strcmpi(ConfirmReject,'Yes')    
    for k = 1:length(ALLEEG)
        % Initialize the analysis for each data
        Restart=0;

        % Loading datasets one by one
        EEG = pop_loadset('filename',ALLEEG(k).filename,'filepath',ALLEEG(k).filepath);
        CURRENTSET = k;

        % ICLabel plugin
        EEG=iclabel(EEG);
        
        % Retrieve results (pre-select all non-brain components)
        Rmdip = [];
        CompsToRej = [];
        AllComps = 1:size(EEG.icaact,1);
        for l=1:size(EEG.icaact,1)
            [~,CompType] = max(EEG.etc.ic_classification.ICLabel.classifications(l,:));
            if CompType ~= 1 ... % Brain component, see EEG.etc.ic_classification.ICLabel.classes
                    ||  EEG.dipfit.model(l).rv > 0.15 % Adding to the pre-selected components each one whose dipole fitting residual variance > 15%
                CompsToRej =[CompsToRej l];
            end

            % Identify dipoles outside brain (taken from pop_multifit())
            if ~isempty(EEG.dipfit.model(l).posxyz)
                if any(sqrt(sum(EEG.dipfit.model(l).posxyz.^2,2)) > 85)
                    Rmdip = [Rmdip l];
                    CompsToRej =[CompsToRej l];
                end
            end
        end
        if ~isempty(Rmdip) && strcmpi(ConfirmAuto,'No')
                fprintf('Out of cortex dipoles identified (usually artifacts):\n  %s \n ',num2str(Rmdip));
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
            
            % Find list of components to keep
            CompsToKeep = find(~ismember(AllComps,CompsToRej));
            
            % Remove the duplicates
            CompsToRej = unique(CompsToRej);

            % Only run this part if semi-automatic rejection process
            if strcmpi(ConfirmAuto,'No')
                % Visualize the results
                pop_viewprops( EEG, 0, 1:size(EEG.icaact,1), {'freqrange', [1 60]}, {}, 1, 'ICLabel' )

                % Move figures all over the screen
                ScreenPos={'northwest','northeast','southeast','southwest'};
                for k=1:length(findobj('type','figure'))        
                    movegui(figure(k),ScreenPos{k})
                end

                % Matrix to integrate in the following uitable
                Response = repmat({false},[size(EEG.icaact,1) 1]);
                Response(unique(CompsToRej))={true};
                to_display = [num2cell(1:size(EEG.icaact,1))', Response];
                CompList = {};
                CompsToRej = [];

                % Select components to reject
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

                    % Find list of components to keep
                    CompsToKeep = find(~ismember(AllComps,CompsToRej));
                end

                % Close all figures
                close all
            
            % If automatic IC removal
            else
                
                % Information
                fprintf('File %s; The following components are rejected: %s\n',...
                    ALLEEG(k).setname,num2str(CompsToRej))
                
                % Ends the loop 
                Restart = 1;

            end
        end

        % Update the study information
        [STUDY ALLEEG] = std_editset( STUDY, ALLEEG, 'commands',...
            {{'index' k 'comps' CompsToKeep }},'updatedat','off','rmclust','off' );
    end

    %% RECOMPUTING COMPONENTS MEASURES
    % Precompute Power Spectra
    [STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, 'components','allcomps','off','recompute','on',...
        'scalp','on','spec','on','specparams',{'specmode' STUDY.SpecMode 'logtrials' 'on'}); 

    % Create preclustering array
    [STUDY,ALLEEG] = std_preclust(STUDY,ALLEEG,[],{ 'spec'  'npca' 10 'norm' 1 ...
        'weight' 1},{ 'scalp' 'npca' 10 'norm' 1 'weight' 1 'abso' 1 },...
        { 'dipoles' 'norm' 1 'weight' 10 });    
    
    % Check for inconsistencies
    [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);
    
    % Saving the STUDY under a new name
    [STUDY EEG] = pop_savestudy(STUDY, EEG, 'filename',[STUDY.name '_Clustered.study'],...
    'filepath',STUDY.filepath);
    
    % No idea why an empty figure is being generated
    close gcf
end

% Restricting frequency range
[STUDY,~,SpectFreqs] = std_readspec(STUDY, ALLEEG,'clusters',1,'freqrange', FreqRange); 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                IC CLUSTERING (MEASURE PROJECTION TOOLBOX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: PROJECT SETUP
% Code taken from : https://sccn.ucsd.edu/wiki/MPT#Tutorial
 
% Domains will be based on anatomical ROI's "super voxels"
% WHAT DOES THIS MEAN? CAN WE SELECT ROI's OR NOT ?? 
if ROI
   roiProjection = pr.regionOfInterestProjection(STUDY.measureProjection.(Measure).object, ...
       STUDY.measureProjection.(Measure).object.getPairwiseFishersZSimilarity, pr.headGrid); 
   roiProjection.makeReport('ROIBasedReport.txt',STUDY.filepath);
end

% read the data from STUDY
% Corrected a bug in the function below at line 67
% There was an inversion between groups/condition position
STUDY.measureProjection.(Measure).object = pr.dipoleAndMeasureOfStudySpec(STUDY, ALLEEG); 

% define HeadGRID based on GUI options 
% Should investigate the arguments !!! 
STUDY.measureProjection.(Measure).headGrid = pr.headGrid(headGridSpacing);
HeadGrid = STUDY.measureProjection.(Measure).headGrid;

% Plotting the headgrid
% figure;STUDY.measureProjection.(Measure).headGrid.plot;

% Creating a source-space neighbouring matrix for TFCE statistics
% Position = HeadGrid.getPosition(HeadGrid.insideBrainCube);
% for m=1:size(Position,1)
%     NodesLoc(m).X = Position(m,1);
%     NodesLoc(m).Y = Position(m,2);
%     NodesLoc(m).Z = Position(m,3);
% end
% ChN = ept_ChN2(NodesLoc, 1);

% Modifying options according to user input
STUDY.measureProjection.option.properties(14).Value = FDRCorrect;
STUDY.measureProjection.option.specFdrCorrection = FDRCorrect;
STUDY.measureProjection.option.standardDeviationOfEstimatedDipoleLocation = stdOfDipoleGaussian;
STUDY.measureProjection.option.numberOfStandardDeviationsToTruncatedGaussaian = numOfStdsToTruncateGaus;
STUDY.measureProjection.option.normalizeInBrainDipoleDenisty = strcmpi(normalizeInBrainDipoleDensity,'on');
STUDY.measureProjection.option.numberOfPermutations = NPermut;
STUDY.measureProjection.option.headGridSpacing = headGridSpacing;
STUDY.measureProjection.option.specMaxCorrelation = SpectMaxCorrel;
STUDY.measureProjection.option.specSignificance = Pval;

% Perform the actual projection (using permutation statistics)
STUDY.measureProjection.(Measure).projection = pr.meanProjection(STUDY.measureProjection.(Measure).object,...
STUDY.measureProjection.(Measure).object.getPairwiseMutualInformationSimilarity, ...
STUDY.measureProjection.(Measure).headGrid, 'numberOfPermutations', ...
NPermut, 'stdOfDipoleGaussian',stdOfDipoleGaussian,'numberOfStdsToTruncateGaussian',...
numOfStdsToTruncateGaus, 'normalizeInBrainDipoleDenisty', ... 
fastif(normalizeInBrainDipoleDensity,'on', 'off'));

%% Plotting the significant projections

% visualize significant voxels 
figure;STUDY.measureProjection.(Measure).projection.plotVoxel(Pval);
SaveFigures(gcf,[ExportPath sprintf('ProjVoxelSigProjection_%s_%sHz',STUDY.name,strrep(num2str(FreqRange),' ','-'))],'k','bmp');

% visualize significant voxels projected in the cortex template 
STUDY.measureProjection.(Measure).projection.plotCortex(Pval);
saveas(gcf,[ExportPath sprintf('ProjCortexSigProjection_%s_%sHz.fig',STUDY.name,strrep(num2str(FreqRange),' ','-'))])
close gcf

%% STEP 2: DOMAINS CREATION

% get the spec and dipole data (dataAndMeasure object) from the STUDY structure
dipoleAndMeasure = STUDY.measureProjection.(Measure).object; 
 
% find out the significance level to use (e.g. corrected by FDR)
% THERE MAY BE ERRORS HERE IF p-FDR = 0 !!!!! WHY ??? 
if STUDY.measureProjection.option.('specFdrCorrection')
   significanceLevel = fdr(STUDY.measureProjection.(Measure).projection.convergenceSignificance(...
STUDY.measureProjection.(Measure).headGrid.insideBrainCube(:)), STUDY.measureProjection.option.('specSignificance'));
else
   significanceLevel = STUDY.measureProjection.option.('specSignificance');
end
maxDomainExemplarCorrelation = STUDY.measureProjection.option.('specMaxCorrelation');

% the command below creates the domains using parameters significanceLevel and maxDomainExemplarCorrelation:
STUDY.measureProjection.(Measure).projection = STUDY.measureProjection.(Measure).projection.createDomain(...
dipoleAndMeasure, maxDomainExemplarCorrelation, significanceLevel);

%% Plotting domains
% visualize domains
% TO DO: set up a semi-automatic loop here which enables adjusting visual
% parameters of the figure parameters ?!!! 

% DOMAINS
figure;STUDY.measureProjection.(Measure).projection.plotVoxelColoredByDomain;
SaveFigures(gcf,[ExportPath sprintf('MeasVoxelColoredByDomain_%s_%sHz',STUDY.name,strrep(num2str(FreqRange),' ','-'))],'k','bmp');

STUDY.measureProjection.(Measure).projection.plotCortexColoredByDomain('plotOptions',{'minimumDomainLight',0.5,'designer','nima'});
saveas(gcf,[ExportPath sprintf('MeasCortexColoredByDomain_%s_%sHz.fig',STUDY.name,strrep(num2str(FreqRange),' ','-'))])
close gcf

% SIGNIFICANCCE
STUDY.measureProjection.(Measure).projection.plotVoxelColoredBySignificance('colormap','jet');
SaveFigures(gcf,[ExportPath sprintf('MeasVoxelColoredBySignificance_%s_%sHz',STUDY.name,strrep(num2str(FreqRange),' ','-'))],'k','bmp');

% MEASURE
STUDY.measureProjection.(Measure).projection.plotVoxelColoredByMeasure(significanceLevel, dipoleAndMeasure,'colormap','ycrcb');
SaveFigures(gcf,[ExportPath sprintf('MeasVoxelColoredByMeasure_%s_%sHz',STUDY.name,strrep(num2str(FreqRange),' ','-'))],'k','bmp');

%% STEP 3: FINDING ICs ASSOCIATED WITH EACH DOMAIN 
% Prepare the export in excel templates, or in txt files (might be easier)

for k = 1:length(STUDY.measureProjection.(Measure).projection.domain) 
    
    mkdir([ExportPath sprintf('Domain%d',k)])
    domainNumber = k;
    CurrentDomain = STUDY.measureProjection.(Measure).projection.domain(domainNumber); % get the domain in a separate variable
    projection  = STUDY.measureProjection.(Measure).projection;
    [dipoleId sortedDipoleDensity orderOfDipoles dipoleDensity dipoleDensityInRegion] = ...
        dipoleAndMeasure.getDipoleDensityContributionToRegionOfInterest(CurrentDomain.membershipCube, projection, cutoffRatio);
     domainICs.(sprintf('Domain%d',k)) = dipoleAndMeasure.createSubsetForId(dipoleId); 
    % here we create a new variable that contain information only for dipoles associates with domain ICs.
    
    %% EXPORTING SESSIOn PROJECTED MEASURES FOR EACH DOMAIN
    headGrid = STUDY.measureProjection.(Measure).headGrid;
    [linearProjectedMeasure sessionConditionCell groupId uniqeDatasetId dipoleDensity] = ...
        dipoleAndMeasure.getMeanProjectedMeasureForEachSession(headGrid,...
        CurrentDomain.membershipCube, projection.projectionParameter);
    
    % Anatomical information 
    [brodmannAreaNumber brodmannAreaDipoleDensityRatio additionalDescriptionForBrodmannArea] = ...
        getBrodmannArea(CurrentDomain, cutoffDensity,  'printDomainLabel', 0);
    [regionLabelOutput domainDipoleMassInAnatomicalRegionOutput] = ...
        getAnatomicalInformation(CurrentDomain, cutoffDensity,  'printDomainLabel', 0);
    for t=1:length(domainDipoleMassInAnatomicalRegionOutput)
        regionLabelOutputFull{t} = [num2str(round(domainDipoleMassInAnatomicalRegionOutput(t)*100,2)) '%-' regionLabelOutput{t}];
    end
    
    % Create subject space with domain 4 as a ROI
    % NOW IDEA WHAT THIS IS !!! 
%     subjectSpace = pr.subjectSpace(STUDY.measureProjection.(Measure).object, headGrid,...
%         STUDY.measureProjection.(Measure).projection.projectionParameter, domain.membershipCube);
%     subjectSpace.plot; 

    % Plotting the domain in voxels
%     figure;CurrentDomain.plotVoxel
%     CurrentDomain.plotCortex('plotOptions',CurrentDomain.color);

    %% STATISTICS

    % Build the design
    % BASED ON THAT WILL DEFINE WHAT TO USE BELOW TO RETRIEVE DATA !! 
    if length(STUDY.group) > 1
        Design.Between = [{'B'};'group';STUDY.group'];
        StatsIdx = 'i';
    end
    if length(STUDY.condition) > 1
        Design.Within = [{'W'};'condition';STUDY.condition'];
        StatsIdx = 'd';
    end

    % Retrieving the data for each condition
    % This comes from: getMeasureDifferenceAcrossGroups
    % go through all sessions and project session dipoles to the given position(s)
    [linearProjectedMeasure,datasetGroupNumber] = ...
        dipoleAndMeasure.getProjectedMeasureForEachSession(CurrentDomain.headGrid,...
        CurrentDomain.membershipCube, CurrentDomain.projectionParameter);            

    % THIS determine whether the analysis will be done on all individual
    % voxels or on the mean of voxels
    AllVoxels = 'each';
    
    if strcmpi(AllVoxels, 'each')
        % add all projected measures from all subjects together, Multiple position for each
        % subject will contribute.
        linearProjectedMeasureReshaped = reshape(linearProjectedMeasure, size(linearProjectedMeasure,1), []);
%         datasetGroupNumber = repmat(datasetGroupNumber, size(linearProjectedMeasure,2),1);
%         datasetGroupNumber =  reshape(datasetGroupNumber, 1, []);
    else % when usePositionProjections = 'mean'
        % average the linearized measure over all the given position, which are in the second
        % dimension.
        linearProjectedMeasureReshaped = squeeze(mean(linearProjectedMeasure,2)); % NEED TO RETRIEVE DATA FOR EACH VOXELS ?? 
    end

%     % separate conditions from linearized form
%     projectedMeasure = dipoleAndMeasure.getSeparatedConditionsForLinearizedMeasure(linearProjectedMeasureForCombinedCondition);
%     projectedMeasure = projectedMeasure{:};
% 
%     % Separate the data in each group and normalizing the data using the 
%     % substracted mean spectrum over all ics and conditions
%     for groupNumber = 1:dipoleAndMeasure.numberOfGroups
%         linearProjectedMeasureGroups{groupNumber} = projectedMeasure(:,datasetGroupNumber==groupNumber)...
%             + dipoleAndMeasure.specMeanOverAllIcAndCondition; % taken from plotMeasureAsArrayInCell
%     end
% 
%     % Building AllData structure
%     % UNSURE IT WILL WORK WITH WITHIN FACTORS ?!
%     SubjFields = {STUDY.datasetinfo.subject};
%     for t = 1:length(STUDY.group)
%         Pos = 1;
%         for p=1:size(linearProjectedMeasureGroups{t},2)
%             for m = 1:length(STUDY.condition)
%                 Fields = [sprintf('Subj%s',SubjFields{p}) '_' STUDY.group{t} '_' STUDY.condition{m}];
%                 AllData.(Fields) = linearProjectedMeasureGroups{t}(:,Pos);
%                 Pos = Pos + 1;
%             end
%         end
%     end

    % Building AllData structure
    SubjFields = {STUDY.datasetinfo.subject};
    for p=1:size(linearProjectedMeasure,3)
        Fields = [sprintf('Subj%s',SubjFields{p}) '_' STUDY.group{datasetGroupNumber(p)} '_' STUDY.condition{1}];
        AllData.(Fields) = squeeze(linearProjectedMeasure(:,:,p));
    end

    % if data is 1D, duplicate the data in a second dimension
    if size(AllData.(Fields),2) == 1 % In case all voxel are averaged
        AllData = structfun(@(x) repmat(x,[1,size(x,1)]),AllData,'UniformOutput',0);
    else % In case voxels are treated individually
        AllData = structfun(@(x) permute(x,[2,1]),AllData,'UniformOutput',0);
    end
    
    % Retrieve the position of cubes inside the domain
    Position = headGrid.getPosition(CurrentDomain.membershipCube);
    for m=1:size(Position,1)
        NodesLocDomain(m).X = Position(m,1);
        NodesLocDomain(m).Y = Position(m,2);
        NodesLocDomain(m).Z = Position(m,3);
    end
    
    % Create a source-space neighbouring matrix of these cubes
    ChN = ept_ChN2(NodesLocDomain,0);
    NodesLocStruct.chanlocs = NodesLocDomain;
    NodesLocStruct.srate = ALLEEG(1).srate;

    % Run the statistics
    [StatsResults,Cluster_Results]=Perm_Ttest(AllData,Design,NodesLocStruct,'N_Permutes',NPermut,...
        'Pval',Pval,'root_folder',pwd,'StatsType',StatsIdx,'ChN',ChN); 

    %% PLOTTING

    % Only plotting the significant clusters
    if ~isempty(Cluster_Results)
        for g = 1:length(Cluster_Results)
            % Finding the significant voxels
            [SigVoxels,SigFreqs] = find(Cluster_Results(g).cluster_locations==1);
            ColorMap = hsv(10+length(Cluster_Results));Pos = 1;figure;subplot(121);            
            for p = 1:length(STUDY.group)
                for f = 1:length(STUDY.condition)
                    linearProjectedMeasureGrp = squeeze(mean(linearProjectedMeasure(:,SigVoxels,datasetGroupNumber==p),2));
                    PlotData{p,f} = linearProjectedMeasureGrp; 
                    
                    % Plotting data curves (Mean +- standard error)
                    hold on
                    Serror=std(PlotData{p,f},0,2) / sqrt( length(PlotData{p,f}));
                    YVal=mean(PlotData{p,f},2); 
                    XVal=1:size(PlotData{p,f},1);
                    shadedErrorBar(XVal,YVal,Serror,'lineprops',{'-bo','MarkerFaceColor',ColorMap(Pos,:),...
                        'color',ColorMap(Pos,:)},'patchSaturation',0.2,'transparent',1);
                    Pos = Pos + 1;
                end
            end
            hold off
            
            % Settings regarding labels, title
            ylabel('10*Log10(\muV^2/Hz)');
            title(sprintf('Mean +- standard error for sig. voxels in domain %d',k),'Color','w')
            set(gca,'xcolor','w','ycolor','w') 
            set(gca,'color','k');
            % Finding the significant clusters
            TempRange = str2num(strrep(Cluster_Results(g).sample_range,' -',''));
            SigClust = TempRange(1):TempRange(end);

            % Adding vertical patches (for each sig. cluster) 
            YLimits = get(gca,'YLim');
            XCoord = [SigClust(1) SigClust(end) SigClust(end) SigClust(1)];
            YCoord = [YLimits(1)  YLimits(1) YLimits(2) YLimits(2)];
            if length(unique(XCoord))>1
                patch(XCoord,YCoord, [1 1 1], 'FaceAlpha', 0.2); 
            else
                xline(unique(XCoord),'Color','w');
            end
            
            % Legend
            hold off;[h, ~, plots] = legend(STUDY.group,'location','best');
            for idx = 1:length(h.String)
                h.String{idx} = ['\color[rgb]{' num2str(plots(idx).Color) '} ' h.String{idx}];
            end
            legend('boxoff');
   
     
            % Brain with significant voxel
            subplot(122)
            plot_dipplot_with_cortex;
            inputColormap = jet();
            statisticalPower = -log(StatsResults.P_Values(SigVoxels,SigFreqs));
            voxelColor = value2color(statisticalPower, inputColormap);
            
            % Atlas
            annotation('textbox',[.67 .5 .5 .5],'String',regionLabelOutputFull,...
                'FitBoxToText','on','FontSize',12,'EdgeColor','w','Color','w');
            
            % Find the significant voxels
            % Taken from pr.meanProjection()
            for l = 1:length(SigVoxels)
                membershipCube = false(size(CurrentDomain.membershipCube));
                TEMP{1} = find(headGrid.xCube == NodesLocDomain(SigVoxels(l)).X);
                TEMP{2} = find(headGrid.yCube == NodesLocDomain(SigVoxels(l)).Y);
                TEMP{3} = find(headGrid.zCube == NodesLocDomain(SigVoxels(l)).Z);
                [C,ia,ib] = intersect(TEMP{1},TEMP{2}, 'stable');
                [C,ia,ib] = intersect(C,TEMP{3}, 'stable');
                TempC(l) = C;
                membershipCube(C) = true;
                pr.plot_head_region(projection.headGrid, membershipCube,...
                'regionColor', voxelColor(l,:), 'showProjectedOnMrs', true, 'projectionAlpha', 0.1);
            end
            membershipCube(TempC) = true;
            
            % Colorbar
            colormap(inputColormap);
            handle = cbar('vert', 1:size(inputColormap,1), [Pval min(min(StatsResults.P_Values(SigVoxels,SigFreqs)))]);
            cbarLabelCharArray = get(handle, 'yticklabel');
            cbarLabel = {};
            for i=1:size(cbarLabelCharArray,1)
                cbarLabel{i} = strtrim(cbarLabelCharArray(i,:));
            end

            % we need fliplr since the order of values in cbar label is actually in decreasing order
            actualPvalues = fliplr(logspace(log10(min(min(StatsResults.P_Values(SigVoxels,SigFreqs)))), log10(Pval), length(cbarLabel)));
            actualPvalues = cellfun(@(x) num2str(round(x,3)),num2cell(actualPvalues),'UniformOutput',0);
            
            set(handle, 'yticklabel',actualPvalues)
            set(handle, 'ycolor', [1 0 0])
            % remove on y axis text and make the colorbar shorter.
            set(handle, 'position', [ 0.9500    0.01    0.0310    0.23])
            % Change background to white
            set(gcf,'color','k');
            
            % Saves the figure in corresponding Domain folder
            saveas(gcf,[ExportPath sprintf('Domain%d',k) '\' sprintf('SourceClust%d_%s',g,STUDY.name)])
            SaveFigures(gcf,[ExportPath sprintf('Domain%d',k) '\' sprintf('SourceClust%d_%s',g,STUDY.name)],'w','bmp');
        end
    end
    
    % Clear data for next domain
    clear NodesLocDomain NodesLocStruct StatsResults Cluster_Results AllData PlotData
        
    %% Creating a .txt file output for each domain
    
    % Initialize variables
    username=getenv('USERNAME');

    % Creating the log file
    fid = fopen([ExportPath sprintf('Domain%d',k) '\' sprintf('ICClustOutput_Domain%d',k) '.txt'],'w');    
    
    % User name
    fprintf(fid,'%s\r\n',['Windows username : ' username]);
    fprintf(fid,'\r\n%s\r\n',['Domain number : ' num2str(k)]);
    fprintf(fid,'\r\n\r\n%s\r\n','------ ANATOMICAL INFORMATION ------');
    fprintf(fid,'%s\r\n','[Dipole mass - anatomical label based on the LPBA40 probabilistic atlas]');
    fprintf(fid,'\r\n%s',regionLabelOutputFull{:});
    SigPosition = headGrid.getPosition(membershipCube);
    
    if ~isempty(Cluster_Results)
        fprintf(fid,'\r\n\r\n%s\r\n','------ SIGNIFICANT RESULTS ------');
        % preparing output for the .txt file
        TValues = StatsResults.Obs(SigVoxels,SigFreqs);
        PValues = StatsResults.P_Values(SigVoxels,SigFreqs);
        for f = 1:size(TValues,1)
           fprintf(fid,'%s\r\n',sprintf('Voxel%d: t = %f, p = %f, Position[XYZ]:%s',...
               f,TValues(f,1),PValues(f,1),num2str(SigPosition(f,:)))); 
        end
    end
    
    % Closing the file
    fclose(fid);
end

% Save outputs for the LOG
Log.Parameters.Measure = Measure;
Log.Parameters.ROI = ROI;
Log.Parameters.headGridSpacing = headGridSpacing;
Log.Parameters.NPermut = NPermut;
Log.Parameters.stdOfDipoleGaussian = stdOfDipoleGaussian;
Log.Parameters.numOfStdsToTruncateGaus = numOfStdsToTruncateGaus;
Log.Parameters.normalizeInBrainDipoleDensity = normalizeInBrainDipoleDensity;
Log.Parameters.cutoffRatio = cutoffRatio;
Log.Parameters.FDRCorrect = FDRCorrect;
Log.Parameters.Pval = Pval;
Log.Parameters.SpectMaxCorrel = SpectMaxCorrel;
Log.Parameters.FreqRange = FreqRange;
Log.Parameters.cutoffDensity = cutoffDensity;
end