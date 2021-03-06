function [Results,Cluster_Results]=Perm_ANOVA(AllData,Design,EEG,varargin)

% This function computes :
% 1) ONE-WAY ANOVA
% 2) MIXED (balanced and unbalanced) ANOVA 
% 3) REPEATED-MEASURES ANOVA (currently 1 WS factor ONLY)
% relying on permutation statistics + TFCE correction for multiple comparisons.
% It uses functions from the FMUT & ept_TFCE-MATLAB toolboxes.
% FMUT : https://github.com/ericcfields/FMUT
% ept_TFCE : https://github.com/Mensen/ept_TFCE-matlab

% Usage:
%    >> [Results,Cluster_Results]=Perm_ANOVA(AllData,Design,EEG,varargin);
%
% Inputs (mandatory):
%   ALLDATA = Structure containing each subject data in separate fields
%   (all conditions, groups, etc together)

%   DESIGN  = Design is a structure containing fields (Between or Within). 
%   Each field is a cell matrix of strings (Nx1) with the:
%   - 1st line containing either a 'B' (Between) or 'W' (Within)
%   - 2nd line containing the name of the factor (e.g. 'Group')
%   - 3rd to last lines containing the name of the levels (e.g. OH, PBO)

% --> CURRENTLY CAN ONLY ACCOMODATE 1-Bs and 1-Ws FACTORS ! 

%   EEG = EEGLab data structure (to retrieve the chanlocs field)

% Inputs (Optional):
%   N_PERMUTES  = [Integer] Number of permutations (def=5000)

%   TFCE  = [Flag] Threshold-free cluster enhancement (def=Yes)

%   PVAL = [Integer] Threshold for significance (def=0.05)

%   E_H = [Integer] TFCE correction parameters (def=[0.666, 1])

%   ROOT FOLDER = [String] Upper folder path used to recompile the c-files 
%   using MEX, in case of errors. Will be bypassed if no argument is passed (def='').

%   CHN = [X-by-X Matrix] Channel/source-space neihbouring matrix generated
%   with the function ept_ChN2.

% Outputs:
%   RESULTS = Data results structure containing the fields:
%   - F_Obs : Observed F-statistics
%   - TFCE_Obs : Observed TFCE-corrected F-statistics
%   - maxTFCE : Maximum TFCE values for each permutation maps 
%   - P_Values : Matrix of P-values
%   - FlagTFCE : Indicate presence/absence of TFCE computation

%   CLUSTER_RESULTS = Results of the clustering procedure 
%   (for details see function ept_calculateClusters)

% Author: Corentin Wicht, LCNS, 2019
% - corentin.wicht@unifr.ch
% - https://github.com/CorentinWicht

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

%% Set Defaults

N_Permutes  = 5000; % default number of permutations
TFCE        = 1; % default threshold-free cluster enhancement
Pval        = 0.05; % default threshold for significance
E_H         = [0.666, 1]; % default parameters of E and H
root_folder = ''; % default upper path containing all files
ChN = []; % default for the channel/source neighbouring matrix
 
% Process Secondary Arguments
if nargin > 4
    for i = 1:2:length(varargin)
        Param = varargin{i};
        Value = varargin{i+1};
        if ~ischar(Param)
          error('Flag arguments must be strings')
        end
        Param = lower(Param);
    
        switch Param
            case 'n_permutes'
                N_Permutes  = Value;
            case 'pval'
                Pval        = Value;
            case 'e_h'
                E_H         = Value;
            case 'root_folder'
                root_folder = Value;
            case 'tfce'
                TFCE = Value;
            case 'chn'
                ChN = Value;
            otherwise
                display (['Unknown parameter setting: ' Param])
        end
    end
end


%% SETTINGS
StringANOVA={'A','B','AB'};
cluster_thresh = zeros(1,3);
Contrasts = {3,4,[3,4]}; % Within, Between, Interaction

% Retrieving all the fields
AllFields = fieldnames(AllData);

% Neighbouring channel matrix calculation (only if not provided as input)
if isempty(ChN)
    ChN = ept_ChN2(EEG.chanlocs, 0);
end

%% STATISTICAL TEST DECISION

% Determining which ANOVA analysis to run
if isfield(Design,'Between') && isfield(Design,'Within')
    AOVTest = 'Mixed ANOVA';
elseif isfield(Design,'Between') && ~isfield(Design,'Within')
    AOVTest = 'One-Way ANOVA';
elseif ~isfield(Design,'Between') && isfield(Design,'Within')
    AOVTest = 'Repeated-Measures ANOVA';
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(AOVTest,'One-Way ANOVA')
    %% 1) ONE-WAY ANOVA  
    % Restructuring the data
    
    % Between-subject factor(s)
    for k=1:length(Design.Between(3:end))
        
        % Position of the current group in the fields
        Idx = find(contains(AllFields,Design.Between(k+2)));
        
        % For each position
        for p=1:length(Idx)
            
            % Store the data in the big matrix
            Data{k}(:,:,p) = AllData.(AllFields{Idx(p)});
        end
    end
    
    % Concatenate the data
    % SIZE must be electrode x time points x (concatenated) subjects
    FullData = [];
    for k=1:size(Data,2)
        FullData = cat(3,FullData,Data{k});
        NSubj(k) = size(Data{k},3);
    end
    
    %% STATISTICS
    % One-Way ANOVA (using FMUT toolbox)
    [F_Obs, F_dist] = perm_crANOVA(FullData,NSubj, N_Permutes);

    %% Calculating the F values and TFCE enhancement of each different permutation
    % For the original statistics
   
    if TFCE
        try
            TFCE_Obs  = ept_mex_TFCE2D(F_Obs, ChN, E_H);
        catch
            if ~isempty(root_folder)

            % CD to the folder containing the ept_ClustRes.c files 
            TempPWD = pwd;
            cd(root_folder);
            EPTFolder = dir(['**/*' 'ept_TFCE-matlab-master/TFCE/Dependencies']);
            EPTFolder = unique({EPTFolder.folder});

            cd(EPTFolder{1})
            mex -setup
            mex ept_mex_TFCE.c -compatibleArrayDims
            mex ept_mex_TFCE2D.c -compatibleArrayDims
            mex ept_mex_TFCE3D.c -compatibleArrayDims

            % CD back to the Data containing upper folder         
            cd(TempPWD);

            % Re-running TFCE algorithm after compiling C-files
            TFCE_Obs = ept_mex_TFCE2D(F_Obs,  ChN, E_H);

            else
                error('Consider providing a root_folder argument (path of upper folder) to recompile the c-files with MEX');
            end
        end

        % For each permutation
        for i = 1:size(F_dist,1)
            
            if TFCE
            % Display
            fprintf('\n Computing permutation %d',i)

            % Squeezing data for each permutation
            F_Perm = squeeze(F_dist(i,:,:));

            % Running TFCE algorithm
            TFCE_Perm = ept_mex_TFCE2D(F_Perm,  ChN, E_H);
            
            else
                TFCE_Perm = T_Perm;
            end
            
            % stores the maximum absolute value
            maxTFCE(i) = max(abs(TFCE_Perm(:)));    
        end

        % Permutation threshold 
        U = round((1-Pval)*N_Permutes); 
        MaxTFCE=sort(maxTFCE);
        cluster_thresh= MaxTFCE(U);
    
    else
        TFCE_Obs = F_Obs;
    end
    
    % Computing P-Values 
    % add observed maximum
    edges = [maxTFCE';max(abs(TFCE_Obs(:)))];

    [~,bin]     = histc(abs(TFCE_Obs),sort(edges));
    P_Values    = 1-bin./(N_Permutes+2);    
    
     %% Calculating clusters

    % Matrix to return as output
    Results.Obs = F_Obs;
    Results.TFCE_Obs = TFCE_Obs;
    Results.maxTFCE = maxTFCE;
    Results.P_Values = P_Values;

    % Will compile the C-files using MEX if not working
    try
        Cluster_Results = ept_calculateClusters(Results, ChN, Pval);
    catch
        if ~isempty(root_folder)

            % CD to the folder containing the ept_ClustRes.c files 
            TempPWD = pwd;
            cd(root_folder);
            EPTFolder = dir(['**/*' 'ept_TFCE-matlab-master/ResultViewer/Dependencies']);
            cd(EPTFolder.folder)
            
            mex -setup
            mex ept_ClusRes.c -compatibleArrayDims
            mex ept_ClusRes3D.c -compatibleArrayDims
            
            % Re-run clustering after compiling C-files
            Cluster_Results = ept_calculateClusters(Results, ChN, Pval);

            % CD back to the Data containing upper folder         
            cd(TempPWD);
        else
            error('Consider providing a root_folder argument (path of upper folder) to recompile the c-files with MEX');
        end
    end
    
    % TEST PLOT
%     Temp = Results.P_Values;
%     Temp(Temp>0.05) = 0;
%     h=imagesc(Results.P_Values);colorbar;colormap('jet(50)'); hold on; 
%     caxis([min(min(Temp)) max(max(Temp))])
%     set(h,'alphadata',Temp~=0)
%     set(gca,'ydir','norm')
%     for k=1:length(Cluster_Results)
%         contour(Cluster_Results(k).cluster_locations,'linecolor','k');
%     end
%     hold off
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



elseif strcmpi(AOVTest,'Mixed ANOVA')
    %% 2) MIXED ANOVA
    % Restructuring the data

    % Between-subject factor(s)
    for k=1:length(Design.Between(3:end))
        % Within-subject factor(s)
        for l=1:length(Design.Within(3:end))
            % Position of the current group/condition in the fields
            BLogical = contains(AllFields,Design.Between(k+2));
            WLogicaL = contains(AllFields,Design.Within(l+2));
            Idx = find(BLogical+WLogicaL==2);
            % For each position
            for p=1:length(Idx)
                % Store the data in the big matrix
                Data{k}(:,:,l,p) = AllData.(AllFields{Idx(p)});
            end
        end
    end

    % Concatenate the data
    % SIZE must be electrode x time points x Factor A x subjects, where
    % subjects are concatenated !!
    for k=1:size(Data,2)-1
        FullData = cat(4,Data{k},Data{k+1});
    end

    %% STATISTICS
    % Unbalanced samples MIXED ANOVA (using FMUT toolbox)
    % The first F_dist permutation is the same as F_Obs !

    % For each Contrasts (both principal effects + interaction)
    for k=1:3 
        [F_Obs.(StringANOVA{k}), F_dist.(StringANOVA{k})] = perm_spANOVA(FullData,...
            [size(Data{1},4), size(Data{2},4)], Contrasts{k},N_Permutes);
    end

    %% Calculating the F values and TFCE enhancement of each different permutation
    % For the original statistics
    for i=1:3
        if TFCE
            try
                TFCE_Obs.(StringANOVA{i})  = ept_mex_TFCE2D(F_Obs.(StringANOVA{i}),  ChN, E_H);
            catch
                if ~isempty(root_folder)

                % CD to the folder containing the ept_ClustRes.c files 
                TempPWD = pwd;
                cd(root_folder);
                EPTFolder = dir(['**/*' 'ept_TFCE-matlab-master/TFCE/Dependencies']);

                cd(EPTFolder.folder)
                mex -setup
                mex ept_mex_TFCE.c -compatibleArrayDims
                mex ept_mex_TFCE2D.c -compatibleArrayDims
                mex ept_mex_TFCE3D.c -compatibleArrayDims

                % CD back to the Data containing upper folder         
                cd(TempPWD);

                % Re-running TFCE algorithm after compiling C-files
                TFCE_Obs.(StringANOVA{i})  = ept_mex_TFCE2D(F_Obs.(StringANOVA{i}),  ChN, E_H);

                else
                    error('Consider providing a root_folder argument (path of upper folder) to recompile the c-files with MEX');
                end
            end
                
        else
            TFCE_Obs.(StringANOVA{i}) = F_Obs.(StringANOVA{i});
        end
    end

    % For each permutation
    for i = 1:size(F_dist.(StringANOVA{1}),1)

        F_Perm = structfun(@(x) squeeze(x(i,:,:)),F_dist,'UniformOutput',0);

        for k=1:3
            
            if TFCE
                TFCE_Perm.(StringANOVA{k})  = ept_mex_TFCE2D(F_Perm.(StringANOVA{k}),  ChN, E_H);
            else
                TFCE_Perm.(StringANOVA{k}) = F_Perm.(StringANOVA{k});
            end

            % stores the maximum absolute value
            % maxTFCE.(StringANOVA{k})(i,1)  = max(max(TFCE_Perm.(StringANOVA{k})));
            maxTFCE.(StringANOVA{k})(i,1) = max(abs(TFCE_Perm.(StringANOVA{k})(:)));  
        end

        % Display
%         fprintf('\n Computing permutation %d',i)
    end

    % Permutation threshold 
    for i=1:3
        U = round((1-Pval)*N_Permutes); 
        MaxTFCE=sort(maxTFCE.(StringANOVA{i}));
        cluster_thresh(i)= MaxTFCE(U);
    end

    %% Calculating the p value from the permutation distribution
    % Including the observed distribution

    % add observed maximum
    edges.A  = [maxTFCE.A;  max(max(TFCE_Obs.A))]; 
    edges.B  = [maxTFCE.B;  max(max(TFCE_Obs.B))];
    edges.AB = [maxTFCE.AB; max(max(TFCE_Obs.AB))];
    
    
    [~,bin.A]      = histc(abs(F_Obs.A),sort(edges.A));
    P_Values.A     = 1-bin.A./(N_Permutes+2);
    [~,bin.B]      = histc(abs(F_Obs.B),sort(edges.B));
    P_Values.B     = 1-bin.B./(N_Permutes+2);
    [~,bin.AB]     = histc(abs(F_Obs.AB),sort(edges.AB));
    P_Values.AB    = 1-bin.AB./(N_Permutes+2);

    %% Calculating clusters

    for i=1:3
        % Matrix to return as output
        Results.(StringANOVA{i}).Obs = F_Obs.(StringANOVA{i});
        Results.(StringANOVA{i}).TFCE_Obs = TFCE_Obs.(StringANOVA{i});
        Results.(StringANOVA{i}).maxTFCE = maxTFCE.(StringANOVA{i});
        Results.(StringANOVA{i}).P_Values = P_Values.(StringANOVA{i});

        % Will mex the C-files if not working
        try
            [Cluster_Results.(StringANOVA{i})] = ept_calculateClusters(Results.(StringANOVA{i}), ChN, Pval);
        catch
            if ~isempty(root_folder)

                % CD to the folder containing the ept_ClustRes.c files 
                TempPWD = pwd;
                cd(root_folder);
                EPTFolder = dir(['**/*' 'ept_TFCE-matlab-master/ResultViewer/Dependencies']);
                cd(EPTFolder.folder)
                
                mex -setup
                mex ept_ClusRes.c -compatibleArrayDims
                mex ept_ClusRes3D.c -compatibleArrayDims

                % Re-running the clustering algorithm
                [Cluster_Results.(StringANOVA{i})] = ept_calculateClusters(Results.(StringANOVA{i}), ChN, Pval);
                
                % CD back to the Data containing upper folder
                cd(TempPWD);
            else
                error('Consider providing a root_folder argument (path of upper folder) to recompile the c-files with MEX');
            end

           [Cluster_Results.(StringANOVA{i})] = ept_calculateClusters(Results.(StringANOVA{i}), ChN, Pval);
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
elseif strcmpi(AOVTest,'Repeated-Measures ANOVA')
    %% 3) Repeated-Measures ANOVA
    
    % Restructuring the data

    % Between-subject factor(s)
    for l=1:length(Design.Within(3:end))
        % Position of the current group/condition in the fields
        Idx = find(contains(AllFields,Design.Within(l+2)));
        % For each position
        for p=1:length(Idx)
            % Store the data in the big matrix
            % SIZE must be electrode x time points x Factor A x subjects
            FullData(:,:,l,p) = AllData.(AllFields{Idx(p)});
        end
    end
    
     %% STATISTICS
    % Repeated-Measures ANOVA (using FMUT toolbox)
    [F_Obs, F_dist] = perm_rbANOVA(FullData,3, N_Permutes);

    %% Calculating the F values and TFCE enhancement of each different permutation
    % For the original statistics
   
    if TFCE
        try
            TFCE_Obs  = ept_mex_TFCE2D(F_Obs, ChN, E_H);
        catch
            if ~isempty(root_folder)

            % CD to the folder containing the ept_ClustRes.c files 
            TempPWD = pwd;
            cd(root_folder);
            EPTFolder = dir(['**/*' 'ept_TFCE-matlab-master/TFCE/Dependencies']);

            cd(EPTFolder.folder)
            mex -setup
            mex ept_mex_TFCE.c -compatibleArrayDims
            mex ept_mex_TFCE2D.c -compatibleArrayDims
            mex ept_mex_TFCE3D.c -compatibleArrayDims

            % CD back to the Data containing upper folder         
            cd(TempPWD);

            % Re-running TFCE algorithm after compiling C-files
            TFCE_Obs = ept_mex_TFCE2D(F_Obs,  ChN, E_H);

            else
                error('Consider providing a root_folder argument (path of upper folder) to recompile the c-files with MEX');
            end
        end

        % For each permutation
        for i = 1:size(F_dist,1)
            if TFCE
            % Display
            fprintf('\n Computing permutation %d',i)

            % Squeezing data for each permutation
            F_Perm = squeeze(F_dist(i,:,:));

            % Running TFCE algorithm
            TFCE_Perm = ept_mex_TFCE2D(F_Perm,  ChN, E_H);

            else
                TFCE_Perm = F_Perm;
            end

            % stores the maximum absolute value
            % maxTFCE(i,1)  = max(max(TFCE_Perm));
            maxTFCE(i) = max(abs(TFCE_Perm(:))); 
        end
        
        % Permutation threshold 
        U = round((1-Pval)*N_Permutes); 
        MaxTFCE=sort(maxTFCE);
        cluster_thresh= MaxTFCE(U);
        
    else
        TFCE_Obs = F_Obs;
    end
    
    % Computing P-Values 
    % add observed maximum
    edges = [maxTFCE';max(abs(TFCE_Obs(:)))];

    [~,bin]     = histc(abs(TFCE_Obs),sort(edges));
    P_Values    = 1-bin./(N_Permutes+2);    
    
     %% Calculating clusters

    % Matrix to return as output
    Results.Obs = F_Obs;
    Results.TFCE_Obs = TFCE_Obs;
    Results.maxTFCE = maxTFCE;
    Results.P_Values = P_Values;
    Results.FlagTFCE = TFCE;

    % Will compile the C-files using MEX if not working
    try
        Cluster_Results = ept_calculateClusters(Results, ChN, Pval);
    catch
        if ~isempty(root_folder)

            % CD to the folder containing the ept_ClustRes.c files 
            TempPWD = pwd;
            cd(root_folder);
            EPTFolder = dir(['**/*' 'ept_TFCE-matlab-master/ResultViewer/Dependencies']);
            cd(EPTFolder.folder)
            
            mex -setup
            mex ept_ClusRes.c -compatibleArrayDims
            mex ept_ClusRes3D.c -compatibleArrayDims
            
            % Re-run clustering after compiling C-files
            Cluster_Results = ept_calculateClusters(Results, ChN, Pval);

            % CD back to the Data containing upper folder         
            cd(TempPWD);
        else
            error('Consider providing a root_folder argument (path of upper folder) to recompile the c-files with MEX');
        end
    end
    
    % TEST PLOT
%     Temp = Results.P_Values;
%     Temp(Temp>0.05) = 0;
%     h=imagesc(Results.P_Values);colorbar;colormap('jet(50)'); hold on; 
%     caxis([min(min(Temp)) max(max(Temp))])
%     set(h,'alphadata',Temp~=0)
%     set(gca,'ydir','norm')
%     for k=1:length(Cluster_Results)
%         contour(Cluster_Results(k).cluster_locations,'linecolor','k');
%     end
%     hold off
    
end
end