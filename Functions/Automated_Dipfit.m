function [EEG,DipolesRegion]= Automated_Dipfit(EEG,varargin)
            
% This function performs automated Dipoles Fitting (DIPFIT2) 
% while locating the dipoles in LBPA40 Atlas. 

% Usage:
%    >> [EEG,DipolesRegion]= Automated_Dipfit(EEG)

% Inputs:
%   EEG       = EEG dataset structure

% Optional:
%   grid     - [cell array] 3 cells corresponding to the x,y,z grid values. 
%               Default is 10 elements between -1 and 1.
%   threshold - threshold value for "true" peak selection (default = 35)

% Outputs:
%   EEG       - Updated EEG structure

% Author: Corentin Wicht, LCNS, 2018

%-------------------------------------------------------------------------%
% DIFPIT
%-------------------------------------------------------------------------%
% Defaults
threshold=35;
xgrid=-85:17:85;
ygrid=-85:17:85;
zgrid=0:17:85;

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
            case 'threshold'
                threshold  = Value;
            case 'grid'
                xgrid=Value{1};
                ygrid=Value{2};
                zgrid=Value{3}; 
            otherwise
                display (['Unknown parameter setting: ' Param])
        end
    end
end


%% STARTS Dipole fitting
% Need to use MNI template to match the ATLAS coordinates! 
DipFitPath=backslash(strcat(strrep(which('eeglab'),'eeglab.m',''),'plugins\dipfit3.0\'));

%  Calculate the invidivualized transform parameters of elect locations
EEG=pop_chanedit(EEG, 'lookup',strcat(DipFitPath,'standard_BEM\\elec\\standard_1005.elc')); 

% Setting DipFit models and preferences
EEG = pop_dipfit_settings(EEG, 'hdmfile',...
    strcat(DipFitPath,'standard_BEM\\standard_vol.mat'),...
    'coordformat','MNI','mrifile',strcat(DipFitPath,'standard_BEM\\standard_mri.mat'),... 
    'chanfile',strcat(DipFitPath,'standard_BEM\\elec\\standard_1005.elc'),...
    'coord_transform',[0 0 0 0 0 -1.5708 1 1 1] ,'chansel',1:EEG.nbchan); 

% Initial fitting - Scanning on a coarse-grained grid
EEG = pop_dipfit_gridsearch(EEG, 1:length(EEG.reject.gcompreject) ,...
    xgrid,ygrid,zgrid,0.15); 
% 0.15 parameter is based on: Wyczesany, Grzybowski, & Kaiser, 2015 + Ferdek et al. 2016 + Hammon et al. 2008

% Automated dipole fitting of selected components
EEG = pop_multifit(EEG, 1:length(EEG.reject.gcompreject) ,'threshold',100,...
    'dipplot','off','plotopt',{'normlen' 'on'},'rmout','on');

% Search for and estimate symmetrically constrained bilateral dipoles
% Source: https://link.springer.com/chapter/10.1007%2F978-3-319-32703-7_22
EEG = fitTwoDipoles(EEG, 'LRR', threshold);
close gcf

% Adding LABELS to IC using LBPA40 ATLAS
% Also identifies IC that are probably artifacts (in skull,
% CSF, etc.)
DipolesLoc=[];
SubCompRej=[];
RemainCompPlot=[];
jj=1;
mm=1;
for kk=1:length(EEG.dipfit.model)
    if isempty(EEG.dipfit.model(kk).posxyz)==0 
        DipolesLoc(jj,1:3)=EEG.dipfit.model(kk).posxyz(1,:);
        RemainCompPlot(jj)=kk;
        jj=jj+1;
    else
        SubCompRej(mm)=kk; % Tag for rejection components identified by DIPFIT as artifacts ([] in ...dipfit.model)
        mm=mm+1;
    end
end

% Labeling dipoles
[Pdipoles,labels] = label_dipoles(DipolesLoc);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% Reduce size of labels
labels = strrep(labels,'Avg152T1 ','');
        
% Tag more components to be rejected, which cannot be mapped to a brain region (=ARTIFACTS)
% Homemade method:
DipolesRegion={};
pp=1;
DipLabelMax=max(Pdipoles');
LowFit=[];
pq=1;
for ttt=1:length(Pdipoles(:,1))
     [M,DipolesMaxIndex] = max(Pdipoles(ttt,:));
    if M==0 
        SubCompRej(mm)=RemainCompPlot(ttt);
        mm=mm+1;
    elseif M<mean(DipLabelMax)
        SubCompRej(mm)=RemainCompPlot(ttt);
        mm=mm+1;
        LowFit(pq)=ttt;
        pq=pq+1;
    else
%         DipolesRegion{pp}=strcat(num2str(pp),')--C',num2str(RemainCompPlot(ttt)),...
%             '-',labels{DipolesMaxIndex}); 
%         pp=pp+1;
    end
    DipolesRegion{ttt}=[num2str(ttt) ')' labels{DipolesMaxIndex}]; 
end
SubCompRej=sort(SubCompRej);

% Bar graph of Label Fitting to Dipoles Location
figure(3)
BarG=bar(DipLabelMax);
hold on
hline(mean(DipLabelMax),':k','MEAN')
xlabel(gca,'Components')
set(gca,'XTick',1:size(Pdipoles,1))
ylabel(gca,'Label Fit')
BarG.FaceColor = 'flat';
for xx=1:length(LowFit)
    BarG.CData(LowFit(xx),:) = [1 0 0];
end
title('Label Fitting to Dipoles Location (red bars are lower than mean)');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% % Remaining components to plot
% RemainCompPlot=[];
% qqq=1;
% for fff=1:length(EEG.reject.gcompreject)
%     if ismember(fff,SubCompRej)==0
%         RemainCompPlot(qqq)=fff;
%         qqq=qqq+1;
%     end
% end

% % Updating IC list
% for mmm=1:length(SubCompRej)
%     EEG.reject.gcompreject(SubCompRej(mmm))=1;
% end
% Updating IC 
% EEG = pop_subcomp( EEG,[],0); 

% Plot dipoles
dipplot(EEG.dipfit.model(RemainCompPlot),'color', { 'b' },'mri',...
    strcat(DipFitPath,'standard_BESA\\avg152t1.mat'),'summary',...
    'on','num','off','dipnames',DipolesRegion,'normlen','on');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% pop_dipplot( EEG,RemainCompPlot ,'color', { 'b' },'mri',...
%     strcat(DipFitPath,'standard_BESA\\avg152t1.mat'),'summary',...
%     'on','dipnames',DipolesRegion,'num','on','normlen','on');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

end