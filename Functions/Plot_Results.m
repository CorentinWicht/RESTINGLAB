% Author: Corentin Wicht, LCNS, 2019
% - corentin.wicht@unifr.ch
% - https://github.com/CorentinWicht

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

%% Function to plot t-test results

function Plot_Results(EEG,Design,Spect,varargin)

%% Set Defaults

Pval            = 0.05; % default threshold for significance
Normalization   = 'Y'; % default  normalisation (= on)
 
% Process Secondary Arguments
if nargin > 3
    for i = 1:2:length(varargin)
        Param = varargin{i};
        Value = varargin{i+1};
        if ~ischar(Param)
          error('Flag arguments must be strings')
        end
        Param = lower(Param);
    
        switch Param
            case 'normalization'
                Normalization  = Value;
            case 'pval'
                Pval           = Value;
            otherwise
                display (['Unknown parameter setting: ' Param])
        end
    end
end


% Thresholding original data with P-val (user defined)
TempP = Spect.TFCE.Statistics.P_Values;
TempP(TempP > Pval) = 0;

% PLOT TEST %%%%%%%%%%%%%%
h =  findobj('type','figure');
figure(length(h)+1); clf
% annotation('textbox',[.35 .5 .5 .5],'String',...
% strrep(Stats.FileName,'_','-'),'FitBoxToText','on','FontSize',20,'EdgeColor','w');
Labels = {EEG.chanlocs.labels};
Dat = Spect.TFCE.Statistics.P_Values;
LogicalDat = TempP;
LogicalDat(TempP>0)=1;
XTStr = round(EEG.xmin*1000):50:round(EEG.xmax*1000,-1);
Ticks = linspace(1,size(Dat,2),size(XTStr,2));
MinVal=min([min(min(min(Spect.data{1}))) min(min(min(Spect.data{2})))]);
MaxVal=max([max(max(max(Spect.data{1}))) max(max(max(Spect.data{2})))]);

% I STOPPED HERE !!!!!!!!!!!

for pp=1:2
    subplot(2,2,pp)
    PlotData = squeeze(mean(permute(TempData.(CompareFN{pp}),[2 3 1]),3));
    imagesc(PlotData)
    title(CurrentDesign{pp+2})
    set(gca,'ydir','norm')
    hcb = colorbar;
    if strcmpi(Normalization,'Y')
        set(get(hcb,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
    else
        set(get(hcb,'XLabel'),'String','Power (\muV)');
    end
    caxis([MinVal MaxVal])
    colormap(gca,'hsv(50)')
    set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
    set(gca,'YTick',pp:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(pp:2:end));
    hold on; contour(logical(LogicalDat),1,'linecolor','k'); hold off;
end

subplot(2,2,3)
h = imagesc(TempP);
set(gca,'ydir','norm')
title('Stats (p-val)')
colorbar
colormap(gca,'hot(50)')
set(h,'alphadata',TempP~=0)
set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
set(gca,'YTick',1:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(1:2:end));

subplot(2,2,4)
DatFDR=Dat;
DatFDR(DatFDR>pN)=0;
hh=imagesc(DatFDR);
set(gca,'ydir','norm')
title('FDR-corrected P-val')
colorbar
colormap(gca,'hot(50)')
set(hh,'alphadata',DatFDR~=0)
set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
set(gca,'YTick',1:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(1:2:end));

end