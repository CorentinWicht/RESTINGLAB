classdef EEGBrowser < CoreBrowser
    properties
        windowWidth = 5;
        normalizeFlag = false;
        % showChannelNumber = false;
        gain = 0.5;
        numberOfChannelsToPlot
        yTickLabel = [];
        yTick = [];
        colormap = 'eegplot';
        colorInCell
        textHandle;
        isEpoched
        dim
        trialSelection = [];
    end
    properties(SetObservable)
        color
        showChannelNumber = false;
    end
    properties(SetAccess=private)
        hRec = [];
    end
    methods
        %% constructor
        function obj = EEGBrowser(EEG, objMaster)
            if nargin < 2, objMaster = -1;end
            obj.dim = size(EEG.data);
            obj.isEpoched = false;
            if length(obj.dim) > 2
                obj.isEpoched = true;
            end
            if obj.isEpoched
                EEG.data = reshape(EEG.data,[obj.dim(1) prod(obj.dim(2:3))]);
            end
            ntimePoints = size(EEG.data,2);
            obj.streamHandle.name = EEG.setname;
            obj.streamHandle.timeStamp = (0:ntimePoints-1)/EEG.srate;
            obj.streamHandle.numberOfChannels = obj.dim(1);
            obj.streamHandle.mmfName = [tempname '.bin'];
            precision = class(EEG.data(1));
            fid = fopen(obj.streamHandle.mmfName,'w');
            delta = ceil(0.2*obj.streamHandle.numberOfChannels);
            for it=1:delta:obj.streamHandle.numberOfChannels
                if it+delta-1 <=obj.streamHandle.numberOfChannels
                    fwrite(fid,EEG.data(it:it+delta-1,:)',precision);
                else
                    fwrite(fid,EEG.data(it:end,:)',precision);
                end
            end
            fclose(fid);
            obj.streamHandle.mmfObj = memmapfile(obj.streamHandle.mmfName,'Format',{precision [length(obj.streamHandle.timeStamp) obj.streamHandle.numberOfChannels] 'x'},'Writable',false);
            % obj.streamHandle.data = EEG.data(:,1:end)';
            obj.streamHandle.event = event;
            obj.streamHandle.samplingRate = EEG.srate;
            try
                obj.streamHandle.label = {EEG.chanlocs.labels};
            catch 
                for it=1:obj.streamHandle.numberOfChannels, obj.streamHandle.label{it} = num2str(it);end
            end
            while length(obj.streamHandle.label) < obj.streamHandle.numberOfChannels
                obj.streamHandle.label{end+1} = num2str(length(obj.streamHandle.label)+1);
            end
            obj.streamHandle.chanlocs = EEG.chanlocs;
            obj.streamHandle.icawinv = EEG.icawinv;
            obj.streamHandle.icachansind = EEG.icachansind;
            if isfield(EEG.event,'type')
                type = {EEG.event.type};
                for ev=1:length(type), if isnumeric(type{ev}), type{ev} = num2str(type{ev});end;end
                latency = cell2mat({EEG.event.latency});
                if obj.isEpoched
                    latency = [latency [1:obj.dim(2):prod(obj.dim(2:3)) prod(obj.dim(2:3))]];
                    type = [type repmat({''},1,obj.dim(3)+1)];
                end
                obj.streamHandle.event = obj.streamHandle.event.addEvent(latency,type);
            end
            obj.addlistener('channelIndex','PostSet',@EEGBrowser.updateChannelDependencies);
            obj.addlistener('showChannelNumber','PostSet',@EEGBrowser.updateChannelDependencies);
            obj.addlistener('timeIndex','PostSet',@EEGBrowser.updateTimeIndexDenpendencies);
            obj.addlistener('color','PostSet',@EEGBrowser.updateColorInCell);
            obj.master = objMaster;
            obj.init();
        end
        %% plot
        function init(obj)
            
            % Do not re-create the gui if it already exists
            if isempty(obj.figureHandle) || ~isvalid(obj.figureHandle)
                init@CoreBrowser(obj);
                set(obj.figureHandle,'Renderer','Painters','CloseRequestFcn',@(src, event)onClose(obj,[], event))
                set(obj.figureHandle,'name','Scroll activities -- pop_eegbrowser()');
                        
                %- topoplot on ButtonDown event
                % set(obj.gObjHandle,'ButtonDownFcn',@gObject_ButtonDownFcn);
                resFolder = [fileparts(which('EEGBrowser.m')) filesep 'resources'];
                imgScaleDown = imread([resFolder filesep 'Gnome-list-remove.svg.png']);
                imgScaleUp = imread([resFolder filesep 'Gnome-list-add.svg.png']);
                helpIcon = imread([resFolder filesep 'Gnome-help-browser.svg.png']);
                selectIcon = imread([resFolder filesep 'selectTrial.png']);
                hDown = uicontrol('Parent', obj.figureHandle, 'Style', 'pushbutton','TooltipString','Scale down','Position',[543+150 53 40 40],'Callback',@onScaleDown,'CData',imgScaleDown);
                hUp = uicontrol('Parent', obj.figureHandle, 'Style', 'pushbutton','TooltipString','Scale up','Position',[543+190 53 40 40],'Callback',@onScaleUp,'CData',imgScaleUp);
                hHelp = uicontrol(obj.figureHandle,'CData',helpIcon,'TooltipString','Help','Position',[543+230 53 40 40],'Callback','web(''https://github.com/aojeda/EEGBrowser#eegbrowser'',''-browser'')');
                set([hDown, hUp hHelp],'Units','Normalized')
                set(obj.figureHandle,'KeyPressFcn',@onKeyPress);
                toolbarHandle = findall(obj.figureHandle,'Type','uitoolbar');
                hcb = uitoggletool(toolbarHandle,'CData',selectIcon,'Separator','on','HandleVisibility','off','TooltipString','Select trial','State','off');
                set(hcb,'OnCallback',@(src,event)enableTrialSelection(obj,hcb,'select'),'OffCallback',@(src, event)enableTrialSelection(obj,hcb,'remove'));
                set(obj.axesHandle,'ButtonDownFcn',@selectTrial);
            end
            
            % Find now cursor index
            [~,t1] = min(abs(obj.streamHandle.timeStamp(obj.timeIndex) - (obj.nowCursor-obj.windowWidth/2)));
            [~,t2] = min(abs(obj.streamHandle.timeStamp(obj.timeIndex) - (obj.nowCursor+obj.windowWidth/2)));
            data = obj.streamHandle.mmfObj.Data.x(obj.timeIndex(t1:t2),obj.channelIndex);
            cla(obj.axesHandle);
            hold(obj.axesHandle,'on');
            obj.gObjHandle = plot(obj.axesHandle,obj.streamHandle.timeStamp(obj.timeIndex(t1:t2)),data);
            for it=1:obj.numberOfChannelsToPlot
                set(obj.gObjHandle(it),'color',obj.color(it,:),'userData',{obj.streamHandle, obj.channelIndex(it)});
            end
            
            sigma = nanstd(obj.streamHandle.mmfObj.Data.x);
            if obj.numberOfChannelsToPlot > 1
                obj.yTick = (1:obj.numberOfChannelsToPlot)*mean(sigma)/obj.gain;
                delta = abs(diff(obj.yTick([2 1])));
                lim = [obj.yTick(1) - delta obj.yTick(end) + delta];
            else
                obj.yTick = mean(data);
                mx = 1.5*max(abs(data));
                lim = obj.yTick+[-mx mx];
            end
            set(obj.axesHandle,'YTickLabel',obj.yTickLabel,'YTick',obj.yTick,'YLim',lim);
            obj.plotThisTimeStamp(obj.nowCursor);
        end
        %%
        function plotThisTimeStamp(obj,nowCursor)
            delta = obj.windowWidth/2;
            if  nowCursor + delta < obj.streamHandle.timeStamp(obj.timeIndex(end)) &&...
                    nowCursor - delta > obj.streamHandle.timeStamp(obj.timeIndex(1))
                newNowCursor = nowCursor;
            elseif nowCursor + delta >= obj.streamHandle.timeStamp(obj.timeIndex(end))
                newNowCursor = obj.streamHandle.timeStamp(obj.timeIndex(end)) - delta;
                if strcmp(get(obj.timerObj,'Running'),'on')
                    stop(obj.timerObj);
                end
            else
                newNowCursor = obj.streamHandle.timeStamp(obj.timeIndex(1)) + delta;
            end
            nowCursor = newNowCursor;
            
            % find now cursor index
            obj.nowCursor = nowCursor;
            [~,t1] = min(abs(obj.streamHandle.timeStamp(obj.timeIndex) - (obj.nowCursor-obj.windowWidth/2)));  
            [~,t2] = min(abs(obj.streamHandle.timeStamp(obj.timeIndex) - (obj.nowCursor+obj.windowWidth/2)));  
            if t1==t2, return;end
            dt = length(t1:t2)/2;
            data = obj.streamHandle.mmfObj.Data.x(obj.timeIndex(t1:t2),obj.channelIndex);
            data = obj.gain*data + ones(2*dt,1)*fliplr(obj.yTick);
            if obj.numberOfChannelsToPlot <= 1
                data = {data};
            elseif obj.numberOfChannelsToPlot == 2
                data = num2cell(data,1)';
            else
                data = num2cell(data,[1 obj.numberOfChannelsToPlot])';
            end
            set(obj.gObjHandle,'XData',obj.streamHandle.timeStamp(obj.timeIndex(t1:t2)),{'YData'},data,{'Color'},obj.colorInCell);           
            xlim(obj.axesHandle,obj.streamHandle.timeStamp(obj.timeIndex([t1 t2])));
            
            if obj.showEvents
                if length(obj.eventObj.latencyInFrame) ~= size(obj.eventColor,1), obj.initEventColor;end
                loc = find(obj.eventLatencyLookUp(obj.eventObj.latencyInFrame)>obj.streamHandle.timeStamp(obj.timeIndex(t1)) &...
                    obj.eventLatencyLookUp(obj.eventObj.latencyInFrame)<obj.streamHandle.timeStamp(obj.timeIndex(t2)));
                if ~isempty(loc)
                    lim = ylim();
                    Nloc = length(loc);
                    hold(obj.axesHandle,'on');
                    set(obj.figureHandle,'CurrentAxes',obj.axesHandle)
                    linesHandler = line(ones(2,1)*obj.eventLatencyLookUp(obj.eventObj.latencyInFrame(loc)),(ones(length(loc),1)*lim)','Parent',obj.axesHandle);
                    textPos = [obj.eventLatencyLookUp(obj.eventObj.latencyInFrame(loc));ones(1,Nloc)*lim(2)*1.01];
                    try delete(obj.textHandle);end %#ok
                    set(linesHandler,{'color'},num2cell(obj.eventColor(loc,:)',[1 3])');
                    obj.textHandle = zeros(length(loc),1);
                    for it=1:Nloc
                        obj.textHandle(it) = text('Position',textPos(:,it),'String',obj.eventObj.label(loc(it)),'Color',obj.eventColor(loc(it),:),...
                            'Parent',obj.axesHandle,'FontSize',9,'FontWeight','bold','Rotation',45);
                    end
                    boundaryEvents = cellfun(@isempty,obj.eventObj.label(loc));
                    set(linesHandler(boundaryEvents),'LineStyle','-.','Color','k','LineWidth',2);
                    set(obj.figureHandle,'CurrentAxes',obj.axesHandle)
                    hold(obj.axesHandle,'off');
                end
            end
            set(obj.timeTexttHandle,'String',['Current latency = ' num2str(obj.nowCursor,4) ' sec']);
            set(obj.sliderHandle,'Value',obj.nowCursor);
        end
        %%
        function plotStep(obj,step)
            delta = obj.windowWidth/2;
            if obj.nowCursor+step+delta < obj.streamHandle.timeStamp(obj.timeIndex(end)) &&...
                    obj.nowCursor+step-delta > obj.streamHandle.timeStamp(obj.timeIndex(1))
                newNowCursor = obj.nowCursor+step;
            elseif obj.nowCursor+step+delta > obj.streamHandle.timeStamp(obj.timeIndex(end))
                newNowCursor = obj.streamHandle.timeStamp(obj.timeIndex(end))-delta;
            else
                newNowCursor = obj.streamHandle.timeStamp(obj.timeIndex(1))+delta;
            end
            obj.plotThisTimeStamp(newNowCursor);
        end
        %%
        function changeColormap(obj,newColormap)
            if nargin < 2, newColormap = '';end
            cmaps = dir(fileparts(which('colormap')));
            if strcmpi(newColormap,'eegplot')
                obj.color = ones(obj.numberOfChannelsToPlot,1)*[0 0 0.4];
            elseif any(ismember({cmaps.name},[newColormap '.m']))
                obj.color = eval([newColormap '(' num2str(obj.numberOfChannelsToPlot) ')']);
            else
                warndlg(['Colormap ''' newColormap ''' is not available. We will use ''lines'' instead.'])
                obj.color = lines(obj.numberOfChannelsToPlot);
                newColormap = 'lines';
            end
            obj.colormap = newColormap;
        end
        %%
        function onClose(obj,~,~)
            delete(obj);
        end
        function delete(obj)
            delete(obj.figureHandle);
            delete(obj.streamHandle.mmfName);
        end
        
        %%
        function obj = changeSettings(obj)
            prompt = {'Scale (use +/- keys to change the scale only)','Channels to plot','Speed (values between 1 an 5)','Page width (in seconds)','Normalize','Show channel number','Show events','Colormap (any MATLAB colormap or eegplot for eegplot color)'};
            defaultVal = {num2str(obj.gain), ['[' num2str(obj.channelIndex) ']'], num2str(obj.speed), num2str(obj.windowWidth), num2str(obj.normalizeFlag),...
                num2str(obj.showChannelNumber), num2str(obj.showEvents), obj.colormap};
            properties = inputdlg(prompt,'Preferences',1,defaultVal);
            
            if isempty(properties)
                return;
            end
            obj.gain = abs(str2num(properties{1}));                 %#ok
            obj.channelIndex = str2num(properties{2});              %#ok
            tmp = str2double(properties{3});
            obj.speed = interp1(1:5,1:5,tmp,'nearest','extrap');
            obj.windowWidth = abs(str2num(properties{4}));
            if obj.windowWidth > obj.streamHandle.timeStamp(obj.timeIndex(end))
                obj.windowWidth = obj.streamHandle.timeStamp(obj.timeIndex(end));
            end
            obj.normalizeFlag = logical(str2num(properties{5}));    %#ok
            obj.showChannelNumber = logical(str2num(properties{6}));%#ok
            obj.showEvents = str2num(properties{7});                %#ok
            obj.changeColormap(properties{8});
            
            figure(obj.figureHandle);
            obj.init();
        end
    end
    %%
    methods(Static)
        function updateChannelDependencies(~,evnt)
            evnt.AffectedObject.numberOfChannelsToPlot = length(evnt.AffectedObject.channelIndex);
            
            evnt.AffectedObject.yTickLabel = cell(evnt.AffectedObject.numberOfChannelsToPlot,1);
            labels = cell(evnt.AffectedObject.numberOfChannelsToPlot,1);
            channels = evnt.AffectedObject.channelIndex;
            
            if evnt.AffectedObject.showChannelNumber
                if evnt.AffectedObject.numberOfChannelsToPlot > 1
                    for jt=fliplr(1:evnt.AffectedObject.numberOfChannelsToPlot), labels{jt} = num2str(channels(evnt.AffectedObject.numberOfChannelsToPlot-jt+1));end
                else
                    labels{1} = num2str(channels);
                end
            else
                if evnt.AffectedObject.numberOfChannelsToPlot > 1
                    if isempty(evnt.AffectedObject.streamHandle.label)
                        for jt=fliplr(1:evnt.AffectedObject.numberOfChannelsToPlot), labels{jt} = num2str(evnt.AffectedObject.numberOfChannelsToPlot-jt+1);end
                    else
                        for jt=fliplr(1:evnt.AffectedObject.numberOfChannelsToPlot), labels{jt} = evnt.AffectedObject.streamHandle.label{channels(evnt.AffectedObject.numberOfChannelsToPlot-jt+1)};end
                    end
                else
                    if ~isempty(evnt.AffectedObject.streamHandle.label)
                        labels = evnt.AffectedObject.streamHandle.label(evnt.AffectedObject.channelIndex);
                    else
                        labels{1} = '';
                    end
                end
            end
            evnt.AffectedObject.yTickLabel = labels;
            evnt.AffectedObject.changeColormap(evnt.AffectedObject.colormap);
        end
        %%
        function updateTimeIndexDenpendencies(~,evnt)
            if evnt.AffectedObject.timeIndex(1) ~= -1
                evnt.AffectedObject.nowCursor = evnt.AffectedObject.streamHandle.timeStamp(evnt.AffectedObject.timeIndex(1)) + 2.5;
                set(evnt.AffectedObject.sliderHandle,'Min',evnt.AffectedObject.streamHandle.timeStamp(evnt.AffectedObject.timeIndex(1)));
                set(evnt.AffectedObject.sliderHandle,'Max',evnt.AffectedObject.streamHandle.timeStamp(evnt.AffectedObject.timeIndex(end)));
                set(evnt.AffectedObject.sliderHandle,'Value',evnt.AffectedObject.nowCursor);
            end
        end
        %%
        function updateColorInCell(~,evnt)
            evnt.AffectedObject.colorInCell = cell(evnt.AffectedObject.numberOfChannelsToPlot,1);
            for it=1:evnt.AffectedObject.numberOfChannelsToPlot
                evnt.AffectedObject.colorInCell{it} = evnt.AffectedObject.color(it,:);
            end
        end
    end
end

%%
function onKeyPress(src,evnt)
obj = src.UserData;
switch evnt.Key
    case 'leftarrow',  plotStep(obj,-obj.step*obj.speed*2);
    case 'rightarrow', plotStep(obj,obj.step*obj.speed*2);
    case 'subtract'
        obj.gain = obj.gain/2;
        plotStep(obj,0);
    case 'add'
        obj.gain = obj.gain*2;
        plotStep(obj,0);
end
end
%%
function onScaleDown(src,evnt)
obj = src.Parent.UserData;
obj.gain = obj.gain/2;
plotStep(obj,0);
end

function onScaleUp(src,evnt)
obj = src.Parent.UserData;
obj.gain = obj.gain*2;
plotStep(obj,0);
end

function enableTrialSelection(obj,evnt,state)
if strcmp(state,'select')
    set(obj.axesHandle,'ButtonDownFcn',@selectTrial);
else
    set(obj.axesHandle,'ButtonDownFcn','');
end
end

function selectTrial(src,evnt)
obj = src.Parent.UserData;
if ~obj.isEpoched
    return;
end
boundary = find(cellfun(@isempty,obj.eventObj.label));
b1 = find(obj.eventLatencyLookUp(obj.eventObj.latencyInFrame(boundary)) < evnt.IntersectionPoint(1),1,'last');
if isempty(b1)
    b1 = 1;
    trial = 1;
else
    trial = b1;
end
if any(ismember(obj.trialSelection,trial))
    return;
end
b2 = b1+1;
obj.trialSelection = sort([obj.trialSelection trial]);
assignin('base','trialSelection',obj.trialSelection);
disp(['Selection: ' num2str(obj.trialSelection)]);
x = obj.eventLatencyLookUp(obj.eventObj.latencyInFrame(boundary(b1)));
w = obj.eventLatencyLookUp(obj.eventObj.latencyInFrame(boundary(b2)));
y = obj.axesHandle.YLim(1);
h = obj.axesHandle.YLim(2);
patch('Parent',obj.axesHandle,'XData',[x w w x],'YData',[y y h h],'FaceColor',[0.75 0.75 0.75],...
    'FaceAlpha',0.5,'lineStyle','none','UserData', trial,'ButtonDownFcn',@removeSel);
end

function removeSel(src,evnt)
obj = src.Parent.Parent.UserData;
obj.trialSelection(obj.trialSelection==src.UserData) = [];
assignin('base','trialSelection',obj.trialSelection);
disp(['Selection: ' num2str(obj.trialSelection)]);
delete(src);
end
