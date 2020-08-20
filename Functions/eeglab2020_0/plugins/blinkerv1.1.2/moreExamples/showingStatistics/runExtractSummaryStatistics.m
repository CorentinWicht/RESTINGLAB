%% Extract the summary statistics for a collection of datasets 
%
% This script assumes that a filelist is provided and that the blink
% structures have already been computed using a function such as
% pop_blinker.  The goal is to put all of the statistics for the collection
% into a single file for display and analysis purposes.
%
% BLINKER extracts blinks and ocular indices from time series. 
% Copyright (C) 2016  Kay A. Robbins, Kelly Kleifgas, UTSA
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Setup to generate summary for the VEP data
% blinkDirInd = 'O:\ARL_Data\VEP\BlinkOutput\AllUnRef';
% blinkFileList = 'O:\ARL_Data\VEP\BlinkOutput\vep_blinkFileInfo';
% typeBlinks = 'AllUnRef';
% excludedTasks = {};
% summaryFile = 'O:\ARL_Data\VEP\BlinkOutput\vep_oddball_ALLUnRef_summary.mat';

%% Setup to generate summary for the Shooter combined data
% blinkDirInd = 'O:\ARL_Data\Shooter\BlinkOutput\AllMastRefCombined';
% blinkFileList = 'O:\ARL_Data\Shooter\BlinkOutput\shooter_blinkFileInfo';
% typeBlinks = 'AllMastRefCombined';
% excludedTasks = {'EC', 'EO'};
% summaryFile = 'O:\ARL_Data\Shooter\BlinkOutput\shooter_AllMastRefCombined_summary.mat';

%% Setup to generate summary for the NCTU_RWN_VDE combined data
experiment = 'NCTU_RWN_VDE';
blinkDir = 'J:\CTAData\NCTU_RWN_VDE\Blinks';
blinkDirInd = [blinkDir filesep 'AllMastRefCombinedWithDate'];
typeBlinks = 'AllMastRefCombined';
summaryFile = [blinkDir filesep experiment ...
              '_AllMastRefCombinedWithDate_summary.mat'];
excludedTasks = {'Pre_EXP_resting', 'Post_EXP_resting'};
blinkFileList = [blinkDir filesep experiment '_blinkFileList.mat'];

%% Load the blink file information and find the information
load(blinkFileList);
numberFiles = length(blinkFiles);
[blinkFilePaths, fileMask] = getBlinkFilePaths(blinkDirInd, blinkFiles, ...
    typeBlinks, excludedTasks);

%% Fill in an empty structure for efficiency
blinkStatisticsSummary(numberFiles) = extractBlinkStatistics();
for k = 1:numberFiles - 1
    blinkStatisticsSummary(k) = extractBlinkStatistics();
end

%% Now read in the individual files and process
mapGood = containers.Map('KeyType', 'char', 'ValueType', 'any');
mapMarginal = containers.Map('KeyType', 'char', 'ValueType', 'any');
nanMask = false(numberFiles, 1);

for k = 1:numberFiles
    clear blinks blinkFits blinkProperties blinkStatistics params;
    if ~fileMask(k)
        continue;
    end

    fprintf('Loading %s...\n', blinkFilePaths{k});
    load (blinkFilePaths{k});
    if ~exist('blinks', 'var')
        fileMask(k) = false;
        warning('---%s does not contain blink structures\n', blinkFilePaths{k});
        continue;
    elseif isnan(blinks.usedSignal) || isempty(blinks.usedSignal)
        nanMask(k) = true;
        warning('---%s does not have blinks\n', blinkFilePaths{k});
        continue;
    elseif ~exist('blinkStatistics', 'var') || isempty(blinkStatistics)
        nanMask(k) = true;
        warning('---%s does not have blinkStatistics\n', blinkFilePaths{k});
        continue;
    end
    blinkStatisticsSummary(k) = blinkStatistics;
    sData = blinks.signalData;
    signalNumbers = cellfun(@double, {sData.signalNumber});
    pos = find(signalNumbers == abs(blinks.usedSignal), 1, 'first');
    theLabel = lower(sData(pos).signalLabel);
    if strcmpi(blinkStatisticsSummary(k).status, 'marginal')
        if isKey(mapMarginal, theLabel)
            theCount = mapMarginal(theLabel);
        else
            theCount = 0;
        end
        theCount = theCount + 1;
        mapMarginal(theLabel) = theCount;    
    elseif strcmpi(blinkStatisticsSummary(k).status, 'good')
        if isKey(mapGood, theLabel)
            theCount = mapGood(theLabel);
        else
            theCount = 0;
        end
        theCount = theCount + 1;
        mapGood(theLabel) = theCount;
        blinkStatisticsSummary(k).status = 'good';
    end
end

%% Save the file for later analysis
save(summaryFile, 'blinkStatisticsSummary', 'mapGood', 'mapMarginal', ...
    'nanMask', 'fileMask', '-v7.3');
