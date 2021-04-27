%% This script shows how to add blink events to an EEG data structure

%% Set up the directories
blinkFileDir = 'J:\CTAData\NCTU_RWN_VDE\Blinks\AllMastRefCombinedWithDate';
blinkFileList = 'J:\CTAData\NCTU_RWN_VDE\NCTU_RWN_VDE_BlinkFileList.mat';
blinkType = 'AllMastRefCombined';
EEGDir = 'J:\CTAData\NCTU_RWN_VDE\Level1\session';
fieldList = {'maxFrame', 'leftZero', 'rightZero', 'leftBase', 'rightBase', ...
             'leftZeroHalfHeight', 'rightZeroHalfHeight'};
%% Load the blink file list
load(blinkFileList);

%% Pick out a dataset from the blinkFiles and load corresponding EEG
blinkFileListPos = 9;
session = blinkFiles(blinkFileListPos).session;
fileName = blinkFiles(blinkFileListPos).fileName;
blinkName = [blinkFiles(blinkFileListPos).blinkFileName '_' blinkType '.mat'];
EEGFile = [EEGDir filesep num2str(session) filesep fileName];
EEG = pop_loadset(EEGFile);
test = load([blinkFileDir filesep blinkName]);

%% Add blink events to the EEG
[EEGNew, blinkSignal] = ...
    addBlinkEvents(EEG, test.blinks, test.blinkFits, ...
                   test.blinkProperties, fieldList);
