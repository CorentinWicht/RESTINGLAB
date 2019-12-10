% clean_rawdata(): a wrapper for EEGLAB to call Christian's clean_artifacts.
%
% Usage:
%   >>  EEG = clean_rawdata(EEG, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst, arg_window, varargin)
%
% ---------------- This is a copy from clean_artifacts --------------------
%
% This function removes flatline channels, low-frequency drifts, noisy channels, short-time bursts
% and incompletely repaird segments from the data. Tip: Any of the core parameters can also be
% passed in as [] to use the respective default of the underlying functions, or as 'off' to disable
% it entirely.
%
% Hopefully parameter tuning should be the exception when using this function -- however, there are
% 3 parameters governing how aggressively bad channels, bursts, and irrecoverable time windows are
% being removed, plus several detail parameters that only need tuning under special circumstances.
%
%   FlatlineCriterion: Maximum tolerated flatline duration. In seconds. If a channel has a longer
%                      flatline than this, it will be considered abnormal. Default: 5
%
%   Highpass :         Transition band for the initial high-pass filter in Hz. This is formatted as
%                      [transition-start, transition-end]. Default: [0.25 0.75].
%
%   ChannelCriterion : Minimum channel correlation. If a channel is correlated at less than this
%                      value to a reconstruction of it based on other channels, it is considered
%                      abnormal in the given time window. This method requires that channel
%                      locations are available and roughly correct; otherwise a fallback criterion
%                      will be used. (default: 0.85)
%
%   LineNoiseCriterion : If a channel has more line noise relative to its signal than this value, in
%                        standard deviations based on the total channel population, it is considered
%                        abnormal. (default: 4)
%
%   BurstCriterion : Standard deviation cutoff for removal of bursts (via ASR). Data portions whose
%                    variance is larger than this threshold relative to the calibration data are
%                    considered missing data and will be removed. According to Chang et al. (2018).
%                    "Evaluation of Artifact Subspace Reconstruction for Automatic EEG Artifact Removal.
%                    Conf Proc IEEE Eng Med Biol Soc. 2018", the recommended value here is 10-100.
%                    For more detail, see https://sccn.ucsd.edu/wiki/Artifact_Subspace_Reconstruction_(ASR)#Comments_to_the_HAPPE_paper.2C_and_how_to_choose_the_critical_parameters
%                    I put the default value of 20 here.
%
%                    (Original description: separated by Makoto, 03/26/2019) The most aggressive value that can
%                    be used without losing much EEG is 3. For new users it is recommended to at
%                    first visually inspect the difference between the original and cleaned data to
%                    get a sense of the removed content at various levels. A quite conservative
%                    value is 5. Default: 5.
%
%   WindowCriterion :  Criterion for removing time windows that were not repaired completely. This may
%                      happen if the artifact in a window was composed of too many simultaneous
%                      uncorrelated sources (for example, extreme movements such as jumps). This is
%                      the maximum fraction of contaminated channels that are tolerated in the final
%                      output data for each considered window. Generally a lower value makes the
%                      criterion more aggressive. Default: 0.25. Reasonable range: 0.05 (very
%                      aggressive) to 0.3 (very lax).
%
%
%
% --------------- Optional inputs that goes for varargin ------------------
%
%   NOTE: The following are detail parameters that may be tuned if one of the criteria does
%   not seem to be doing the right thing. These basically amount to side assumptions about the
%   data that usually do not change much across recordings, but sometimes do.
%
%   'ChannelCriterionMaxBadTime': This is the maximum tolerated fraction of the recording duration 
%                                 during which a channel may be flagged as "bad" without being
%                                 removed altogether. Generally a lower (shorter) value makes the
%                                 criterion more aggresive. Reasonable range: 0.15 (very aggressive)
%                                 to 0.6 (very lax). Default: 0.5.
%
%   'BurstCriterionRefMaxBadChns': If a number is passed in here, the ASR method will be calibrated based 
%                                  on sufficiently clean data that is extracted first from the
%                                  recording that is then processed with ASR. This number is the
%                                  maximum tolerated fraction of "bad" channels within a given time
%                                  window of the recording that is considered acceptable for use as
%                                  calibration data. Any data windows within the tolerance range are
%                                  then used for calibrating the threshold statistics. Instead of a
%                                  number one may also directly pass in a data set that contains
%                                  calibration data (for example a minute of resting EEG).
%
%                                  If this is set to 'off', all data is used for calibration. This will 
%                                  work as long as the fraction of contaminated data is lower than the
%                                  the breakdown point of the robust statistics in the ASR
%                                  calibration (50%, where 30% of clearly recognizable artifacts is a
%                                  better estimate of the practical breakdown point).
%
%                                  A lower value makes this criterion more aggressive. Reasonable
%                                  range: 0.05 (very aggressive) to 0.3 (quite lax). If you have lots
%                                  of little glitches in a few channels that don't get entirely
%                                  cleaned you might want to reduce this number so that they don't go
%                                  into the calibration data. Default: 0.075.
%
%   'BurstCriterionRefTolerances': These are the power tolerances outside of which a channel in a
%                                  given time window is considered "bad", in standard deviations
%                                  relative to a robust EEG power distribution (lower and upper
%                                  bound). Together with the previous parameter this determines how
%                                  ASR calibration data is be extracted from a recording. Can also be
%                                  specified as 'off' to achieve the same effect as in the previous
%                                  parameter. Default: [-3.5 5.5].
%
%   'WindowCriterionTolerances': These are the power tolerances outside of which a channel in the final
%                                output data is considered "bad", in standard deviations relative
%                                to a robust EEG power distribution (lower and upper bound). Any time
%                                window in the final (repaired) output which has more than the
%                                tolerated fraction (set by the WindowCriterion parameter) of channel
%                                with a power outside of this range will be considered incompletely 
%                                repaired and will be removed from the output. This last stage can be
%                                skipped either by setting the WindowCriterion to 'off' or by taking
%                                the third output of this processing function (which does not include
%                                the last stage). Default: [-3.5 7].
%
%   'FlatlineCriterion': Maximum tolerated flatline duration. In seconds. If a channel has a longer
%                        flatline than this, it will be considered abnormal. Default: 5
%
%   'NoLocsChannelCriterion': Criterion for removing bad channels when no channel locations are
%                             present. This is a minimum correlation value that a given channel must
%                             have w.r.t. a fraction of other channels. A higher value makes the
%                             criterion more aggressive. Reasonable range: 0.4 (very lax) - 0.6
%                             (quite aggressive). Default: 0.45.
%
%   'NoLocsChannelCriterionExcluded': The fraction of channels that must be sufficiently correlated with
%                                     a given channel for it to be considered "good" in a given time
%                                     window. Applies only to the NoLocsChannelCriterion. This adds
%                                     robustness against pairs of channels that are shorted or other
%                                     that are disconnected but record the same noise process.
%                                     Reasonable range: 0.1 (fairly lax) to 0.3 (very aggressive);
%                                     note that increasing this value requires the ChannelCriterion
%                                     to be relaxed in order to maintain the same overall amount of
%                                     removed channels. Default: 0.1.
%
%
%
% ----------------- Optional inputs for clean_flatlines() -----------------
%
%   'MaxAllowedJitter' : Maximum tolerated jitter during flatlines. As a multiple of epsilon.
%                        Default: 20
%
%
%
% ----------------- Optional inputs for clean_drifts() -----------------
%
%   'Attenuation' : stop-band attenuation, in db (default: 80)
%
%
%
% ----------------- Optional inputs for celan_channels() ------------------
%   
%   'CleanChannelsWindowLength' : Length of the windows (in seconds) for which correlation is computed; ideally
%                                 short enough to reasonably capture periods of global artifacts or intermittent 
%                                 sensor dropouts, but not shorter (for statistical reasons). Default: 5.
%
%   'NumSamples' : Number of RANSAC samples. This is the number of samples to generate in the random
%                  sampling consensus process. The larger this value, the more robust but also slower 
%                  the processing will be. Default: 50.
%
%   'SubsetSize' : Subset size. This is the size of the channel subsets to use for robust reconstruction, 
%                  as a fraction of the total number of channels. Default: 0.25.
%
%
%
% --------------- Optional inputs for celan_channels_nolocs() -------------
%
%   'LineNoiseAware' : Whether the operation should be performed in a line-noise aware manner. If enabled,
%                      the correlation measure will not be affected by the presence or absence of line 
%                      noise (using a temporary notch filter). Default: true.
%
%
%
% -------------------- Optional inputs for clean_asr() --------------------
%
%   The following are detail parameters that usually do not have to be tuned. If you cannot get
%   the function to do what you want, you might consider adapting these better to your data.
%
%   'ASR_WindowLength' : Length of the statistcs window, in seconds. This should not be much longer 
%                        than the time scale over which artifacts persist, but the number of samples in
%                        the window should not be smaller than 1.5x the number of channels. Default:
%                        max(0.5,1.5*Signal.nbchan/Signal.srate);
%
%   'ASR_StepSize' : Step size for processing. The reprojection matrix will be updated every this many
%                    samples and a blended matrix is used for the in-between samples. If empty this will
%                    be set the WindowLength/2 in samples. Default: []
%
%   'MaxDimensions' : Maximum dimensionality to reconstruct. Up to this many dimensions (or up to this 
%                     fraction of dimensions) can be reconstructed for a given data segment. This is
%                     since the lower eigenvalues are usually not estimated very well. Default: 2/3.
%
%   'ReferenceWindowLength' : Granularity at which EEG time windows are extracted
%                             for calibration purposes, in seconds. Default: 1.
%
%   'UseGPU' : Whether to run on the GPU. This makes sense for offline processing if you have a a card with
%              enough memory and good double-precision performance (e.g., NVIDIA GTX Titan or K20). 
%              Note that for this to work you need to a) have the Parallel Computing toolbox and b) remove 
%              the dummy gather.m file from the path. Default: false
%
%   'availableRAM_GB' : Specify the amount of available RAM in GB. Default:[] (empty). It is mandatory to enter
%                       a fixed value here to stabilize the output length from the subsequent window rejection.
%                       It is because ASR's behavior and cleaning result depends on the available RAM size.
%                       By default, the amount of available RAM is obtained by getFreePhysicalMemorySize function
%                       in hlp_memfree(). However, it returns smaller amount than it appears on OS's System Monitor.
%                       It's because definition of 'free' and 'available' are different due to reserving and caching
%                       etc, and they are always 'available' >> 'free'. In our test, >100 iteration of the same
%                       process showed +/- 15% of RAM availability (it depends on background processes etc), and
%                       the final data length changed within the range of +/- 2%. By specifying the fixed number
%                       here (in GB), the output length became identical across iterations. (03/26/2019 Makoto.)
%
%
%
% ------------------- Optional inputs for clean_windows() -----------------
%
%   The following are detail parameters that usually do not have to be tuned. If you can't get
%   the function to do what you want, you might consider adapting these to your data.
%
%   'CleanWindowsWindowLength'  : Window length that is used to check the data for artifact content. This is 
%                                 ideally as long as the expected time scale of the artifacts but not shorter 
%                                 than half a cycle of the high-pass filter that was used. Default: 1.
%
%   'WindowOverlap' : Window overlap fraction. The fraction of two successive windows that overlaps.
%                     Higher overlap ensures that fewer artifact portions are going to be missed (but
%                     is slower). (default: 0.66)
% 
%   'MaxDropoutFraction' : Maximum fraction that can have dropouts. This is the maximum fraction of
%                          time windows that may have arbitrarily low amplitude (e.g., due to the
%                          sensors being unplugged). (default: 0.1)
%
%   'MinCleanFraction' : Minimum fraction that needs to be clean. This is the minimum fraction of time
%                        windows that need to contain essentially uncontaminated EEG. (default: 0.25)
%
%   
%   The following are expert-level parameters that you should not tune unless you fully understand
%   how the method works.
%
%   'TruncateQuantile' : Truncated Gaussian quantile. Quantile range [upper,lower] of the truncated
%                        Gaussian distribution that shall be fit to the EEG contents. (default: [0.022 0.6])
%
%   'CleanWindowsStepSizes' : Grid search stepping. Step size of the grid search, in quantiles; separately for
%                             [lower,upper] edge of the truncated Gaussian. The lower edge has finer stepping
%                             because the clean data density is assumed to be lower there, so small changes in
%                             quantile amount to large changes in data space. (default: [0.01 0.01])
%
%   'ShapeRange' : Shape parameter range. Search range for the shape parameter of the generalized
%                  Gaussian distribution used to fit clean EEG. (default: 1.7:0.15:3.5)
%
%

% Author: Makoto Miyakoshi and Christian Kothe, SCCN,INC,UCSD
%
% History:
% 03/26/2019 Makoto and Chiyuan. Supported 'availableRAM_GB'. GUI switched to GUIDE-made. Options for ASR supported.
% 05/13/2014 ver 1.2 by Christian. Added better channel removal function (uses locations if available).
% 07/16/2013 ver 1.1 by Makoto and Christian. Minor update for help and default values.
% 06/26/2013 ver 1.0 by Makoto. Created.

% Copyright (C) 2013, Makoto Miyakoshi and Christian Kothe, SCCN,INC,UCSD
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function cleanEEG = clean_rawdata(EEG, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst, arg_window, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Disable functions if requested. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ischar(arg_flatline)
    if strcmp(arg_flatline, 'off')
        disp('flatchan rej disabled.');
    else
        error('Invalid input detected.')
    end
elseif arg_flatline == -1
    arg_flatline = 'off';
    disp('flatchan rej disabled.'); 
end
        
if ischar(arg_highpass)
    if strcmp(arg_highpass, 'off')
        disp('highpass disabled.');
    else
        error('Invalid input detected.')
    end
elseif arg_highpass == -1
    arg_highpass = 'off';
    disp('highpass disabled.'); 
end

if ischar(arg_channel)
    if strcmp(arg_channel, 'off')
        disp('badchan rej disabled.');
    else
        error('Invalid input detected.')
    end
elseif arg_channel == -1
    arg_channel = 'off';
    disp('badchan rej disabled.'); 
end

if ischar(arg_noisy)
    if strcmp(arg_noisy, 'off')
        disp('nose-based rej disabled.');
    else
        error('Invalid input detected.')
    end
elseif arg_noisy == -1
    arg_noisy = 'off';
    disp('noise-based rej disabled.'); 
end

if ischar(arg_burst)
    if strcmp(arg_burst, 'off')
        disp('burst clean disabled.');
    else
        error('Invalid input detected.')
    end
elseif arg_burst == -1
    arg_burst = 'off';
    disp('burst clean disabled.'); 
end

if ischar(arg_window)
    if strcmp(arg_window, 'off')
        disp('bad window rej disabled.'  );
    else
        error('Invalid input detected.')
    end
elseif arg_window == -1
    arg_window = 'off';
    disp('bad window rej disabled.'); 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build the main input cells. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mainInputCells{1 } = 'FlatlineCriterion';
mainInputCells{2 } = arg_flatline;
mainInputCells{3 } = 'Highpass';
mainInputCells{4 } = arg_highpass;
mainInputCells{5 } = 'ChannelCriterion';
mainInputCells{6 } = arg_channel;
mainInputCells{7 } = 'LineNoiseCriterion';
mainInputCells{8 } = arg_noisy;
mainInputCells{9 } = 'BurstCriterion';
mainInputCells{10} = arg_burst;
mainInputCells{11} = 'WindowCriterion';
mainInputCells{12} = arg_window;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Combine the main and optional inputs. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
combinedInput = [mainInputCells varargin{:}];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Christian's main function. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Welcome to the labyrinth of Trismegistus who wrote BCILAB, LSL, and SNAP.
cleanEEG = clean_artifacts(EEG, combinedInput);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare the main parameter log. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ischar(arg_flatline)
    if strcmp(arg_flatline, 'off')
        arg_flatline_log = '''off''';
    end
else
   arg_flatline_log = num2str(arg_flatline);
end
      
if isnumeric(arg_highpass)
    if     length(arg_highpass) == 1 % -1
        arg_highpass_log = num2str(arg_highpass);
    elseif length(arg_highpass)==2 % [low high]
        arg_highpass_log = ['[' num2str(arg_highpass(1)) ' ' num2str(arg_highpass(2)) ']'];
    end
else % off
    arg_highpass_log = '''off''';
end

if ischar(arg_channel)
    if strcmp(arg_channel, 'off')
        arg_channel_log = '''off''';
    end
else
   arg_channel_log = num2str(arg_channel);
end

if ischar(arg_noisy)
    if strcmp(arg_noisy, 'off')
        arg_noisy_log = '''off''';
    end
else
   arg_noisy_log = num2str(arg_noisy);
end

if ischar(arg_burst)
    if strcmp(arg_burst, 'off')
        arg_burst_log = '''off''';
    end
else
   arg_burst_log = num2str(arg_burst);
end

if ischar(arg_window)
    if strcmp(arg_window, 'off')
        arg_window_log = '''off''';
    end
else
   arg_window_log = num2str(arg_window);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare optional parameter log. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checking varargin structure. https://www.mathworks.com/matlabcentral/answers/160576-passing-arguments-to-varargin
if length(varargin) == 1
    vararginContents = varargin{:};
else
    vararginContents = varargin;
end

optionalInput_log = '';
for itemIdx = 1:length(vararginContents)
    if itemIdx == 1
        optionalInput_log = [char(39) vararginContents{itemIdx} char(39)];
    else
        if mod(itemIdx,2)==0
            
            % Obtain the numeric item.
            currentItem = vararginContents{itemIdx};
            
            % If multiple entries of numeric item, separatation with one space and square bracket are added.
            if length(currentItem)>1
                for currentItemIdx = 1:length(currentItem)
                    if currentItemIdx == 1
                        multipleInputString = ['[' num2str(currentItem(currentItemIdx))];
                    else
                        multipleInputString = [multipleInputString ' ' num2str(currentItem(currentItemIdx))];
                    end
                end
                multipleInputString = [multipleInputString ']'];
                currentItem = multipleInputString;
            else
                currentItem = num2str(currentItem);
            end
        else
            currentItem = [char(39) vararginContents{itemIdx} char(39)];
        end
        optionalInput_log = [optionalInput_log ', ' currentItem];
    end
end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate the clean_rawdata() log. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clean_rawdataLog = sprintf('EEG = clean_rawdata(EEG, %s, %s, %s, %s, %s, %s, %s);',...
                                                    arg_flatline_log,...
                                                    arg_highpass_log,...
                                                    arg_channel_log,...
                                                    arg_noisy_log,...
                                                    arg_burst_log,...
                                                    arg_window_log,...
                                                    optionalInput_log);
% Update EEG.
cleanEEG.etc.clean_rawdata_log = clean_rawdataLog;