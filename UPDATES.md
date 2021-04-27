25.09.2020

• v0.62

Fixed error when importing parameters in SubjectsGUI.mlapp trough StudyGUI.mlapp
Fixed error when loading parameters in PreprocessingGUI.mlapp (seconds for asleep epochs rejection)
Implemented AreasList with channels labels instead of numbers
All plugins are now up-to-date and have been tested (11.12.2019)
Implemented frequency resolution smaller than integer (e.g. 0.5)
Solved the problem of the 0 frequency bin when using pwelch (i.e. DC Offset : https://github.com/sccn/eeglab/issues/101)
Separated the channels rejection from interpolation and re-referencing, following the pipeline from: https://www.frontiersin.org/articles/10.3389/fnins.2018.00097/full. Now bad channels interpolation and robust referencing are performed at the end of the script (after optional ICA).

• v0.62.1
DipFit.bmp figures are now saved as .fig and not .bmp (useless)
Implemented Fieldtrip MonteCarlo Permutations with Max Cluster correction for statistics on frequency bins (3rd figures)
Created a -dev branch in which I updated all the plugins (i.e. EEGLAB v.2020.0 + updated plugins)
Temporarilly removed MPT toolbox due to incomptability issues with EEGLAB v.2020.0 and study functions:
https://github.com/sccn/eeglab/issues/198
https://github.com/sccn/eeglab/issues/184

• v0.62.2
Cleaned the code in MAIN.m
Corrected a major bug in MAIN.m where the wrong datasets (i.e. the ones where rejected channels were not interpolated and placed back in the EEG.data field) were imported for ICA cleaning and Power Spectrum Computation
Changed the mandatory export datasets to be the "Artifacts Rejected" in the EEGParamGUI.mlapp GUI instead of "Filtered" (see error below).
Updated plugins: DIPFIT v.3.4 & FieldTrip-lite v.20200922

• v0.63
Updated EEGLAB and plugins to latest versions
Replaced MPT toolbox (deprecated) by EEGLAB updated versions (mp_clustering & NIMA)
Fixed an issue related to participant_load variable and use of "contains" when the "pattern" is longer than the "str" argument.

• TO DO
Finish the MPT source localisation script by EEGLAB new implementation 