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

• TO DO
Test Spectrum Interpolation instead of CleanLine : https://github.com/fieldtrip/fieldtrip/blob/master/preproc/ft_preproc_dftfilter.m (need to test it!; https://www.sciencedirect.com/science/article/abs/pii/S1053811919300266?via%3Dihub)
Check the source neighbouring matrix
Include in GUI possibility to adjust the IC rejection criteria for MPT
Finish the MPT source localisation script

• Problems to solve in the future:

Find a way to implement ASR-interpolation for BLINKER
Find a way to run the sleep detection algorithm