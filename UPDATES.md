20.07.2020

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



• TO DO
Try using Fieldtrip MonteCarlo Permutations with Max Cluster correction for statistics on frequency bins (3rd figures)
- Problem with std.stat = https://github.com/sccn/eeglab/issues/184
Check the source neighbouring matrix
Include in GUI possibility to adjust the IC rejection criteria for MPT
Finish the MPT source localisation script
Implement a try-catch for excel template generation! The software crashes if the name of sheets is changed and files already exist


• Problems to solve in the future:

Find a way to implement ASR-interpolation for BLINKER
Find a way to run the sleep detection algorithm