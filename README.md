# RESTINGLAB
Developed at the Laboratory for Neurorehabilitation Science, University of Fribourg, CH.
Author: Corentin Wicht (corentin.wicht@unifr.ch); Contributor : Christian Mancini (christian.mancini@unifr.ch)

---------------------------------------------------------------------------------------------------------------
RESTINGLAB is an open-source EEGLAB-based (Delorme and Makeig, 2004) standalone software for (semi)automated resting-state BioSemi EEG data pre-processing.

Currently the software is in beta version which means it may still contains errors. 

Preprocessing dependencies:
-	i) CleanLine (sinusoidal, line noise; (see Mullen, 2012) 
-	ii) Artifact Subspace Reconstruction (ASR: non-stationary, high-amplitude bursts; see Mullen et al., 2015; Chang et al., 2018) 
-	iii) BLINKER (detection of eye blinks; see Kleifges et al., 2017). 
-	Independent Components Analysis (ICA) based on the following algorithms: AMICA for the decomposition (Palmer et al., 2008) and ICLabel for the automated classification of artifacts-containing independent components (Pion-Tonachini et al., 2019)
-	Channel(s) interpolation using multiquadric interpolation relying on radial basis functions (see Jäger et al., 2016; Janin, 2018; Buhmann and Jäger, 2019) 

Analyses dependencies:
-	MST toolbox: EEG microstates analyses (Poulsen et al., 2018)
-	Measure Projection Toolbox: probabilistic EEG source localisation (Bigdely-Shamlo et al., 2013)
-	Ept_TFCE and Factorial Mass Univariate ERP Toolbox: Mass univariate/multivariate analyses (Groppe et al., 2011; Mensen and Khatami, 2013)




UPDATES
---------------------------------------------------------------------------------------------------------------

•	v0.62
- Fixed error when importing parameters in SubjectsGUI.mlapp trough StudyGUI.mlapp
- Fixed error when loading parameters in PreprocessingGUI.mlapp (seconds for asleep epochs rejection)
- Implemented AreasList with channels labels instead of numbers
- All plugins are now up-to-date and have been tested (11.12.2019)
- Implemented frequency resolution smaller than integer (e.g. 0.5)
- Solved the problem of the 0 frequency bin when using pwelch (i.e. DC Offset : https://github.com/sccn/eeglab/issues/101)
- Separated the channels rejection from interpolation and re-referencing, following the pipeline from: https://www.frontiersin.org/articles/10.3389/fnins.2018.00097/full. Now bad channels interpolation and robust referencing are performed at the end of the script (after optional ICA). 

•	v0.63 (TO DO)
- Use triggers to restrain data length (!)
- Check the source neighbouring matrix
- Find a solution for statistics on frequency bins (3rd figures)
- Include in GUI possibility to adjust the IC rejection criteria for MPT
- Finish the MPT source localisation script
- Implement a try-catch for excel template generation! The software crashes if the name of sheets is changed and files already exist
- Error when opening 2.EEG parameters (when parameters are already loaded)

•	Problems to solve in the future:
- Find a way to implement ASR-interpolation for BLINKER
- Implement the method for 3D electrode localisation 
- Find a way to run the sleep detection algorithm

