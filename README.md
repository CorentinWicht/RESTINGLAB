# RESTINGLAB
Developed at the Laboratory for Neurorehabilitation Science, University of Fribourg, CH.
Author: Corentin Wicht (corentin.wicht@unifr.ch); Contributor : Christian Mancini (christian.mancini@unifr.ch)

---------------------------------------------------------------------------------------------------------------
RESTINGLAB is an open-source EEGLAB-based (Delorme and Makeig, 2004) standalone software for (semi)automated resting-state BioSemi EEG data pre-processing.

Currently the software is in beta version which means it may still contains errors. 

Preprocessing dependencies:
•	i) CleanLine (sinusoidal, line noise; (see Mullen, 2012) 
•	ii) Artifact Subspace Reconstruction (ASR: non-stationary, high-amplitude bursts; see Mullen et al., 2015; Chang et al., 2018) 
•	iii) BLINKER (detection of eye blinks; see Kleifges et al., 2017). 
•	Independent Components Analysis (ICA) based on the following algorithms: AMICA for the decomposition (Palmer et al., 2008) and ICLabel for the automated classification of artifacts-containing independent components (Pion-Tonachini et al., 2019)
•	Channel(s) interpolation using multiquadric interpolation relying on radial basis functions (see Jäger et al., 2016; Janin, 2018; Buhmann and Jäger, 2019) 

Analyses dependencies:
•	MST toolbox: EEG microstates analyses (Poulsen et al., 2018)
•	Measure Projection Toolbox: probabilistic EEG source localisation (Bigdely-Shamlo et al., 2013)
•	Ept_TFCE and Factorial Mass Univariate ERP Toolbox: Mass univariate/multivariate analyses (Groppe et al., 2011; Mensen and Khatami, 2013)




UPDATES
---------------------------------------------------------------------------------------------------------------

•	v0.62
- Fixed error when importing parameters in SubjectsGUI.mlapp trough StudyGUI.mlapp
- Fixed error when loading parameters in PreprocessingGUI.mlapp (seconds for asleep epochs rejection)

•	v0.63 (TO DO)
- Implement AreasList with channels labels instead of numbers
- Try the 0.5 frequency resolution, also implement a method to have any frequency resolution
- Check the correctedness of the source neighbouring matrix
- Find a solution for statistics on frequency bins (3rd figures)
- Test ASR-interpolation for BLINKER
- Replace ASR 1.0 with ASR 2.0 to fix the problem of randomness?
- Include in GUI possibility to adjust the IC rejection criteria for MPT?
- Finish the MPT source localisation script

