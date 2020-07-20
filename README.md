# RESTINGLAB

Developed at the Laboratory for Neurorehabilitation Science, University of Fribourg, CH.
Author: Corentin Wicht (corentin.wicht@unifr.ch); Contributor : 

---------------------------------------------------------------------------------------------------------------
RESTINGLAB is an open-source EEGLAB-based (Delorme and Makeig, 2004) standalone software for (semi)automated resting-state BioSemi EEG data pre-processing.

Currently the software is in beta version which means it may still contains errors. 

## Dependencies
| PLUGINS | Description |
| ------ | ------ |
| [EEGLAB v14.1.2b](https://github.com/sccn/eeglab) | Importing the .set EEG files | 
| [NMD v.2.00](http://www.physics.lancs.ac.uk/research/nbmphysics/diats/tfr/) | Computing the Morlet wavelet-based time-frequency decomposition |
| [FMUT v.0.5.1](https://github.com/ericcfields/FMUT) | Computation of permutation-based statistics |
| [ept_TFCE](https://github.com/Mensen/ept_TFCE-matlab) | Computation of permutation-based statistics and TFCE correction |

Isolated functions:
* [Timerwaitbar v1.02](https://ch.mathworks.com/matlabcentral/fileexchange/55985-timer-waitbar) (upgraded)
* [bluewhitered v1.00](https://ch.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered)
* [limo_FDR from the LIMO Toolbox](https://github.com/LIMO-EEG-Toolbox/limo_tools)
* [natsort v2.10](https://ch.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort)
* [saveeph](https://sites.google.com/site/cartoolcommunity/files)
* [textprogressbar v1.00](https://ch.mathworks.com/matlabcentral/fileexchange/28067-text-progress-bar)
* [topoplotIndie](https://www.mikexcohen.com/)
* [vis_artifacts](https://github.com/sccn/clean_rawdata/blob/master/vis_artifacts.m)
* [vline v1.00](https://ch.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline)

The dependencies are already included in the [functions/Dependencies](functions/Dependencies) folder.

## Authors
[**Corentin Wicht**](https://www.researchgate.net/profile/Wicht_Corentin)\
*SNSF Doc.CH PhD student*\
*corentin.wicht@unifr.ch, corentinw.lcns@gmail.com*\
*[Laboratory for Neurorehabilitation Science](https://www3.unifr.ch/med/spierer/en/)*\
*University of Fribourg, Switzerland*

[**Christian Mancini**](https://www.researchgate.net/profile/Christian_Mancini)\
*SNSF Doc.CH PhD student*\
*corentin.wicht@unifr.ch, corentinw.lcns@gmail.com*\
*[Laboratory for Neurorehabilitation Science](https://www3.unifr.ch/med/spierer/en/)*\
*University of Fribourg, Switzerland*
 (christian.mancini@unifr.ch)

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

See the [LICENSE.md](LICENSE.md) file for details

## Acknowledgements
PD Dr. Lucas Spierer, President of the Section of Medicine (Faculty of Science and Medicine, University of Fribourg, Switzerland), Professor of
Neurology and Director of the [Laboratory for Neurorehabilitation Science (LNS, UNIFR)](https://www3.unifr.ch/med/spierer/en/) provided substantial support and advices regarding theoretical conceptualization as well as access to the workplace and the infrastructure required to successfully complete the project. Additionally, [Hugo Najberg](https://github.com/HugoNjb) and [Dr. Michael Mouthon](https://www3.unifr.ch/med/fr/section/personnel/all/people/3229/6a825) provided valuable advices regarding programming issues in MATLAB and technical support.

## Fundings
This research was supported by [Swiss National Science Foundation](http://www.snf.ch/fr/Pages/default.aspx) grants:
* [#P0LAP1_181689](http://p3.snf.ch/project-181689) to Corentin Wicht
* [#320030_175469](http://p3.snf.ch/project-175469) to PD Dr. Lucas Spierer

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
