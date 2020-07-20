# RESTINGLAB

RESTINGLAB is an open-source EEGLAB-based (Delorme and Makeig, 2004) standalone software for (semi)automated resting-state BioSemi EEG data pre-processing.

*Currently the software is in beta version which means it may still contains errors.* 


## Dependencies
| PLUGINS | Description |
| ------ | ------ |
| [EEGLAB v14.1.2b](https://github.com/sccn/eeglab) | Importing the .set EEG files | 
| [FMUT v0.5.1](https://github.com/ericcfields/FMUT) | Computation of permutation-based statistics |
| [ept_TFCE](https://github.com/Mensen/ept_TFCE-matlab) | Computation of permutation-based statistics and TFCE correction |
| [Automatic Human Sleep Stage Scoring Using Deep Neural Networks](https://github.com/alexander-malafeev/feature-based-sleep-scoring) | |
| [LBPA40 atlas v2011-04-28](https://resource.loni.usc.edu/resources/atlases-downloads/) | |
| [CubeHelix v2.0](https://github.com/DrosteEffect/CubeHelix) | Color scheme generator |
| [mColonFolder v1.6.0](https://ch.mathworks.com/matlabcentral/fileexchange/29854-multiple-colon) | Color scheme generator |
| [Gramm](https://github.com/piermorel/gramm v2.0) | Color scheme generator |
| [PrepPipeline v0.55.3](http://vislab.github.io/EEG-Clean-Tools/) | Color scheme generator |

Isolated functions:
[ShadedErrorBar XXX](https://ch.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar)
XXX


| EEGLAB EXTENSIONS | Description |
| ------ | ------ |
| [AMICA v1.5](https://github.com/japalmer29/amica) | The Adaptive Mixture Independent Component Analysis (AMICA) toolbox provides the best IC decomposition | 
| [MPT v1.661](https://sccn.ucsd.edu/wiki/MPT) |  Probabilistic approach to EEG source comparison and multi-subject inference | 
| [BLINKER v1.1.2](http://vislab.github.io/EEG-Blinks/) | BLINKER  is an automated pipeline for detecting eye blinks in EEG and calculating various properties of these blinks | 
| [ASR v2.0](https://github.com/sccn/clean_rawdata) | ASR (automated subspace removal) detects and rejects or removes high-amplitude non-brain ('artifact') activity (produced by eye blinks, muscle activity, sensor motion, etc.) by comparing its structure to that of known artifact-free reference data | 
| [CleanLine v1.04](https://github.com/sccn/cleanline) | This plugin adaptively estimates and removes sinusoidal (e.g. line) noise from your ICA components
or scalp channels using multi-tapering and a Thompson F-statistic | 
| [EEGBrowser v1.0 ](https://github.com/aojeda/EEGBrowser) | Enhanced visualization for continuous EEG recordings | 
| [fitTwoDipoles v0.01](https://link.springer.com/chapter/10.1007%2F978-3-319-32703-7_22) | Routine for automated recommendation of ICs that may be best fit with a position-symmetric dual-dipole model | 
| [ICLabel v1.2.4](https://sccn.ucsd.edu/wiki/ICLabel) | Automatic independent component (IC) classifcation based on the ICLabel project's classifier | 
| [MST v1.0](https://github.com/atpoulsen/Microstate-EEGlab-toolbox) | Tee toolbox includes ability to run microstate analysis on both ERP and spontaneous (e.g. resting state) EEG | 
| [Viewprops v1.5.4](https://sccn.ucsd.edu/wiki/Viewprops) | Enhanced visualization of independent components (ICA) | 

The dependencies are already included in the [functions/Dependencies](functions/Dependencies) folder.

## Authors
[**Corentin Wicht**](https://www.researchgate.net/profile/Wicht_Corentin)\
*SNSF Doc.CH PhD student*\
*corentin.wicht@unifr.ch, corentinw.lcns@gmail.com*\
*[Laboratory for Neurorehabilitation Science](https://www3.unifr.ch/med/spierer/en/)*\
*University of Fribourg, Switzerland*

[**Christian Mancini**](https://www.researchgate.net/profile/Christian_Mancini)\
*Research assistant*\
*christian.mancini@unifr.ch*\
*[Laboratory for Cognitive and Neurological Sciences](https://www3.unifr.ch/med/annoni/en/)*\
*University of Fribourg, Switzerland*

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

See the [LICENSE.md](LICENSE.md) file for details

## Acknowledgements
PD Dr. Lucas Spierer, President of the Section of Medicine (Faculty of Science and Medicine, University of Fribourg, Switzerland), Professor of
Neurology and Director of the [Laboratory for Neurorehabilitation Science (LNS, UNIFR)](https://www3.unifr.ch/med/spierer/en/) provided substantial support and advices regarding theoretical conceptualization as well as access to the workplace and the infrastructure required to successfully complete the project. Additionally, [Hugo Najberg](https://github.com/HugoNjb) and [Dr. Michael Mouthon](https://www3.unifr.ch/med/fr/section/personnel/all/people/3229/6a825) provided valuable advices regarding programming issues in MATLAB and technical support.

--------------
**WHAT ABOUT JMA???**
--------------

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
