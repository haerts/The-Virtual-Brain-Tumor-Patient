# The-Virtual-Brain-Tumor-Patient
Code used for postprocessing analyses in manuscript "Modeling brain dynamics in brain tumor patients using The Virtual Brain" (Aerts et al. 2018; eNeuro, 5 (3) ENEURO.0083-18.2018; DOI: https://doi.org/10.1523/ENEURO.0083-18.2018).

(1) Run TVBii code (for TVBii code see Schirner et al. 2018 eLife, https://elifesciences.org/articles/28927, https://github.com/BrainModes/The-Hybrid-Virtual-Brain)
> TVBii_initiate.sh
> TVBii_run.sh

(2) Perform parameter space exploration to identify optimal model parameters
> PSE.m

(3) Sanity checks model parameters:
> PSE_SanChecks.m
> SanChecks_Ji.R: evaluate relationship between Ji, in-strength and ROI size
> ModelFit_SanChecks.R: evaluate model fit

(4) Postprocessing:
> GTA_SC.m: compute structural network topology measures
> descriptives.R: descriptive analyses
> regression_models.R: assess associations between fitted model parameters, structural network topology measures and cognitive performance measures
