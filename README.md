# CompTox-PFASHalfLife

This repository accompanies the manuscript "A Machine Learning Model to Estimate Toxicokinetic Half-lives of Per- and Polyfluoro-alkyl Substances (PFAS) in Multiple Species" which is currently submitted to a journal for peer-review.

The major steps of the workflow for this study included training dataset assembly, predictor set assembly, model construction, and model application.

All analyses were performed using the freely available R statistical software platform (v4.3.0)

Use CurrentScripts/PFAS_QSAR_MakePredictions_060622.R to make predictions for novel chemicals and species

Dataset assembly is described in detail in the Supplemental Information (S1.1). 

## Supplemental Information

* S1_Dawson et al._ML PFAS_HL_052622.docx
* S2_Dawson et al._ML PFAS_HL_053022.xlsx
* S3_Dawson et al._ML PFAS_HL_060622.zip

## Abstract

Per- and polyfluoroalkyl substances (PFAS) are a diverse group of man-made chemicals that are known to be widely distributed in the environment and are commonly found in body tissues. While understanding the toxicokinetics of PFAS is necessary to determine potential risk to public health, most are uncharacterized. Long toxicokinetic half-lives (t½) for some compounds in conjunction with frequent exposure pose a potential for bioaccumulation. We used an ensemble machine learning method, random forest, to model the existing in vivo measured half-lives (t½) across 4 species (human, monkey, rat, mouse) and 11 PFAS. Chemical descriptors included physico-chemical properties and protein binding coefficients. In addition, two types of surrogates for renal transporter affinity were considered: 1) physiological descriptors of kidney geometry and 2) chemical similarity to endogenous, non-fluorinated chemicals that might themselves be transporter substrates. We developed a 4-Bin classification model (Bin 1: <12 hours; Bin 2:12 hours – 1 week; Bin 3: 1 week - 2 months; Bin 4: > 2 months) of t½ that uses 16 predictors. The classification model had an average accuracy of 86.4% ± 10.5% across species and chemicals.  In contrast, y-randomized training data had an average accuracy of 31.7% ± 12.7%. When applied to USEPA’s largest list of 6603 PFAS, 4221 compounds were estimated to be within domain of the model. For these 4221 chemicals, human t½ was predicted to be distributed such that 52% were classified in Bin 4, 12% were classified in Bin 3, and 35% were classified in Bin 2. For a large class of chemicals with public health and environmental concerns this machine learning model synthesizes the limited available data to allow both tentative extrapolation as well as prioritize targets for additional investigation.

## Authors

* Daniel E. Dawson
* Christopher S. Lau [lau.christopher@epa.gov]
* Prachi Pradeep
* Risa R. Sayre [sayre.risa@epa.gov]
* Richard S. Judson [judson.richard@epa.gov]
* Rogelio Tornero-Velez [tornero-velez.rogelio@epa.gov]
* John F. Wambaugh [wambaugh.john@epa.gov] 

##

To replicate the results, clone this archive and then use R to run the scripts
in /CurrentScripts in order:
* 1_PFAS_Dataset_building.R
* 2_PFAS_ML_ModelBuilding.R
* 3_PFAS_ML_MakePredictions.R
* 4_PFAS_ML_ModelBuilding_NoLogD.R
* 5_PFAS_ML_MakePredictions_NoLogD.R
* 6_PFAS_ML_Figs_Restricted_by_AD.R
* 7_ClassyFire_of_DomainChemicals.R
* 8_Make_Fig3.R
* 9_ReportedResults.R
