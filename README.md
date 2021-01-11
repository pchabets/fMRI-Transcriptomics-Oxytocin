# fMRI-Transcriptomics-Oxytocin
Differential oxytocin mRNA expression in brain areas affected by intranasal oxytocin.


Scripts:
1)subcortex_cortex_samples.R - extract subcortical and cortical samples
2)inoxt_fmri_gene_expression.R - data preprocessing and differential gene expression analysis
3)inoxt_frmi_gene_correlations - data preprocessing and (differential) correlation analysis
4)male_only_analysis - data preprocessing and differential gene expression analysis for male-only

Files (in 'data' folder):
1)corticalSamplesOnly.csv - output from script 1)
2)SubcorticalSamplesOnly.csv - output from script 1)
3)reannotation_ALL.txt - output from running ReAnnotator pipeline (done May 5th 2020)
4)reannotated_probes.csv - reformatted and selected for single probe to gene reannotations (probes with reannotation to multiple genes that are not synonyms are deleted)
5)EMOTION_OXT-PBO_P.nii - p-statistic map emotional processing (Grace et al. 2018)
6)SOCIAL_OXT-PBO_P.nii - p-statistic map social processing (Grace et al. 2018)
7)OXT-PBO_ALL_P.nii - p-statistic map all-tasks (Grace et al. 2018)
8)MALE_OXT-PBO_P.nii - p-statistic map all-tasks, male-only (Grace et al. 2018)


Steps to run scripts and reproduce results:
1) download all needed scripts and the 'data' folder holding all needed files from this repository and store them in a local folder
2) download all donor microarray data from https://human.brain-map.org/static/download and place in the same local folder in the 'data' directory (that holds all the other needed files). The six donor folders

-normalized_microarray_donor9861
-normalized_microarray_donor10021
-normalized_microarray_donor12876
-normalized_microarray_donor14380
-normalized_microarray_donor15496
-normalized_microarray_donor15697

should be renamed respectively:

-normalized_microarray_donor01
-normalized_microarray_donor02
-normalized_microarray_donor03
-normalized_microarray_donor04
-normalized_microarray_donor05
-normalized_microarray_donor06

3) in the beginning of each script, set the working directory to the 'data' directory in your local folder that contains all the needed files.
4) install all needed packages before running each script.
5) run the scripts.
