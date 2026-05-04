# GSE329659_scripts
Scripts to analyze the data in GEO accession number GSE329659

R_preprocess_samples.R - script to load SpaceRanger output as Seurat objects and process each of the four samples to the point of having UMAPs.

R_boger_logNorm_all_samples.R - original script to merge the four samples, integrate with Harmony, transfer labels from the reference dataset, and generate the first round of figures. 

R_revisions_placenta.R - script generated to perform the revisions requested from reviewer comments. This includes running CCA Integration (replaces Harmony in generating figures) and correcting annotations of Macrophages that express CD4 as T-cells. All figures generated in the manuscript were generated from this script. RDS file in GEO repository was generated in this script.
