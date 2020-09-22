# MiR_Combo_Targeting
Companion code for the paper: "Network potential identifies therapeutic miRNA cocktails in Ewing sarcoma"
You can find the pre-print for this analysis here.https://www.biorxiv.org/content/10.1101/854695v1
If you have any questions regarding how to run this pipeline on your device/ use it for your research, 
feel free to shoot me an email at davis.weaver@case.edu. 

## Here is the general flow: 
1. Main_EWS_miR (by calling the functions in NP_pipeline_MiR and NP_Funcs_MiR) converts the processed RNA sequencing data and the biogrid PPI into cell signaling networks and does in-silico repression experiments to generate a list of protein targets. This step made use of ~20 cores of an HPC. Attempting to do it on a laptop or desktop would likely take a period of days to weeks. 
2. miR_Analysis_Ewing (by calling the functions in miR_Analysis_Functions) processes the results of Main_EWS_miR and evaluates miR cocktails to identify optimal 3 miR combinations. 
3. miR_Figs_Ewing produces all the figures and tables from the paper using the output of both miR_Analysis_Ewing and Main_EWS_miR.

## Dependencies:
Package dependencies are at the top of each script. 
Each script also requires data files to run. 

1. Main_EWS_miR.R
    * rld_counts.csv
    * BIOGRID-ORGANISM-Homo_sapiens-3.5.171.tab2.txt
2. miR_Analysis_Ewing.R
    * EwingsDatasetFinal.csv
    * all_miRNA_targets.rda
    * Census_allMon_2019.csv
    * Housekeeping_GenesHuman.csv
3. miR_Figs_Ewing
    * Protein_RNA_correlation.csv
 
Create a folder in the folder where you store the code called "data_files". Download these 7 dependencies do that folder and then it should (hopefully) all go. 
