#! /bin/sh
###########################################################
### MARKER DISCOVER AND CELLULAR LOCATION ANNOTATION
###########################################################
# STEP1 
Rscript tumor_specific_markers_V1.0.R -r /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/rds_list.txt -o /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/ -c Plasma  -s /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Data/Cell_Surface_Protein_Atlas_S2_File.txt -p /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Data/Human_Protein_Atlas_subcellular_location.txt 

###########################################################
### GTEX TISSUE SPECIFICITY TEST
###########################################################
OUT_DIR=/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/
mkdir -p ${OUT_DIR}
####################################################
### step0 parse the GTEX table to get meta.data
####################################################
perl /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Results/V1/step0.get_surface_genes_metadata.pl > ${OUT_DIR}/GTEX_metadata_subset_tumor_cells_DEG.txt
###########################################################################################
### step1 parse the GTEX table to get expressions only for the marker genes
### cd to the folder containing DE_genes_filtered_surface_3DB.txt
### please do not change the file name of "GTEX_expression_data_subset_tumor_cells_DEG.txt"
#########################################################################################
SCRIPT_DIR=/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/GTEX
GTEX_DIR=/diskmnt/Datasets/GTEX_tpm
CELL_TYPE=Plasma
perl ${SCRIPT_DIR}/step1.get_surface_genes_expression_data.pl ${OUT_DIR}/${CELL_TYPE}_specific_DEG_with_surface_annotations_from_3DB_gene_list.txt ${GTEX_DIR} > ${OUT_DIR}/GTEX_expression_data_subset_tumor_cells_DEG.txt
grep captured ${OUT_DIR}/GTEX_expression_data_subset_tumor_cells_DEG.txt > ${OUT_DIR}/not_captured_GTEX_DEG.txt
sed -e s/#//g -i ${OUT_DIR}/not_captured_GTEX_DEG.txt
################################################################################################
### step2 get statistically tissue-specific markers
### step2.analyze_expression_change.pl will find tissue-specific markers: for the query gene, which tissue has the significant high expression. It includes 1) removes the outlier tissues: Statistics::Outliers 2) do the t-test for the tissues with similar high expression 3) multi-test correction to get fdr
################################################################################################
#LIB_DIR=/diskmnt/results/lyao/denali_lyao/Melanoma/ HOW TO ADD THIS AS A VARIABLE
/diskmnt/results/lyao/denali_lyao/LabCode_Tools/anaconda3/bin/perl ${SCRIPT_DIR}/step2.analyze_expression_change.pl ${OUT_DIR} > ${OUT_DIR}/GTEX_tissue_specific_DEG.txt

###########################################################
###HPA RNA test
###########################################################
/diskmnt/results/lyao/denali_lyao/LabCode_Tools/anaconda3/bin/perl /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/HPA_RNA/find_outlier_expressions_t_test.pl ${OUT_DIR}/${CELL_TYPE}_specific_DEG_with_surface_annotations_from_3DB.txt /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Data/rna_tissue_consensus.tsv > ${OUT_DIR}/HPA_RNA_tissue_specific_DEG.txt

###########################################################
###HPA Protein test
###########################################################
python3 HPA_Protein_filtering.py -p /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Data/HPA_normal_tissue.tsv -o ${OUT_DIR} -t /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Data/HPA_Tissue_type_matching.txt -g ${OUT_DIR}/Plasma_specific_DEG_with_surface_annotations_from_3DB.txt

############################################################
###Adding tissue specificty prediction from GTEX and HPA
#############################################################
Rscript adding_tissue_specificity_annotations.R -o /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test -c Plasma -g lymphocytes -r bone_marrow,thymus,spleen,appendix,tonsil,lymph_node -p BM_and_lymphoid_tissues 

