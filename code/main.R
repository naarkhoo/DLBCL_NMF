rm(list = ls())

# conver the matrix data to gct format. The matrix should be genes x perturbation
source("/home/superiois/Downloads/projectx3/DLBCL_NMF/code/DISSECTOR_lib.v3.R")
data1 = read.csv("/home/superiois/Downloads/projectx3/DLBCL_NMF/data/DLBCL_OE_L1000_NMF_TA.OE002_A375_72H.csv", sep = ",")
data1 = t(data1)
data1 = as.data.frame(data1)
colnames(data1) = as.matrix(data1[1,])
data1 = data1[-1,]
CNMF.write.gct.2(data1, descs = ",", filename = "/home/superiois/Downloads/projectx3/DLBCL_NMF/data/DLBCL_OE_L1000_NMF_TA.OE002_A375_72H.gct") 


source("/home/superiois/Downloads/projectx3/DLBCL_NMF/code/DISSECTOR_lib.v3.R")

DISSECTOR_create_components.v1(            #  Project an input dataset into components using NMF  (A ~ W x H)
	input_dataset          = "/home/superiois/Downloads/projectx3/DLBCL_NMF/data/DLBCL_OE_L1000_NMF_TA.OE002_A375_72H.gct",
	input_normalization    = "rank",        # Normalization for the input dataset: "rank"
	number_of_comp         = 7,             # Number of components to use in the matrix decomposition
	method                 = "NMF",         # Method for matrix factorization: NMF or NMF_offset (IMF under construction)                                  
	gene_subset            = "all-genes",   # Universe of genes to consider for matrix decomposition: "gene-sets", "all-genes"
	gene_sets_files        = NULL,          # If gene_subset = "gene-sets" GMT file with gene sets
	gene_sets              = NULL,          # If gene_subset = "gene-sets" then name of the specific gene set(s) in gene_sets_file to use
	output_plots           = "/home/superiois/Downloads/projectx3/DLBCL_NMF/plots/compound_A375_7.PLOTS.pdf", # Output plots file (PDF, W, H and A)
	output_W_dataset       = "/home/superiois/Downloads/projectx3/DLBCL_NMF/outputs/compound_A375_7.W.gct",     # Output W matrix (GCT file)
	output_H_dataset       = "/home/superiois/Downloads/projectx3/DLBCL_NMF/outputs/compound_A375_7.H.gct",     # Output H matrix (GCT file)
	output_H_w_dataset     = "/home/superiois/Downloads/projectx3/DLBCL_NMF/outputs/compound_A375_7.H_w.gct")   # Output W-derived H matrix (GCT file)


source("/home/superiois/Downloads/projectx3/DLBCL_NMF/code/DISSECTOR_lib.v3.R")
DISSECTOR_make_heatmap_of_matrix.v1(
	input_dataset            = "/home/superiois/Downloads/projectx3/DLBCL_NMF/outputs/compound_A375_7.H_w.gct", # Input dataset (GCT).
	annot_file               = c("/home/superiois/Downloads/projectx3/DLBCL_NMF/data/compound_classes_A375_7.v3.csv", "Sample", "Type", T),
	transpose_data           = F,           # Transpose input matrix
	append_annot             = T,           # Append annotation to column names
	sorting_method           = "HC",        # Sorting method for cols inside phenotype: MDS (Multi_dimensinal Scaling) or HC (Hiererachical Clustering)  
	output_plot_landscape    = "/home/superiois/Downloads/projectx3/DLBCL_NMF/plots/A375_DLBCL_n978x323.H_w_LPLOT.pdf",  # Output (PDF) file
	output_plot_portrait     = "/home/superiois/Downloads/projectx3/DLBCL_NMF/plots/A375_DLBCL_n978x323.H_w_PPLOT.pdf")  # Output (PDF) file




source("/home/superiois/Downloads/projectx3/DLBCL_NMF/code/DISSECTOR_lib.v3.R")
DISSECTOR_project_dataset.v1(     # Project a dataset in the space defined by a W matrix
	input_dataset           = "/home/superiois/Downloads/projectx3/DLBCL_NMF/data/rnaseq.v3.gct",  # Input dataset (GCT) (RNA seq should be the input)
	input_normalization     = "rank",     # Normalization for the input dataset: "rank"
	normalize_after_match   = T,        # Normalize input dataset after matching with rows of W
	input_W_dataset         = "/home/superiois/Downloads/projectx3/DLBCL_NMF/outputs/compound_A375_7.W.gct",  # my W Input W matrix (GCT)
	W_normalization         = "none",         # Normalization for W                                            
	output_H_dataset        = "/home/superiois/Downloads/projectx3/DLBCL_NMF/outputs/DLBCL_OE_L1000_NMF_TA.OE002_PC3_72H_H.gct")



source("/home/superiois/Downloads/projectx3/DLBCL_NMF/code/DISSECTOR_lib.v3.R")
DISSECTOR_make_heatmap_of_matrix.v1(
	input_dataset            = "/home/superiois/Downloads/projectx3/DLBCL_NMF/outputs/DLBCL_OE_L1000_NMF_TA.OE002_PC3_72H_H.gct", # Input dataset (GCT).
	annot_file               = c("/home/superiois/Downloads/projectx3/DLBCL_NMF/data/CCLE_Kim_Tamayo_96AffyCEL_FAZES_Oct12.v2.txt","CCLE_name", "Site.Primary", F),
	transpose_data           = F,           # Transpose input matrix
	append_annot             = T,           # Append annotation to column names
	sorting_method           = "HC",        # Sorting method for cols inside phenotype: MDS (Multi_dimensinal Scaling) or HC (Hiererachical Clustering)  
	output_plot_landscape    = "/home/superiois/Downloads/projectx3/DLBCL_NMF/plots/ccle_DLBCL_H_proj_LPLOT.v1.pdf",  # Output (PDF) file
	output_plot_portrait     = "/home/superiois/Downloads/projectx3/DLBCL_NMF/plots/ccle_DLBCL_H_proj_PPLOT.v1.pdf")  # Output (PDF) file



source("/home/superiois/Downloads/projectx3/DLBCL_NMF/code/DISSECTOR_lib.v3.R")
DISSECTOR_extract_gene_sets.v1(            #  Extract/Define gene sets from a W matrix
	input_dataset      	 = "/home/superiois/Downloads/projectx3/DLBCL_NMF/outputs/compound_A375_7.W.gct",          # Input W matrix (GCT file)
	min.genes.per.comp 	 = 5,          # Minimun number of genes per component gene set
	max.genes.per.comp 	 = 15,         # Maximun number of genes per component gene set
	output_gene_sets_file 	 = "/home/superiois/Downloads/projectx3/DLBCL_NMF/outputs/DLBCL_classes_W_gene_sets.v1.gmt", # Output GMT file with gene sets
	output_plots          	 = "/home/superiois/Downloads/projectx3/DLBCL_NMF/plots/DLBCL_classes_A375_K7.v1.pdf")     # Output PDF file with W and H plots




### under process
source("/home/superiois/Downloads/projectx3/DLBCL_NMF/code/DISSECTOR_lib.v3.R")
source("/home/superiois/Downloads/projectx3/DLBCL/Dissector/FS.library.v8.6.R")

   DISSECTOR_annotate_components.v1(
      input_dataset         = "/home/superiois/Downloads/projectx3/DLBCL/Dissector/A375_7/A375_2/data/DLBCL_OE_L1000_NMF_TA.OE002_PC3_72H_H.gct",
      directory             = "/home/superiois/Downloads/projectx3/DLBCL/Dissector/A375_7/A375_2/annotation",                            
      identifier            = "Comp_annot.v1",
      feature.type.files    =  list("ACHILLES"     = "/home/superiois/Downloads/projectx3/DLBCL/Dissector/Achilles_v2.4.1.rnai.Gs.gct"),
                                    "RPPA"         = "/home/superiois/Downloads/projectx3/DLBCL/Dissector/RPPA.dat.gct")#,
                                    "MUT_CNA"      = "/home/superiois/Downloads/projectx3/DLBCL/Dissector/RNAseqMUTs_CNA_20130729.gct",  
                                    "EXP_PATHWAYS" = "/home/superiois/Downloads/projectx3/DLBCL/Dissector/CCLE_MSigDB_plus_oncogenic.PATHWAYS.v2.gct",
                                    "EXP_GENES"    = "/home/superiois/Downloads/projectx3/DLBCL/Dissector/rnaseq.v3.gct"),
       feature.directions   = c(0, 1, 1, 1, 1),                                         
       n.markers            = 15,                         # Number of top hits shown in the heatmaps
       n.perm               = 2,                                         
       char.scaling         = 0.575,                   # Character scaling for heatmaps
       locs.table.file      = "/home/superiois/Downloads/projectx3/DLBCL/Dissector/hgnc_downloads.txt",
       log.table.file       = "/home/superiois/Downloads/projectx3/DLBCL/Dissector/annot_comp.log.txt",                            
       min.thres            = 10,
       character.scaling    = 0.65,
       phen.table           = NULL,
       phenotypes           = NULL)  





