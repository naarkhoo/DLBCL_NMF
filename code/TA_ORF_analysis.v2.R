##### Target Accelerator Reference 4 Pathways ORF Project
#
# Create components for group of 113 ORFs corresponding to 4 pathways (WNT, MAPK, NFKB and AKT)

   source("~/CGP2014/TA/TA_PATHWAYS/DISSECTOR_lib.v3.R")   

   # Make components

   DISSECTOR_explore_component_creation.v1(
      input_dataset             = "~/CGP2014/TA/TA_PATHWAYS/TAref4pathways_OE_L1000_NMF_TA.OE003_A375_72H_20131012.gct", # Input A dataset (GCT)
      input_normalization       = "rank",                                      # Normalization for the input dataset: "rank"
      k.min = 2,                        # Range of components: minimum
      k.max = 10,                        # Range of components: maximum
      k.incr = 1,                       # Range of components: increment                                         
      number_of_runs            = 30,              # Number of runs to explore in the matrix decomposition
      method                    = "NMF",          # Method for matrix factorization: NMF or NMF_offset (IMF under construction)
      gene_subset               = "all-genes",                                 # Universe of genes to consider for matrix decomposition: "gene-sets", "all-genes"
      gene_sets_files           = NULL,
      gene_sets                 = NULL,
      normalize_after_selection = F,                                                  # If gene_subset = "gene-sets," normalize after selection the gene subset
      output_plots              = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/TAref4pathways_A375.EXPLORE.PLOTS.v1.pdf")  # Output plots with results

   source("~/CGP2014/TA/TA_PATHWAYS/DISSECTOR_lib.v3.R")

   DISSECTOR_create_components.v1(            #  Project an input dataset into components using NMF  (A ~ W x H)
      input_dataset          = "~/CGP2014/TA/TA_PATHWAYS/TAref4pathways_OE_L1000_NMF_TA.OE003_A375_72H_20131012.gct", # Input A dataset (GCT)
      input_normalization    = "rank",        # Normalization for the input dataset: "rank"
      number_of_comp         = 7,             # Number of components to use in the matrix decomposition
      method                 = "NMF",         # Method for matrix factorization: NMF or NMF_offset (IMF under construction)                                  
      gene_subset            = "all-genes",   # Universe of genes to consider for matrix decomposition: "gene-sets", "all-genes"
      gene_sets_files        = NULL,         # If gene_subset = "gene-sets" GMT file with gene sets
      gene_sets              = NULL,          # If gene_subset = "gene-sets" then name of the specific gene set(s) in gene_sets_file to use
      output_plots           = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/TAref4pathways_A375.PLOTS.v1.pdf", # Output plots file (PDF, W, H and A)
      output_W_dataset       = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/TAref4pathways_A375.W.v1.gct",     # Output W matrix (GCT file)
      output_H_dataset       = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/TAref4pathways_A375.H.v1.gct",     # Output H matrix (GCT file)
      output_H_w_dataset     = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/TAref4pathways_A375.H_w.v1.gct")   # Output W-derived H matrix (GCT file)

# Plot the H matrix

   source("~/CGP2014/TA/TA_PATHWAYS/DISSECTOR_lib.v3.R")

   DISSECTOR_make_heatmap_of_matrix.v1(
      input_dataset            = "~/CGP2014/TA/TA_PATHWAYS/Analysis_V1/TAref4pathways_A375.H_w.v1.gct", # Input dataset (GCT).
      annot_file               = c("~/CGP2014/TA/TA_PATHWAYS/TARef_v2_Yash_20130507.txt", "Gene.Symbol", "Pathway.Cellular.Process", T),  
      transpose_data           = F,           # Transpose input matrix
      append_annot             = T,           # Append annotation to column names
      sorting_method           = "HC",        # Sorting method for cols inside phenotype: MDS (Multi_dimensinal Scaling) or HC (Hiererachical Clustering)  
      output_plot_landscape    = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/TAref4pathways_A375.H_w.v1_LPLOT.v1.pdf",  # Output (PDF) file
      output_plot_portrait     = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/TAref4pathways_A375.H_w.v1_PPLOT.v1.pdf")  # Output (PDF) file

# Project a dataset (CCLE) in the space of the components (as defined by the W matrix)

   source("~/CGP2014/TA/TA_PATHWAYS/DISSECTOR_lib.v3.R")

   DISSECTOR_project_dataset.v1(     # Project a dataset in the space defined by a W matrix
      input_dataset           = "~/CGP2014/TA/TA_PATHWAYS/rnaseq.v3.gct",  # Input dataset (GCT)
      input_normalization     = "rank",     # Normalization for the input dataset: "rank"
      normalize_after_match   = T,        # Normalize input dataset after matching with rows of W
      input_W_dataset         = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/TAref4pathways_A375.W.v1.gct",
      W_normalization         = "none",         # Normalization for W                                            
      output_H_dataset        = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/rnaseq.v3_TAref4pathways_A375.H_proj.v1.gct",   # Output dataset H (GCT)
      output_W_dataset        = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/TAref4pathways_A375.W_NORM.v1.gct")  # Output dataset normalized W (GCT)

# Plot the projected dataset (H-like matrix)

   source("~/CGP2014/TA/TA_PATHWAYS/DISSECTOR_lib.v3.R")

   DISSECTOR_make_heatmap_of_matrix.v1(
      input_dataset            = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/rnaseq.v3_TAref4pathways_A375.H_proj.v1.gct", # Input dataset (GCT).
      annot_file               = c("~/CGP2014/TA/TA_PATHWAYS/CCLE_Kim_Tamayo_96AffyCEL_FAZES_Oct12.v2.txt","CCLE_name", "Site.Primary", F),
      transpose_data           = F,           # Transpose input matrix
      append_annot             = T,           # Append annotation to column names
      sorting_method           = "HC",        # Sorting method for cols inside phenotype: MDS (Multi_dimensinal Scaling) or HC (Hiererachical Clustering)  
      output_plot_landscape    = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/rnaseq.v3_TAref4pathways_A375.H_proj.v1_LPLOT.v1.pdf",  # Output (PDF) file
      output_plot_portrait     = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/rnaseq.v3_TAref4pathways_A375.H_proj.v1_PPLOT.v1.pdf")  # Output (PDF) file


   source("~/CGP2013/TA/DISSECTOR_lib.v1.R")

   DISSECTOR_generate_association_matrix.v1(  #  generate the association matrix between the columns (perturbation/states)
      input_dataset            = "~/CGP2014/TA/TA_PATHWAYS/Analysis_V1/TAref4pathways_A375.H_w.v1.gct", # Input dataset (GCT).
      input_dataset2           = "~/CGP2014/TA/TA_PATHWAYS/Analysis_V1/TAref4pathways_A375.H_w.v1.gct", # Input dataset (GCT).
      association_type         = "columns",     # Type of association: between "columns" or between "rows"                                            
      annot_file               = c("~/CGP2014/TA/TA_PATHWAYS/TARef_v2_Yash_20130507.txt", "Gene.Symbol", "Pathway.Cellular.Process", T),
      annot_file2              = c("~/CGP2014/TA/TA_PATHWAYS/TARef_v2_Yash_20130507.txt", "Gene.Symbol", "Pathway.Cellular.Process", T),  
      association_metric       = "IC",   # Association metric: "IC" (Information Coefficient), "COR" (Pearson), "SPEAR" (Spearman)
      output_assoc_matrix_file = "~/CGP2014/TA/TA_PATHWAYS/Analysis_V1/TAref4pathways_A375.H_w_ASSOC.v1.gct", 
      output_assoc_plot        = "~/CGP2014/TA/TA_PATHWAYS/Analysis_V1/TAref4pathways_A375.H_w_ASSOC.v1.pdf")  # Output (PDF) assoc plot

# Find top associations for each ORFs (nearest neighbors)

   source("~/CGP2014/TA/TA_PATHWAYS/DISSECTOR_lib.v3.R")

   DISSECTOR_find_top_associations.v1(   #  Find top associations to each column state 
      input_dataset            = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/TAref4pathways_A375.H_w.v1.gct",
      annot_file               = c("~/CGP2014/TA/TA_PATHWAYS/TARef_v2_Yash_20130507.txt", "Gene.Symbol", "Pathway.Cellular.Process", T),  
      n.top                    = 35,                    # Number of top (UP) matches to display
      n.bottom                 = 10,                    # Number of top (DOWN) matches to display
      output_top_assocs_plot   = "~/CGP2014/TA/TA_PATHWAYS/TAref4pathways_A375.H_TOP_ASSOC.v1.offset.pdf") # Output (PDF) file with top associations

# Annotate the components using genomic datasets

   source("~/CGP2014/TA/TA_PATHWAYS/DISSECTOR_lib.v3.R")

#   DISSECTOR_annotate_components.v1(  #  Match genomic features to each components profile
#      input_dataset         = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/rnaseq.v3_TAref4pathways_A375.H_proj.v1.gct", # Input dataset (GCT).
#      directory             = "~/CGP2013/TA/TA_PATHWAYS/Analysis_v1/",                            
#      identifier            = "Comp_annot.v1",                                                            # string or prefix to identify this analysis
#      feature.type.files    =  list("ACHILLES"     = "~/CGP2014/TA/TA_PATHWAYS/Achilles_v2.4.1.rnai.Gs.gct",
#                                    "RPPA"         = "~/CGP2014/TA/TA_PATHWAYS/RPPA.dat.gct",
#                                    "MUT_CNA"      = "~/CGP2014/TA/TA_PATHWAYS/RNAseqMUTs_CNA_20130729.gct",  
#                                    "EXP_PATHWAYS" = "~/CGP2014/TA/TA_PATHWAYS/CCLE_MSigDB_plus_oncogenic.PATHWAYS.v2.gct",
#                                    "EXP_GENES"    = "~/CGP2014/TA/TA_PATHWAYS/rnaseq.v3.gct"),
#       feature.dir         = c(0, 1, 1, 1, 1),                                         
#       n.markers            = 25,                         # Number of top hits shown in the heatmaps
#       n.perm               = 2,                                         
#       char.scaling         = 0.575,                   # Character scaling for heatmaps
#       locs.table.file      = "~/CGP2014/TA/TA_PATHWAYS/hgnc_downloads.txt",
#       log.table.file       = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/annot_comp.log.txt",                            
#       min.thres            = 10,
#       character.scaling    = 0.65,
#       phen.table           = NULL,
#       phenotypes           = NULL)  


# Produce multi-dimensional scaling biplot projection

   source("~/CGP2014/TA/TA_PATHWAYS/DISSECTOR_lib.v3.R")

   DISSECTOR_produce_MDS_and_network.v1(    # Make MDS projection and network from a GCT file
      input_dataset          = "~/CGP2014/TA/TA_PATHWAYS/Analysis_V1//TAref4pathways_A375.H_w.v1.gct",
      col_annot_file         = c("~/CGP2014/TA/TA_PATHWAYS/TARef_v2_Yash_20130507.txt", "Gene.Symbol", "Pathway.Cellular.Process", T),  
      row_annot_file         = NULL,      # Rows annotation in format c(file, name_column, annot_column, use_prefix)
      transpose_data         = F,         # Transpose input matrix
      append_col_annot       = F,         # Append annotation to column names
      append_row_annot       = F,         # Append annotation to column names
      plot_row_objects       = F,
      plot_row_objects_names = T,
      plot_row_dimensions    = T,
      plot_col_objects       = T,
      plot_col_objects_names = F,
      plot_col_dimensions    = F,
      col_object_size        = 0.25,
      row_object_size        = "auto",
      output_plot            = "~/CGP2014/TA/TA_PATHWAYS/TAref4pathways_A375.H_w.v1.MDS.pdf", 
      MDS_movie_file         = NULL)

### ---------------------------------------------------------------------------------
#
# Project clique drugs into ORF space

   source("~/CGP2014/TA/TA_PATHWAYS/DISSECTOR_lib.v3.R")

   DISSECTOR_project_dataset.v1(     # Project a dataset in the space defined by a W matrix
      input_dataset            = "~/CGP2013/TA/CC/cliques/clique_compound_classes_n2873x978.H_w.v1.gct", # Input dataset (GCT).
      input_normalization     = "rank",     # Normalization for the input dataset: "rank"
      normalize_after_match   = T,        # Normalize input dataset after matching with rows of W
      input_W_dataset         = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/TAref4pathways_A375.W.v1.gct",
      W_normalization         = "none",         # Normalization for W                                            
      output_H_dataset        = "~/CGP2014/TA/TA_PATHWAYS/Analysis_v1/clique_compound_classes_n2873x978.H_TAref4pathways_A375_proj.v1.gct",   # Output dataset H (GCT)
      output_W_dataset        = NULL)  # Output dataset normalized W (GCT)

