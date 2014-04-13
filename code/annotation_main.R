

source("/medpop/mpg-psrl/Dissector/DISSECTOR_lib.v3.R")
source("/medpop/mpg-psrl/Dissector/FS.library.v8.6.R")


   DISSECTOR_annotate_components.v1(
      input_dataset         = "/medpop/mpg-psrl/Dissector/A375_7/A375_2/data/DLBCL_OE_L1000_NMF_TA.OE002_PC3_72H_H.gct",
      directory             = "/medpop/mpg-psrl/Dissector/A375_7/A375_2/annotation",                            
      identifier            = "Comp_annot.v1",
      feature.type.files    =  list("ACHILLES"     = "/medpop/mpg-psrl/Dissector/Achilles_v2.4.1.rnai.Gs.gct",
                                    "RPPA"         = "/medpop/mpg-psrl/Dissector/RPPA.dat.gct",
                                    "MUT_CNA"      = "/medpop/mpg-psrl/Dissector/RNAseqMUTs_CNA_20130729.gct",  
                                    "EXP_PATHWAYS" = "/medpop/mpg-psrl/Dissector/CCLE_MSigDB_plus_oncogenic.PATHWAYS.v2.gct",
                                    "EXP_GENES"    = "/medpop/mpg-psrl/Dissector/rnaseq.v3.gct"),
       feature.directions   = c(0, 1, 1, 1, 1),                                         
       n.markers            = 15,                         # Number of top hits shown in the heatmaps
       n.perm               = 2,                                         
       char.scaling         = 0.575,                   # Character scaling for heatmaps
       locs.table.file      = "/medpop/mpg-psrl/Dissector/hgnc_downloads.txt",
       log.table.file       = "/medpop/mpg-psrl/Dissector/annot_comp.log.txt",                            
       min.thres            = 10,
       character.scaling    = 0.65,
       phen.table           = NULL,
       phenotypes           = NULL)  




