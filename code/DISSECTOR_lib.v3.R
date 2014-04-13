   library(scatterplot3d)
   library(MASS)
   library(RColorBrewer)
#   library(rgl) 
#   library(bpca)
#   library(smacof)


   DISSECTOR_produce_MDS_and_network.v1 <- function(
      #
      # Make MDS projection and network from a GCT file
      #
      input_dataset,                      # Input dataset (GCT). This is e.g. the original dataset A or the H matrix
      col_annot_file         = NULL,      # Columns annotation file (TXT, optional) in format c(file, name_column, annot_column, use_prefix)
      row_annot_file         = NULL,      # Rows annotation file (TXT, optional) in format c(file, name_column, annot_column, use_prefix)
      transpose_data         = F,         # Transpose input matrix
      append_col_annot       = F,         # Append annotation to column names
      append_row_annot       = F,         # Append annotation to column names
      plot_row_objects       = T,
      plot_row_objects_names = F,
      plot_row_dimensions    = F,
      plot_col_objects       = T,
      plot_col_objects_names = T,
      plot_col_dimensions    = T,
      col_object_size        = "auto",
      row_object_size        = "auto",
      output_plot,                                                    
      MDS_movie_file         = NULL)

  {

      set.seed(5209761)
       
      # Read input dataset

      dataset.1 <- MSIG.Gct2Frame(filename = input_dataset)
      H <- data.matrix(dataset.1$ds)
      print(dim(H))

      if (transpose_data == T) H <- t(H)

      width <- ceiling(ncol(H)/100)
      if (width < 11) width <- 11
      pdf(file=output_plot, height=8.5, width=width)
   
      # Read column annotation file
 
      if (!is.null(col_annot_file)) {
         annot.table <- read.table(col_annot_file[[1]], header=T, sep="\t", skip=0, colClasses = "character")
         column.list <- annot.table[, col_annot_file[[2]]]
         annot.list <- annot.table[, col_annot_file[[3]]]
         column.set <- vector(length=ncol(H), mode="character")
         if (col_annot_file[[4]] == T) {
            for (i in 1:ncol(H)) {
               column.set[i] <- strsplit(colnames(H)[i], split="_")[[1]]
            }
         } else {
            column.set <- colnames(H)
         }
         locs <- match(column.set, column.list)
         column.class <- annot.list[locs]
         column.class[is.na(column.class)] <- ""
         for (k in 1:length(column.class)) column.class[k] <- substr(column.class[k], 1, 20)
         all.col.classes <- unique(column.class)
         if (append_col_annot == T) colnames(H) <- paste(colnames(H), " (", column.class, ") ", sep="")
      } else {
         column.class <- rep(" ", ncol(H))
         all.col.classes <- " "
      }

      # Read row annotation file
 
      if (!is.null(row_annot_file)) {
         annot.table <- read.table(row_annot_file[[1]], header=T, sep="\t", skip=0, colClasses = "character")
         column.list <- annot.table[, row_annot_file[[2]]]
         annot.list <- annot.table[, row_annot_file[[3]]]
         row.set <- vector(length=nrow(H), mode="character")
         if (row_annot_file[[4]] == T) {
            for (i in 1:nrow(H)) {
               row.set[i] <- strsplit(row.names(H)[i], split="_")[[1]]
            }
         } else {
            row.set <- row.names(H)
         }
         locs <- match(row.set, column.list)
         row.class <- annot.list[locs]
         row.class[is.na(row.class)] <- ""
         for (k in 1:length(row.class)) row.class[k] <- substr(row.class[k], 1, 20)
         all.row.classes <- unique(row.class)
         if (append_row_annot == T) row.names(H) <- paste(row.names(H), " (", row.class, ") ", sep="")
      } else {
         row.class <- rep(" ", nrow(H))
         all.row.classes <- " "
      }
 
      # Color map
 
      mycol <- vector(length=512, mode = "numeric")
      for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
      for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
      mycol <- rev(mycol)
      max.cont.color <- 511
      mycol <- c(mycol,
              "black",  # brewer.pal(9, "YlGn")[1],     # Missing feature color (light yellow)
              mycol[256 - 75],                          # Binary feature's 0's color (blue)
              mycol[256 + 220])                         # Binary feature's 1's color (red)
      cex.axis = 1
      phen.col <- c("steelblue", brewer.pal(7, "Set1"), brewer.pal(7, "Pastel1"), brewer.pal(7, "Dark2"), brewer.pal(7, "Paired"),
                 brewer.pal(7, "Pastel2"), brewer.pal(8, "Accent"), brewer.pal(8, "Set2"), brewer.pal(11, "Spectral"), brewer.pal(12, "Set3"),
                 sample(c(brewer.pal(9, "Blues"), brewer.pal(9, "Reds"), brewer.pal(9, "Oranges"), brewer.pal(9, "Greys"),
                 brewer.pal(9, "Purples"), brewer.pal(9, "Greens"))))
      row.names <- row.names(H)
      col.names <- colnames(H)
      row.names(H) <- paste("    ", row.names, sep="")
      colnames(H) <- paste("    ", col.names, sep="")

      bpca.1 <- bpca(t(H), method='hj', lambda.end=3)

      col.objects <- bpca.1$coord$objects
      row.objects <- bpca.1$coord$variables
      row.objects <- 1.10 * (max(col.objects)/ max(row.objects)) * row.objects

      open3d(windowRect = c(10, 10, 1200, 1200), zoom=0.75)

      if (col_object_size == "auto") {
         col.radius <- 0.15 + 7.5 / nrow(col.objects) 
      } else {
         col.radius <- col_object_size
      }
      for (i in 1:nrow(col.objects)) {
         obj.col <- phen.col[match(column.class, all.col.classes)]
         if (plot_col_objects == T) spheres3d(col.objects[i,1], col.objects[i,2], col.objects[i,3], radius= col.radius, color=obj.col[i])
         if (plot_col_objects_names == T) text3d(1.2*col.objects[i,1], 1.2*col.objects[i,2], 1.2* col.objects[i,3],
                                              texts=row.names(col.objects)[i], color="blue", cex=1.5)
         if (plot_col_dimensions == T)  lines3d(c(0,col.objects[i,1]), c(0, col.objects[i,2]), c(0, col.objects[i,3]), lwd=1, color="blue")
      }
      if (row_object_size == "auto") {
        row.radius <- 0.15 + 7.5 / nrow(row.objects) 
      } else {
        row.radius <- row_object_size
      }
      for (i in 1:nrow(row.objects)) {
         row.col <- phen.col[match(row.class, all.row.classes)]
         if (plot_row_objects == T) spheres3d(row.objects[i,1], row.objects[i,2], row.objects[i,3], radius= row.radius, color=row.col[i])
         if (plot_row_objects_names == T)  text3d(1.075*row.objects[i,1], 1.075*row.objects[i,2], 1.075*row.objects[i,3],
                                               texts=row.names(row.objects)[i], color="red", cex=1.5) 
         if (plot_row_dimensions == T)  lines3d(c(0,row.objects[i,1]), c(0, row.objects[i,2]), c(0, row.objects[i,3]), lwd = 1, color="red")
      }

      s <- strsplit(movie_file, split="/")
      movie.name <- s[[1]][length(s[[1]])]
      dir <- paste(s[[1]][seq(1, length(s[[1]])-1)], collapse="/")
      dir <- paste(dir, "/", sep="")

#   my.movie3d(spin3d(axis=c(1,1,1), rpm=2), duration=5, fps = 10, movie = movie_file, dir = dir,
#          convert = TRUE, clean = TRUE, verbose=TRUE, top = TRUE, type = "gif", startTime = 0)

      dev.off()
   }

   DISSECTOR_annotate_components.v1 <- function(
      #
      #  Match genomic features to each components profile
      #
      input_dataset,                           # Input dataset (GCT). This is e.g. the H matrix
      directory,                               # Directory where to produce results                                          
      identifier            = "Comp_annot.v1", # string or prefix to identify this analysis
      feature.type.files,                      # List of feature typers and correponding file e.g.
                                               # list("ACHILLES"     = "~/CGP2013/Distiller/Achilles_v2.4.1.rnai.Gs.gct",
                                               #  "MUT_CNA"      = "~/CGP2013/Distiller/RNAseqMUTs_CNA_20130729.gct",  
                                               #  "EXP_PATHWAYS" = "~/CGP2013/Distiller/CCLE_MSigDB_plus_oncogenic.PATHWAYS.v2.gct",
                                               #  "EXP_GENES"    = "~/CGP2013/CCLE/rnaseq.v3.gct",
                                               #  "RPPA"         = "~/CGP2013/Distiller/RPPA.dat.gct")
       feature.directions,                     # c(0, 1, 1, 1, 1),                                         
       n.markers            = 25,                         # Number of top hits shown in the heatmaps
       n.perm               = 10,                                         
       char.scaling         = 0.575,                   # Character scaling for heatmaps
       locs.table.file      = NULL,
       log.table.file,
       min.thres            = 10,
       character.scaling    = 0.65,
       phen.table           = NULL,
       phenotypes           = NULL)            # list(list("ALL"))
     
   {

      version <- identifier
      target.file <- input_dataset
      dataset.1 <- MSIG.Gct2Frame(filename = target.file)
      
      H <- data.matrix(dataset.1$ds)
      targets <- row.names(H)
      target.types <- rep("COMP", length(targets))
      target.dir <- rep(1, length(targets))
      feature.types <- names(feature.type.files)
      n.targets <- length(targets)
      n.f.types <- length(feature.types)
      
      if (!is.null(phen.table)) {
         samples.table <- read.delim(phen.table, header=T, row.names=1, sep="\t", skip=0)
         table.names <- row.names(samples.table)
       }
      log.table <- NULL

      for (d in 1:n.targets) {  # loop over targets

         target.dir.name <- paste(directory, targets[d], sep="")
         dir.create(target.dir.name, showWarnings=FALSE)

         for (f in 1:n.f.types) {   # loop over feature types
        
            t.dir <- ifelse (length(target.dir) == 1, target.dir[[1]], target.dir[[d]])
           
            dir <- ifelse(xor(t.dir, feature.directions[[f]]), "negative", "positive")

            print(paste("Reading target file:", target.file))
            dataset.1 <- MSIG.Gct2Frame(filename = target.file)
            sample.names.1 <- dataset.1$names
            
            print(paste("Reading features file:", feature.type.files[[f]]))
            dataset.2 <- MSIG.Gct2Frame(filename = feature.type.files[[f]])
            sample.names.2 <- dataset.2$names

            if (is.null(phen.table)) {
               phen.set <- 1
               phen.selected <- NULL
               results.file.pdf <- paste(directory, targets[d], "/", targets[d], "_vs_", feature.types[[f]], "_", version, ".pdf", sep="")
               results.file.txt <- paste(directory, targets[d], "/", targets[d], "_vs_", feature.types[[f]], "_", version, ".txt", sep="")               
            } else {
              results.file.pdf <- paste(directory, targets[d], "/", phen, "_", targets[d], "_vs_", feature.types[[f]], "_", version, ".pdf", sep="")
              results.file.txt <- paste(directory, targets[d], "/", phen, "_", targets[d], "_vs_", feature.types[[f]], "_", version, ".txt", sep="")               
              if (length(phenotypes) == 1) {
                  phen.set <- unlist(phenotypes[[1]])
                  phen.selected <- "1"
               } else {
                  phen.set <- unlist(phenotypes[[d]])
                  phen.selected <- "1"
               }
             }
            for (p in 1:length(phen.set)) {  # loop over phenotypes
               if (phen.set == 1) {
                  phen <- NULL
                  print(paste("Working on: target=", d, " feature=", f))
                  overlap <- intersect(sample.names.1, sample.names.2)
                  doc.string <- c(ifelse(length(overlap) < min.thres, "NO", "YES"), length(overlap), targets[d], dir, feature.types[[f]],
                                  version, target.file, feature.type.files[[f]])
                } else {
                  phen <- phen.set[p]
                  print(paste("Working on: target=", d, " feature=", f, " phen=", phen))
                  sample.names.3 <- table.names[samples.table[,  phen] == "1"]
                  overlap <- intersect(sample.names.1, intersect(sample.names.2, sample.names.3))
                  doc.string <- c(ifelse(length(overlap) < min.thres, "NO", "YES"), length(overlap), phen, targets[d], dir, feature.types[[f]],
                                  version, target.file, feature.type.files[[f]])
               }

               log.table <- rbind(log.table, doc.string)
               
               # Only run analysis if there are at least min.thres samples
         
               if (length(overlap) < min.thres) next
             
               SCREENER.v1(
                  ds1 =                            target.file,
                  target.name =                    targets[d],
                  ds2 =                            feature.type.files[[f]],
                  n.markers =                      n.markers,             
                  n.perm =                         n.perm,           
                  permutation.test.type =          "standard",  
                  n.boot =                         200,         
                  seed =                           2345971, 
                  assoc.metric.type =              "IC",      
                  direction =                      dir,
                  sort.target =                    TRUE,           
                  results.file.pdf =               results.file.pdf,
                  results.file.txt =               results.file.txt,
                  sort.columns.inside.classes =    F,
                  locs.table.file =                locs.table.file,
                  consolidate.identical.features = "identical",
                  cons.features.hamming.thres =    0,   
                  save.matched.dataset =           F,  
                  produce.aux.histograms =         F,
                  produce.heat.map =               T,
                  produce.mds.plots =              T,
                  character.scaling =              character.scaling,
                  mds.plot.type =                  "smacof",
                  n.grid =                         25,
                  phen.table =                     phen.table,
                  phen.column =                    phen,
                  phen.selected =                  phen.selected)
            
              } # loop over phenotypes
  
          } # loop over feature types
    
       } # loop over targets

       # Save log records
      
       if (phen.set == 1) {
          header <- c("Run:", "Samples:", "target:", "Direction:", "Feature.type:", "Version:", "target.file:", "Feature.type.file:")
        } else {
          header <- c("Run:", "Samples:", "Phenotype:", "target:", "Direction:", "Feature.type:", "Version:", "target.file:", "Feature.type.file:")
        }
        colnames(log.table) <- header
        write.table(log.table, file=log.table.file, quote=F, col.names = T, row.names = F, append = F, sep="\t")
    }


   DISSECTOR_find_top_associations.v1 <- function(
      #
      #  Find top associations to each column state 
      #
      input_dataset,                    # Input dataset (GCT). This is e.g. the original dataset A or the H matrix
      annot_file = NULL,                # Gene nnotation file (TXT, optional) in format c(file, name_column, annot_column, use_prefix)
      n.top = 35,                       # Number of top (UP) matches to display
      n.bottom = 10,                    # Number of top (DOWN) matches to display
      output_top_assocs_plot)           # Output (PDF) file with top associations
   {

   set.seed(5209761)
            
   # Read input dataset

   dataset.1 <- MSIG.Gct2Frame(filename = input_dataset)
   H <- data.matrix(dataset.1$ds)
   print(dim(H))
   k.comp <- nrow(H)

   # Read annotation file

   if (!is.null(annot_file)) {
      annot.table <- read.table(annot_file[[1]], header=T, sep="\t", skip=0, colClasses = "character")
      gene.list <- annot.table[, annot_file[[2]]]
      annot.list <- annot.table[, annot_file[[3]]]
      gene.set <- vector(length=ncol(H), mode="character")

      if (annot_file[[4]] == T) {
         for (i in 1:ncol(H)) {
            gene.set[i] <- strsplit(colnames(H)[i], split="_")[[1]]
         }
      } else {
         gene.set <- colnames(H)
      }
      locs <- match(gene.set, gene.list)
      gene.class <- annot.list[locs]
      for (k in 1:length(gene.class)) gene.class[k] <- substr(gene.class[k], 1, 10)
      all.classes <- unique(gene.class)
      colnames(H) <- paste(colnames(H), " (", gene.class, ") ", sep="")
   }
   
   Ht <- t(H)

   mycol <- vector(length=512, mode = "numeric")   # Red/Blue "pinkogram" color map
   for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol <- rev(mycol)
   ncolors <- length(mycol)

   pdf(file=output_top_assocs_plot, height=8.5, width=11)

   for (i in 1:nrow(Ht)) {
      nf <- layout(matrix(c(1, 2, 3, 4), 4, 1, byrow=T), 1, c(4, ceiling(n.top/2) + 2, ceiling(n.bottom/2) + 2, 3),  FALSE)
      ind <- order(Ht[i,], decreasing=T)
      Ht.temp <- Ht[, ind]
      mutinf.m <- NULL
      for (k in 1:nrow(Ht.temp)) mutinf.m <- c(mutinf.m, IC.v1(Ht.temp[i,], Ht.temp[k,]))
      ind <- order(mutinf.m, decreasing=T)
      mutinf.m <- signif(mutinf.m[ind], 3)
      Ht.temp <- Ht.temp[ind,]
      gene.class.temp <- gene.class[ind]
     
      # Target profile
    
      V <- Ht.temp[1,]
      gc <- gene.class.temp[1]
      V <- V*ncolors
      par(mar = c(1, 20, 4, 6))    
      image(1:length(V), 1:1, as.matrix(V), zlim = c(0, ncolors), col=mycol, axes=FALSE, main=row.names(Ht)[i], sub = "", cex.main = 2, xlab= "", ylab="")

      if (!is.null(annot_file)) {
         col.classes <- c(brewer.pal(7, "Set1"), brewer.pal(7, "Pastel1"), brewer.pal(7, "Dark2"), brewer.pal(7, "Paired"), brewer.pal(7, "Pastel2"),
                 brewer.pal(8, "Accent"), brewer.pal(8, "Set2"), brewer.pal(11, "Spectral"), brewer.pal(12, "Set3"),
                 sample(c(brewer.pal(9, "Blues"), brewer.pal(9, "Reds"), brewer.pal(9, "Oranges"), brewer.pal(9, "Greys"),
                          brewer.pal(9, "Purples"), brewer.pal(9, "Greens"))))
         gc <- gene.class[ind]
         cols <- col.classes[match(gc, all.classes)]
      } else {
         cols = "black"
      }
      mtext(row.names(Ht)[i], at=1, side = 2, cex=0.50, col=cols, line=0, las=1, font=2, family="")
      axis(4, at=1+0.2, labels="        IC  ", adj= 0.5, tick=FALSE, las = 1, cex.axis=0.85, line=-1, font=2, family="")
      axis(4, at=1, labels="   (Information ", adj= 0.5, tick=FALSE, las = 1, cex.axis=0.75, line=-1, font=2, family="")
      axis(4, at=1-0.2, labels="    Coefficient) ", adj= 0.5, tick=FALSE, las = 1, cex.axis=0.75, line=-1, font=2, family="")     
      axis(3, at=1:length(V), labels=colnames(Ht.temp), adj= 0.5, tick=FALSE, las = 1, cex.axis=ifelse(k.comp >= 25, 0.6, 0.85),
           font.axis=1, line=-1, font=2, family="")
     
      # HIGHEST IC ASSOCIATIONS

      V <- Ht.temp[2:(n.top+1),]
      gc <- gene.class.temp[2:(n.top+1)]
      row.names(V) <- paste(row.names(V), " (", gc, ") ", sep="")
      muti <- mutinf.m[2:(n.top+1)]
      V <- V*ncolors
      V <- apply(V, MARGIN=2, FUN=rev)
      gc <- rev(gc)
      muti <- rev(muti)
      par(mar = c(1, 20, 3, 6))
      image(1:dim(V)[2], 1:dim(V)[1], t(V), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="HIGHEST IC ASSOCIATIONS", sub = "", xlab= "", ylab="")     
      mtext(row.names(V), at=1:nrow(V), side = 2, cex=0.50, col=col.classes[match(gc, all.classes)], line=0, las=1, font=2, family="")
      axis(4, at=1:nrow(V), labels=muti, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.85, font.axis=1, line=0, font=2, family="")
      axis(3, at=1:ncol(V), labels=colnames(V), adj= 0.5, tick=FALSE, las = 1, cex.axis=ifelse(k.comp >= 25, 0.6, 0.85),
           font.axis=1, line=-1, font=2, family="")

      # LOWEST IC ASSOCIATIONS
     
      V <- Ht.temp[seq(nrow(Ht.temp) - n.bottom + 1, nrow(Ht.temp)),]
      gc <- gene.class.temp[seq(nrow(Ht.temp) - n.bottom + 1, nrow(Ht.temp))]
      row.names(V) <- paste(row.names(V), " (", gc, ") ", sep="")
      muti <- mutinf.m[seq(nrow(Ht.temp) - n.bottom + 1, nrow(Ht.temp))]
      V <- V*ncolors
      V <- apply(V, MARGIN=2, FUN=rev)
      gc <- rev(gc)
      muti <- rev(muti)     
      par(mar = c(2, 20, 3, 6))     
      image(1:dim(V)[2], 1:dim(V)[1], t(V), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="LOWEST IC ASSOCIATIONS", sub = "", xlab= "", ylab="")     
      mtext(row.names(V), at=1:nrow(V), side = 2, cex=0.50, col=col.classes[match(gc, all.classes)], line=0, las=1, font=2, family="")    
      axis(4, at=1:nrow(V), labels=muti, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.85, font.axis=1, line=0, font=2, family="")
      axis(3, at=1:ncol(V), labels=colnames(V), adj= 0.5, tick=FALSE, las = 1, cex.axis=ifelse(k.comp >= 25, 0.6, 0.85),
           font.axis=1, line=-1, font=2, family="")

      par.mar <- par("mar")
      par(mar = c(4, 45, 1, 5))
      leg.set <- seq(0, 1, 0.01)
      image(1:length(leg.set), 1:1, as.matrix(leg.set), zlim=c(0, 1), col=mycol, axes=FALSE, main=paste("Legend"), font=2, family="")
      ticks <- seq(0, 1, 0.1)
      tick.cols <- rep("black", 5)
      tick.lwd <- 2
      locs <- NULL
      for (k in 1:length(ticks)) locs <- c(locs, which.min(abs(ticks[k] - leg.set)))
      axis(1, at=locs, labels=ticks, adj= 0.5, tick=T, cex=0.7, cex.axis=1, line=0, font=2, family="")
      mtext("Component Amplitude", cex=0.85, side = 1, line = 2.5, outer=F)     
      par(mar = par.mar)

   }
   dev.off()
 }

   DISSECTOR_make_heatmap_of_matrix.v1 <- function(
      #
      # Make heatmap of matrix from GCT file
      #
      input_dataset,           # Input dataset (GCT). This is e.g. the original dataset A or the H matrix
      annot_file = NULL,       # Phenotype annotation file (TXT, optional) in format c(file, name_column, annot_column, use_prefix)
      transpose_data = F,      # Transpose input matrix
      append_annot = F,        # Append annotation to column names
      sorting_method = "MDS",  # Sorting method for columns inside a phenotype: MDS (Multi_dimensinal Scaling) or HC (Hiererachical Clustering)
      output_plot_landscape,   # Output (PDF) plot in landscape format
      output_plot_portrait)    # Output (PDF) plot in portrait format
  {

   set.seed(5209761)
       
   # Read input dataset

   dataset.1 <- MSIG.Gct2Frame(filename = input_dataset)
   H <- data.matrix(dataset.1$ds)
   print(dim(H))

   if (transpose_data == T) H <- t(H)

   width <- ceiling(ncol(H)/100)
   if (width < 11) width <- 11
   pdf(file=output_plot_landscape, height=8.5, width=width)
   
   # Read annotation file
 
   if (!is.null(annot_file)) {
      annot.table <- read.table(annot_file[[1]], header=T, sep="\t", skip=0, colClasses = "character")
      column.list <- annot.table[, annot_file[[2]]]
      annot.list <- annot.table[, annot_file[[3]]]
      column.set <- vector(length=ncol(H), mode="character")

      if (annot_file[[4]] == T) {
         for (i in 1:ncol(H)) {
            column.set[i] <- strsplit(colnames(H)[i], split="_")[[1]]
         }
      } else {
         column.set <- colnames(H)
      }

#print(column.set)
#print(annot.list)



      locs <- match(column.set, column.list)
print("c.list")
print(cbind(column.set, column.list))


      column.class <- annot.list[locs]

      column.class[is.na(column.class)] <- "UNLABELED"
      for (k in 1:length(column.class)) column.class[k] <- substr(column.class[k], 1, 20)
      all.classes <- unique(column.class)




      if (append_annot == T) colnames(H) <- paste(colnames(H), " (", column.class, ") ", sep="")
   } else {
      column.class <- rep(" ", ncol(H))
      all.classes <- " "
   }
  
   # Color map

   mycol <- vector(length=512, mode = "numeric")
   for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol <- rev(mycol)
   max.cont.color <- 511
   mycol <- c(mycol,
              "black",  # brewer.pal(9, "YlGn")[1],     # Missing feature color (light yellow)
              mycol[256 - 75],                          # Binary feature's 0's color (blue)
              mycol[256 + 220])                         # Binary feature's 1's color (red)
   cex.axis = 1
   phen.col <- c(brewer.pal(7, "Set1"), brewer.pal(7, "Pastel1"), brewer.pal(7, "Dark2"), brewer.pal(7, "Paired"), brewer.pal(7, "Pastel2"),
                 brewer.pal(8, "Accent"), brewer.pal(8, "Set2"), brewer.pal(11, "Spectral"), brewer.pal(12, "Set3"),
                 sample(c(brewer.pal(9, "Blues"), brewer.pal(9, "Reds"), brewer.pal(9, "Oranges"), brewer.pal(9, "Greys"),
                          brewer.pal(9, "Purples"), brewer.pal(9, "Greens"))))
 
    # Normalize and apply color map

   cutoff <- 2
   for (i in 1:nrow(H)) {
      x <- H[i,]
      locs.non.na <- !is.na(x)
      x.nonzero <- x[locs.non.na]
      x.nonzero2 <- (x.nonzero - mean(x.nonzero))/sd(x.nonzero)         
      x.nonzero2[x.nonzero2 > cutoff] <- cutoff
      x.nonzero2[x.nonzero2 < - cutoff] <- - cutoff      
      s <- strsplit(row.names(H)[i], "_")[[1]]
      suffix <- s[length(s)]
      if (suffix == "MUT" | suffix == "AMP" | suffix == "DEL" | suffix == "AMP_2" | suffix == "AMP_3" | suffix == "DEL_2" | suffix == "DEL_3" |
          suffix == "all" | length(table(x.nonzero)) == 2) {  # Binary feature
         H[i,locs.non.na] <- x.nonzero + max.cont.color + 2   # binary feature colors
       } else {
         H[i, locs.non.na] <- x.nonzero2
         H[i, locs.non.na] <- ceiling(max.cont.color * (H[i,locs.non.na] + cutoff)/(2*cutoff))
         H[i, locs.non.na] <- ifelse (H[i, locs.non.na] > max.cont.color, max.cont.color, H[i, locs.non.na])
       }
      H[i, is.na(x)] <- max.cont.color + 1 # missing feature color 
    }

   # Horizontal/Landscape versions
   
   # Plot sorted by phenotype ------------------------------------------------------------------------------

   nf <- layout(matrix(c(1, 2, 3), nrow=3, ncol=1, byrow=T), 1, c(1, 6.5, 0.75), FALSE)

   H_orig <- H
   column.class_orig <- column.class
   ind <- order(column.class, decreasing=F)
   H <- H[, ind]
   column.class <- column.class[ind]
   all.classes <- unique(column.class)

   V1.phen <- match(column.class, all.classes)
   par(mar = c(1, 10, 2, 6))
   image(1:length(V1.phen), 1:1, as.matrix(V1.phen), col=phen.col[1:max(V1.phen)], axes=FALSE, main="Phenotype", sub = "", xlab= "", ylab="")
   axis(2, at=1:1, labels="Phenotype", adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)


   leg.txt <- all.classes
   for (i in 1:length(leg.txt)) leg.txt[i] <- substr(leg.txt[i], 1, 25)
   boundaries <- NULL
   for (i in 2:length(column.class)) {
      if (column.class[i] != column.class[i-1]) boundaries <- c(boundaries, i-1)
    }


   boundaries <- c(boundaries, length(column.class))    

   locs.bound <- c(boundaries[1]/2, boundaries[2:length(boundaries)] - (boundaries[2:length(boundaries)] - boundaries[1:(length(boundaries)-1)])/2)
   for (i in 1:length(leg.txt)) {
      text.size <- 0.05 + 7/ifelse(nchar(leg.txt[i]) < 6, 6, nchar(leg.txt[i]))
      text(locs.bound[i], 1, labels=leg.txt[i], adj=c(0.5, 1), srt=90, cex=text.size)
    }
 
  # Sort columns inside each phenotypic class 


   for (k in all.classes) {
      if (sum(column.class == k) <= 1) next;
      V1 <- H[, column.class == k]
      dist.matrix <- dist(t(V1))
      if (sorting_method == "MDS") {
         s <- smacofSym(dist.matrix, ndim=1)
         ind <- order(s$conf, decreasing=T)
       } else if (sorting_method == "HC") {
         HC <- hclust(dist.matrix, method="ward")
        ind <- HC$order
       } else {
          stop(paste("ERROR: unknown sorting method:", sorting_method))
       }
      V1 <- V1[ , ind]
      H[, column.class == k] <- V1
   }

   V2 <- apply(H, MARGIN=2, FUN=rev)

   lower.space <-  ceiling(4 + 100/nrow(H))
   par(mar = c(lower.space, 10, 2, 6))

   image(1:ncol(V2), 1:nrow(V2), t(V2), zlim = c(0, max.cont.color + 3), col=mycol, axes=FALSE, main="Matrix Sorted by Phenotype",
         sub = "", xlab= "", ylab="", cex.main=1)


   cex.rows <- 0.20 + 200/(nrow(V2) * max(nchar(row.names(V2))) + 200)



   axis(2, at=1:nrow(V2), labels=row.names(V2), adj= 0.5, tick=FALSE, las = 1, cex.axis=cex.rows, font.axis=1, line=-1)

    print("boundries is ")
	print(boundaries)
  for (i in 1:(length(boundaries)-1)) 
{
print(i)
print("x is ")
print(c(boundaries[i]+0.5, boundaries[i]+0.5))
print("y is ")
print(c(0.5, nrow(V2) + 0.5))

lines(c(boundaries[i]+0.5, boundaries[i]+0.5), c(0.5, nrow(V2) + 0.5), lwd=2, lty=1, col="black")
}
    print("a-here")

   if (!is.null(annot_file)) {
      cols2 <- phen.col[match(column.class, all.classes)]
   } else {
      cols2 = "black"
   }


   cex.cols <- 0.20 + 200/(ncol(V2) * max(nchar(colnames(V2))) + 200)
   mtext(colnames(V2), at=1:ncol(V2), side = 1, cex=cex.cols, col=cols2, line=0, las=3, font=2, family="")
   
   # Legend

   par(mar = c(3, 35, 1, 6))


   leg.set <- seq(-cutoff, cutoff, 0.05)

   image(1:length(leg.set), 1:1, as.matrix(leg.set), zlim=c(-cutoff, cutoff), col=mycol, axes=FALSE, main="Matrix Standardized Profile",
       sub = "", xlab= "", ylab="",font=2, family="", mgp = c(0, 0, 0), cex.main=0.8)


   ticks <- seq(-cutoff, cutoff, 0.5)
   tick.cols <- rep("black", 5)
   tick.lwd <- 1
   locs <- NULL
   for (k in 1:length(ticks)) locs <- c(locs, which.min(abs(ticks[k] - leg.set)))
   axis(1, at=locs, labels=ticks, adj= 0.5, tick=T, cex=0.6, cex.axis=0.6, line=0, font=2, family="", mgp = c(0.1, 0.1, 0.1))

   # Plot sorted by columns and rows ------------------------------------------------------------------------------




   nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(8, 1), FALSE)

   hc <- hclust(dist(t(H)), "complete")
   ind1 <- hc$order
   hc2 <- hclust(dist(H), "complete")
   ind2 <- hc2$order
   H <- H[ind2, ind1]
   column.class <- column.class[ind1]

   V3 <- apply(H, MARGIN=2, FUN=rev)
   lower.space <-  ceiling(4 + 100/nrow(H))
   par(mar = c(lower.space, 10, 2, 6))
   image(1:ncol(V3), 1:nrow(V3), t(V3), zlim = c(0, max.cont.color + 3), col=mycol, axes=FALSE, main="Matrix Sorted (Rows and Columns)",
         sub = "", xlab= "", ylab="", cex.main=1)
   cex.rows <- 0.20 + 200/(nrow(V3) * max(nchar(row.names(V3))) + 200)
   axis(2, at=1:nrow(V3), labels=row.names(V3), adj= 0.5, tick=FALSE, las = 1, cex.axis=cex.rows, font.axis=1, line=-1)   
   cex.cols <- 0.20 + 200/(ncol(V3) * max(nchar(colnames(V3))) + 200)   
#   axis(1, at=1:ncol(V), labels=colnames(V), adj= 0.5, tick=FALSE, las = 3, cex.axis=cex.cols, font.axis=1, line=-1) # Add sample names        

   if (!is.null(annot_file)) {
      cols <- phen.col[match(column.class, all.classes)]
   } else {
      cols = "black"
   }
   cex.cols <- 0.20 + 200/(ncol(V3) * max(nchar(colnames(V3))) + 200)
   mtext(colnames(V3), at=1:ncol(V3), side = 1, cex=cex.cols, col=cols, line=0, las=3, font=2, family="")

   # Legend

   par(mar = c(3, 35, 1, 6))
   leg.set <- seq(-cutoff, cutoff, 0.05)
   image(1:length(leg.set), 1:1, as.matrix(leg.set), zlim=c(-cutoff, cutoff), col=mycol, axes=FALSE, main="Matrix Standardized Profile",
       sub = "", xlab= "", ylab="",font=2, family="", mgp = c(0, 0, 0), cex.main=0.8)
   ticks <- seq(-cutoff, cutoff, 0.5)
   tick.cols <- rep("black", 5)
   tick.lwd <- 1
   locs <- NULL
   for (k in 1:length(ticks)) locs <- c(locs, which.min(abs(ticks[k] - leg.set)))
   axis(1, at=locs, labels=ticks, adj= 0.5, tick=T, cex=0.6, cex.axis=0.6, line=0, font=2, family="", mgp = c(0.1, 0.1, 0.1))

   # Plot in original order ------------------------------------------------------------------------------

   H <- H_orig
   column.class <- column.class_orig
   
   nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(8, 1), FALSE)

   V4 <- apply(H, MARGIN=2, FUN=rev)
   lower.space <-  ceiling(4 + 100/nrow(H))
   par(mar = c(lower.space, 10, 2, 6))
   image(1:ncol(V4), 1:nrow(V4), t(V4), zlim = c(0, max.cont.color + 3), col=mycol, axes=FALSE, main="Matrix in Original Order",
         sub = "", xlab= "", ylab="", cex.main=0.8)
   cex.rows <- 0.20 + 200/(nrow(V4) * max(nchar(row.names(V4))) + 200)   
   axis(2, at=1:nrow(V4), labels=row.names(V4), adj= 0.5, tick=FALSE, las = 1, cex.axis=cex.rows, font.axis=1, line=-1)
   cex.cols <- 0.20 + 200/(ncol(V4) * max(nchar(colnames(V4))) + 200)   



   if (!is.null(annot_file)) {
      cols <- phen.col[match(column.class, all.classes)]
   } else {
      cols = "black"
   }
   cex.cols <- 0.20 + 200/(ncol(V4) * max(nchar(colnames(V4))) + 200)
   mtext(colnames(V4), at=1:ncol(V4), side = 1, cex=cex.cols, col=cols, line=0, las=3, font=2, family="")

   # Legend

   par(mar = c(3, 32, 1, 6))
   leg.set <- seq(-cutoff, cutoff, 0.05)
   image(1:length(leg.set), 1:1, as.matrix(leg.set), zlim=c(-cutoff, cutoff), col=mycol, axes=FALSE, main="Matrix Standardized Profile",
       sub = "", xlab= "", ylab="",font=2, family="", mgp = c(0, 0, 0), cex.main=0.8)
   ticks <- seq(-cutoff, cutoff, 0.5)
   tick.cols <- rep("black", 5)
   tick.lwd <- 1
   locs <- NULL
   for (k in 1:length(ticks)) locs <- c(locs, which.min(abs(ticks[k] - leg.set)))
   axis(1, at=locs, labels=ticks, adj= 0.5, tick=T, cex=0.6, cex.axis=0.6, line=0, font=2, family="", mgp = c(0.1, 0.1, 0.1))

   dev.off()

   # Vertical/Portrait versions
   
   # Plot sorted by phenotype ------------------------------------------------------------------------------

   height <- ceiling(ncol(H)/50)
   if (height < 11) height <- 14
   pdf(file=output_plot_portrait, height=height, width=8.5)
   
   nf <- layout(matrix(c(1, 2, 0, 3), nrow=2, ncol=2, byrow=T), c(2, 5), c(12, 1), FALSE)
   par(mar = c(3, 2, 4, 2))
   image(1:1, 1:length(V1.phen), t(as.matrix(V1.phen)), col=phen.col[1:max(V1.phen)], axes=FALSE, main="Phenotype", sub = "", xlab= "", ylab="")
   for (i in 1:length(leg.txt)) {
      text.size <- 0.05 + 7/ifelse(nchar(leg.txt[i]) < 6, 6, nchar(leg.txt[i]))
      text(1, locs.bound[i], labels=leg.txt[i], adj=c(0.5, 0), srt=0, cex=text.size)
    }
   right.space <-  ceiling(2 + 100/nrow(H))
   par(mar = c(3, 8, 4, right.space))
   V2 <- apply(V2, MARGIN=2, FUN=rev)
   image(1:nrow(V2), 1:ncol(V2), V2, zlim = c(0, max.cont.color + 3), col=mycol, axes=FALSE, main="Matrix Sorted by Phenotype",
         sub = "", xlab= "", ylab="", cex.main=1)
   cex.rows <- 0.20 + 150/(nrow(V2) * max(nchar(row.names(V2))) + 200)

   for (i in 1:(length(boundaries)-1)) lines(c(0.5, ncol(V2) + 0.5), c(boundaries[i]+0.5, boundaries[i]+0.5), lwd=2, lty=1, col="black")   
   axis(3, at=1:nrow(V2), labels=row.names(V2), adj= 0.5, tick=FALSE, las = 1, cex.axis=cex.rows, font.axis=1, line=-1)   
   cex.cols <- 0.20 + 200/(ncol(V2) * max(nchar(colnames(V2))) + 200)
   mtext(colnames(V2), at=1:ncol(V2), side = 2, cex=cex.cols, col=cols2, line=0, las=1, font=2, family="")

   # Legend

   par(mar = c(3, 5, 1, 25))
   leg.set <- seq(-cutoff, cutoff, 0.05)
   image(1:length(leg.set), 1:1, as.matrix(leg.set), zlim=c(-cutoff, cutoff), col=mycol, axes=FALSE, main="Matrix Standardized Profile",
       sub = "", xlab= "", ylab="",font=2, family="", mgp = c(0, 0, 0), cex.main=0.8)
   ticks <- seq(-cutoff, cutoff, 0.5)
   tick.cols <- rep("black", 5)
   tick.lwd <- 1
   locs <- NULL
   for (k in 1:length(ticks)) locs <- c(locs, which.min(abs(ticks[k] - leg.set)))
   axis(1, at=locs, labels=ticks, adj= 0.5, tick=T, cex=0.6, cex.axis=0.6, line=0, font=2, family="", mgp = c(0.1, 0.1, 0.1))

   # Plot sorted by columns and rows ------------------------------------------------------------------------------
   
   nf <- layout(matrix(c(1, 2), nrow=2, ncol=1, byrow=T), 1, c(12, 1), FALSE)
   right.space <-  ceiling(2 + 100/nrow(H))
   par(mar = c(3, 8, 4, right.space))
   V3 <- apply(V3, MARGIN=2, FUN=rev)
   image(1:nrow(V3), 1:ncol(V3), V3, zlim = c(0, max.cont.color + 3), col=mycol, axes=FALSE,  main="Matrix Sorted (Rows and Columns)",
         sub = "", xlab= "", ylab="", cex.main=1)
   cex.rows <- 0.20 + 150/(nrow(V3) * max(nchar(row.names(V3))) + 200)      
   axis(3, at=1:nrow(V3), labels=row.names(V3), adj= 0.5, tick=FALSE, las = 1, cex.axis=cex.rows, font.axis=1, line=-1)
   cex.cols <- 0.20 + 200/(ncol(V3) * max(nchar(colnames(V3))) + 200)   
   mtext(colnames(V3), at=1:ncol(V3), side = 2, cex=cex.cols, col=cols, line=0, las=1, font=2, family="")

   # Legend

   par(mar = c(3, 5, 1, 25))
   leg.set <- seq(-cutoff, cutoff, 0.05)
   image(1:length(leg.set), 1:1, as.matrix(leg.set), zlim=c(-cutoff, cutoff), col=mycol, axes=FALSE, main="Matrix Standardized Profile",
       sub = "", xlab= "", ylab="",font=2, family="", mgp = c(0, 0, 0), cex.main=0.8)
   ticks <- seq(-cutoff, cutoff, 0.5)
   tick.cols <- rep("black", 5)
   tick.lwd <- 1
   locs <- NULL
   for (k in 1:length(ticks)) locs <- c(locs, which.min(abs(ticks[k] - leg.set)))
   axis(1, at=locs, labels=ticks, adj= 0.5, tick=T, cex=0.6, cex.axis=0.6, line=0, font=2, family="", mgp = c(0.1, 0.1, 0.1))

   # Plot in original order ------------------------------------------------------------------------------
   
   nf <- layout(matrix(c(1, 2), nrow=2, ncol=1, byrow=T), 1, c(12, 1), FALSE)
   right.space <-  ceiling(2 + 100/nrow(H))
   par(mar = c(3, 8, 4, right.space))
   V4 <- apply(V4, MARGIN=2, FUN=rev)
   image(1:nrow(V4), 1:ncol(V4), V4, zlim = c(0, max.cont.color + 3), col=mycol, axes=FALSE, main="Matrix in Original Order",
         sub = "", xlab= "", ylab="", cex.main=1)
   cex.rows <- 0.20 + 150/(nrow(V4) * max(nchar(row.names(V4))) + 200)         
   axis(3, at=1:nrow(V4), labels=row.names(V4), adj= 0.5, tick=FALSE, las = 1, cex.axis=cex.rows, font.axis=1, line=-1)
   cex.cols <- 0.20 + 200/(ncol(V4) * max(nchar(colnames(V4))) + 200)      
   mtext(colnames(V4), at=1:ncol(V4), side = 2, cex=cex.cols, col=cols, line=0, las=1, font=2, family="")

   # Legend

   par(mar = c(3, 5, 1, 25))   
   leg.set <- seq(-cutoff, cutoff, 0.05)
   image(1:length(leg.set), 1:1, as.matrix(leg.set), zlim=c(-cutoff, cutoff), col=mycol, axes=FALSE, main="Matrix Standardized Profile",
       sub = "", xlab= "", ylab="",font=2, family="", mgp = c(0, 0, 0), cex.main=0.8)
   ticks <- seq(-cutoff, cutoff, 0.5)
   tick.cols <- rep("black", 5)
   tick.lwd <- 1
   locs <- NULL
   for (k in 1:length(ticks)) locs <- c(locs, which.min(abs(ticks[k] - leg.set)))
   axis(1, at=locs, labels=ticks, adj= 0.5, tick=T, cex=0.6, cex.axis=0.6, line=0, font=2, family="", mgp = c(0.1, 0.1, 0.1))

   dev.off()
   
 }

   DISSECTOR_project_dataset.v1 <- function(
      #
      # Project a dataset in the space defined by a W matrix
      #
      input_dataset,                    # Input dataset (GCT)
      input_normalization = "rank",     # Normalization for the input dataset: "rank"
      normalize_after_match = T,        # Normalize input dataset after matching with rows of W
      input_W_dataset,                  # Input W matrix (GCT)
      W_normalization = "none",         # Normalization for W                                            
      output_H_dataset,                 # Output dataset H (GCT)
      output_W_dataset = NULL)          # Output dataset normalized W (GCT)                                            
  {

   set.seed(5209761)
       
   # Read input dataset

   dataset.1 <- MSIG.Gct2Frame(filename = input_dataset)
   m <- data.matrix(dataset.1$ds)
   print(dim(m))

   if (normalize_after_match == F) {  # Normalize input data here before matching with W
     if (input_normalization == "rank") {
         print("A before normalization:")
         print(m[1:5, ])

         max.n <- 10000
         for (i in 1:ncol(m)) m[,i] <- (max.n - 1) * (rank(m[,i]) - 1) /(nrow(m) - 1) + 1
 
         print("A after normalization:")
         print(m[1:5, ])
   
      } else if (input_normalization == "none") {   

      } else {
         stop(paste("ERROR: unknown input normalization:", input_normalization))  
      }
   }
   
   dataset.2 <- MSIG.Gct2Frame(filename = input_W_dataset)
   W <- data.matrix(dataset.2$ds)
   print(dim(W))

   k.comp <- ncol(W)

  # match gene lists

  overlap <- intersect(row.names(W), row.names(m))
  print(paste("overlap:", length(overlap)))
  locs1 <- match(overlap, row.names(m))
  locs2 <- match(overlap, row.names(W))
  m <- m[locs1,]
  W <- W[locs2,]

  ind <- order(row.names(W))
  W <- W[ind,]
  m <- m[ind,]

  if (normalize_after_match == T) {  # Normalize input data here after matching with W
     if (input_normalization == "rank") {
         print("A before normalization:")
         print(m[1:5, ])

         max.n <- 10000
         for (i in 1:ncol(m)) m[,i] <- (max.n - 1) * (rank(m[,i]) - 1) /(nrow(m) - 1) + 1
 
         print("A after normalization:")
         print(m[1:5, ])
   
      } else if (input_normalization == "none") {   

      } else {
         stop(paste("ERROR: unknown input normalization:", input_normalization))  
      }
   }

  # Normalize W

  if (W_normalization == "equal_sum") {
     norm.factors <- apply(W, MARGIN=2, FUN=sum)
     for (i in 1:ncol(W)) {
        W[, i] <- 1000*W[, i]/norm.factors[i]
      }
   }

   print("W:")
   print(W[1:5, ])
   print(" m:")
   print(m[1:5, ])
   
   H <- matrix(0, nrow=k.comp, ncol= ncol(m), dimnames=list(colnames(W), colnames(m)))
   for (i in 1:ncol(H)) H[, i] <- nnls.fit(W, m[, i], wsqrt=1, eps=0, rank.tol=1e-07)

   print("H:")
   print(H[1:5, ])
   
   # Save H matrix

   write.gct.2(gct.data.frame = H, descs = row.names(H), filename = output_H_dataset)

   # Save (optional) normalized W
   
   if(!(is.null(output_W_dataset))) write.gct.2(gct.data.frame = W, descs = row.names(W), filename = output_W_dataset)

 }

   DISSECTOR_generate_association_matrix.v1 <- function(
      #
      #  For an input file generate the association matrix between the columns (perturbation/states)
      #
      input_dataset,                    # Input dataset (GCT). This is e.g. an original dataset A or the H matrix
      input_dataset2 = NULL,            # Second input dataset (GCT). This is e.g. an original dataset A or the H matrix                                       
      association_type = "columns",     # Type of association: between "columns" or between "rows"
      exclude_suffix = F,               # Exclude suffix (sub-string after last "_") from row/column names
      annot_file = NULL,                # Column or row annotation file (TXT, optional) in format c(file, name_column, annot_column, use_prefix)
      annot_file2 = NULL,               # Column or row annotation file (TXT, optional) in format c(file, name_column, annot_column, use_prefix)      
      association_metric = "IC",        # Association metric: "IC" (Information Coefficient), "ICR" (IC ranked), "COR" (Pearson), "SPEAR" (Spearman)
      output_assoc_matrix_file,         # Output (GCT) file with association matrix
      output_assoc_plot)                # Output (PDF) file with association plot
   {

   set.seed(5209761)
   
   # Read input dataset

   dataset.1 <- MSIG.Gct2Frame(filename = input_dataset)
   H <- data.matrix(dataset.1$ds)
   print(paste("Dimensions dataset1:", dim(H)))
   if (association_type == "rows") H <- t(H)

   size <- ncol(H)/50
   if (size < 11) size <- 11
   pdf(file=output_assoc_plot, height=size, width=size)

   s <- strsplit(input_dataset, split="/")
   file.name <- s[[1]][length(s[[1]])]

   if (exclude_suffix == T) {
      row.names.H <- vector(length=nrow(H), mode="character")
      for (i in 1:nrow(H)) {
         temp <- unlist(strsplit(row.names(H)[i], split="_"))
         row.names.H[i] <- paste(temp[1:(length(temp)-1)], collaps="_")
     }
   } else {
      row.names.H <- row.names(H)
   }
 
   if (!is.null(input_dataset2)) {
      dataset.2 <- MSIG.Gct2Frame(filename = input_dataset2)
      H2 <- data.matrix(dataset.2$ds)
      print(paste("Dimensions dataset2:", dim(H2)))

     s <- strsplit(input_dataset2, split="/")
     file.name2 <- s[[1]][length(s[[1]])]

      if (association_type == "rows") H2 <- t(H2)

      if (exclude_suffix == T) {
         row.names.H2 <- vector(length=nrow(H2), mode="character")
         for (i in 1:nrow(H2)) {
            temp <- unlist(strsplit(row.names(H2)[i], split="_"))
            row.names.H2[i] <- paste(temp[1:(length(temp)-1)], collaps="_")
         }
      } else {
         row.names.H2 <- row.names(H2)
      }

   } else {
      H2 <- H
      row.names.H2 <- row.names.H
      file.name2 <- file.name
    }
   
   # Read annotation file

   if (!is.null(annot_file)) {
      annot.table <- read.table(annot_file[[1]], header=T, sep="\t", skip=0, colClasses = "character")
      gene.list <- annot.table[, annot_file[[2]]]
      annot.list <- annot.table[, annot_file[[3]]]
      gene.set <- vector(length=ncol(H), mode="character")
      if (annot_file[[4]] == T) {
         for (i in 1:ncol(H)) {
            gene.set[i] <- strsplit(colnames(H)[i], split="_")[[1]]
         }
      } else {
         gene.set <- colnames(H)
      }
      locs <- match(gene.set, gene.list)
      gene.class <- annot.list[locs]
      for (k in 1:length(gene.class)) gene.class[k] <- substr(gene.class[k], 1, 10)
      all.classes <- unique(gene.class)
      colnames(H) <- paste(colnames(H), " (", gene.class, ") ", sep="")
    }
   
   if (!is.null(annot_file2)) {
      annot.table2 <- read.table(annot_file2[[1]], header=T, sep="\t", skip=0, colClasses = "character")
      gene.list2 <- annot.table2[, annot_file2[[2]]]
      annot.list2 <- annot.table2[, annot_file2[[3]]]
      gene.set2 <- vector(length=ncol(H2), mode="character")
      if (annot_file2[[4]] == T) {
         for (i in 1:ncol(H2)) {
            gene.set2[i] <- strsplit(colnames(H2)[i], split="_")[[1]]
         }
      } else {
         gene.set2 <- colnames(H2)
      }
      locs <- match(gene.set2, gene.list2)
      gene.class2 <- annot.list2[locs]
      for (k in 1:length(gene.class2)) gene.class2[k] <- substr(gene.class2[k], 1, 10)
      all.classes2 <- unique(gene.class2)
      colnames(H2) <- paste(colnames(H2), " (", gene.class2, ") ", sep="")
    } else if (!is.null(annot_file)) {
       all.classes2 <- all.classes
       gene.class2 <- gene.class
    }

   # Define overlapping set

   overlap <- intersect(row.names.H, row.names.H2)
   print(paste("Size of overlap space:", length(overlap)))
   locs1 <- match(overlap, row.names.H)
   locs2 <- match(overlap, row.names.H2)   
   H <- H[locs1,]
   H2 <- H2[locs2,]

   # Signatures association plot
       
   nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(8, 1), FALSE)

   mycol <- vector(length=512, mode = "numeric")   # Red/Blue "pinkogram" color map
   for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol <- rev(mycol)
   ncolors <- length(mycol)

   col.classes <- c(brewer.pal(7, "Set1"), brewer.pal(7, "Pastel1"), brewer.pal(7, "Dark2"), brewer.pal(7, "Paired"), brewer.pal(7, "Pastel2"),
                 brewer.pal(8, "Accent"), brewer.pal(8, "Set2"), brewer.pal(11, "Spectral"), brewer.pal(12, "Set3"),
                 sample(c(brewer.pal(9, "Blues"), brewer.pal(9, "Reds"), brewer.pal(9, "Oranges"), brewer.pal(9, "Greys"),
                          brewer.pal(9, "Purples"), brewer.pal(9, "Greens"))))

   assoc.matrix <- matrix(0, nrow=ncol(H), ncol=ncol(H2), dimnames = list(colnames(H), colnames(H2)))
   for (i in 1:ncol(H)) {
      for (j in 1:ncol(H2)) {
         if (association_metric == "IC") {
            assoc.matrix[i, j] <- IC.v1(H[,i], H2[,j])
         } else if (association_metric == "ICR") {
            assoc.matrix[i, j] <- IC.v1(rank(H[,i]), rank(H2[,j]))
         } else if (association_metric == "COR") {
            assoc.matrix[i, j] <- cor(H[,i], H2[,j], method = "pearson")
         } else if (association_metric == "SPEAR") {
            assoc.matrix[i, j] <- cor(H[,i], H2[,j], method = "spearman")
         } else {
            stop(paste("ERROR: unknown association metric:", association_metric))
        }
       }
   } 
   dist.matrix <- as.dist(1 - assoc.matrix)

   hc <- hclust(dist(dist.matrix), "ward")
   hc2 <- hclust(dist(t(dist.matrix)), "ward")
       
#   HC <- hclust(dist.matrix, method="ward")
   assoc.matrix <- assoc.matrix[hc$order, hc2$order]

   write.gct.2(gct.data.frame = assoc.matrix, descs = row.names(assoc.matrix), filename = output_assoc_matrix_file)
 
   assoc.matrix <- ceiling(ncolors * (assoc.matrix + 1)/2)
   V <- apply(assoc.matrix, MARGIN=2, FUN=rev)
   par(mar = c(3, 8, 9, 4))
   image(1:dim(V)[2], 1:dim(V)[1], t(V), main = "", zlim = c(0, ncolors), col=mycol,
         axes=FALSE, xlab= "", ylab="") 

   mtext("Association Matrix", cex=1.3, side = 3, line = 7, outer=F)   
   mtext(file.name2, cex=1, side = 3, line = 5, outer=F)      
   mtext(file.name, cex=1, side = 2, line = 5, outer=F)

   
   if (!is.null(annot_file)) {
      cols <- col.classes[match(gene.class[hc$order], all.classes)]
   } else {
      cols <- "black"
   }
   if (!is.null(annot_file2)) {
      cols2 <- col.classes[match(gene.class2[hc2$order], all.classes2)]
   } else {
      cols2 <- cols
   }
   cex.rows <- 0.15 + 180/(nrow(V) * max(nchar(row.names(V))) + 200)
   cex.cols <-cex.rows
   mtext(row.names(V), at=1:nrow(V), side = 2, cex=cex.rows, col=rev(cols), line=0, las=1, font=2, family="")
   mtext(colnames(V), at=1:ncol(V), side = 3, cex=cex.cols, col=cols2, line=0, las=3, font=2, family="")

   # Legend

   par(mar = c(3, 25, 1, 5))
   leg.set <- seq(-1, 1, 0.01)
   image(1:length(leg.set), 1:1, as.matrix(leg.set), zlim=c(-1, 1), col=mycol, axes=FALSE, main=paste("Association Metric [", association_metric, "]", sep=""),
       sub = "", xlab= "", ylab="",font=2, family="", mgp = c(0, 0, 0), cex.main=0.8)
   ticks <- seq(-1, 1, 0.25)
   tick.cols <- rep("black", 5)
   tick.lwd <- c(1,1,2,1,1)
   locs <- NULL
   for (k in 1:length(ticks)) locs <- c(locs, which.min(abs(ticks[k] - leg.set)))
   axis(1, at=locs, labels=ticks, adj= 0.5, tick=T, cex=0.6, cex.axis=0.6, line=0, font=2, family="", mgp = c(0.1, 0.1, 0.1))
#   mtext(paste("Association Metric [", association_metric, "]", sep=""), cex=0.3, side = 1, line = 3.5, outer=F)
   
   dev.off()
 }

   DISSECTOR_extract_gene_sets.v1 <- function(
      #
      #  Extract/Define gene sets from an existing W matrix  (A ~ W x H)
      #
      input_dataset,                    # Input W matrix dataset (GCT)
      min.genes.per.comp = 50,          # Minimun number of genes per component gene set
      max.genes.per.comp = 150,         # Maximun number of genes per component gene set
      output_gene_sets_file,            # Output (GMT) file with gene sets
      output_plots)                     # Output (PDF) file with W and H plots
  {

   set.seed(5209761)
       
   pdf(file=output_plots, height=8.5, width=11)
    
   # Read input dataset

   dataset.1 <- MSIG.Gct2Frame(filename = input_dataset)
   W <- data.matrix(dataset.1$ds)
   print(dim(W))

   k.comp <- ncol(W)
    
   z.vec <- vector(length=length(k.comp), mode="numeric")
   for (i in 1:k.comp) {

      W.order <- order(W[,i], decreasing=T)
      W.vals <- W[W.order,i]
      W.genes <- row.names(W)[W.order]
      W.index <- seq(1, length(W.order))
      start <- ceiling(0.25*length(W.index))
      end <- ceiling(0.75*length(W.index))         
      x <- W.index[start:end]
      y <- W.vals[start:end]
      model <- lm(y ~ x)
      a <- model[[1]][1]
      b <- model[[1]][2]
      y.mod <- a + b*W.index
      env <- 0.20*mad(y)
      y.up <- y.mod + env
      y.dn <- y.mod - env
      dist <- (W.vals - y.up < 0)
      z  <- match("TRUE", dist)
      z.vec[i] <- z
      par(mar = c(10,10,10,10))
      plot(W.index, W.vals, pch=20, cex=0.8, col="gray", main=paste("Component ", colnames(W)[i], " Profile"),
           sub=paste("Component ", colnames(W)[i], "  cutoff z=", z, "(genes with rank less than z are included in gene set)"),
           xlab="Gene Rank", ylab="Component Amplitude (W column)")  # plot metagene amplitudes
      points(W.index[seq(1, z)], W.vals[seq(1, z)], pch=20, cex=0.8, col="black")
      points(W.index, y.up, pch=20, type="l", col="blue", lwd=1)  # plot metagene amplitudes
      points(W.index, y.mod, pch=20, type="l", col="blue", lwd=2)  # plot metagene amplitudes
      points(W.index, y.dn, pch=20, type="l", col="blue", lwd=1)  # plot metagene amplitudes            
      lines(c(z, z), range(W.vals), type="l", col = "red", lwd=2)
      lines(c(start, start), range(W.vals), type="l", col = "green", lwd=1)
      lines(c(end, end), range(W.vals), type="l", col = "green", lwd=1)            
    }
   z.vec <- ifelse(z.vec < min.genes.per.comp, min.genes.per.comp, z.vec)
   z.vec <- ifelse(z.vec > max.genes.per.comp, max.genes.per.comp, z.vec)

   # Produce gene sets (each with top.n genes) to represent each component

   gene.sets <- matrix(NA, nrow=k.comp, ncol=max(z.vec))
   for (i in 1:k.comp) {
      ind <- order(W[, i], decreasing=T)
      sorted.genes <- row.names(W)[ind]
      gene.sets[i, 1:z.vec[i]] <- sorted.genes[1:z.vec[i]]
    }

   # Save them as a "GMT" gene set file

   for (i in 1:k.comp) {
      row.header <- paste("C", i, "_", k.comp, sep="")
      output.line <- paste(gene.sets[i, 1:z.vec[i]], collapse="\t")
      output.line <- paste(row.header, row.header, output.line, sep="\t")
      if (i == 1) {
         write(noquote(output.line), file = output_gene_sets_file, append = F, ncolumns = length(z.vec[i]) + 2)
      } else {
        write(noquote(output.line), file = output_gene_sets_file, append = T, ncolumns = length(z.vec[i]) + 2)
      }
    }
    dev.off()
 }

   DISSECTOR_create_components.v1 <- function(
      #
      #  Project an input dataset into components using NMF                                            
      #
      input_dataset,                    # Input GCT dataset A where the matrix decomposition takes place (A ~ W x H)
      input_normalization = "rank",     # Normalization for the input dataset: "rank"
      number_of_comp,                   # Number of components to use in the matrix decomposition
      method = "NMF",                   # Method for matrix factorization: NMF or NMF_offset (IMF under construction)
      gene_subset = "all-genes",        # Universe of genes to consider for matrix decomposition: "gene-sets", "all-genes"
      gene_sets_files = NULL,           # If gene_subset = "gene-sets" GMT files with gene sets
      gene_sets = NULL,                 # If gene_subset = "gene-sets" then name of the specific gene set(s) in gene_sets_file to use
      normalize_after_selection = T,    # If gene_subset = "gene-sets," normalize after selection the gene subset
      preprojection_dataset = NULL,     # Save pre-projection input dataset in this file
      output_plots,                     # Output PDF file with W and H plots
      output_W_dataset,                 # Output GCT file with W matrix
      output_H_dataset,                 # Output GCT file with H matrix
      output_H_w_dataset)               # Output GCT file with W-derived H matrix                                              
  {

   set.seed(5209761)
    
   pdf(file=output_plots, height=8.5, width=11)

   comp.names <- paste("C", seq(1, number_of_comp), "_", number_of_comp, sep="")
    
   # Read expression dataset

   dataset.1 <- MSIG.Gct2Frame(filename = input_dataset)
   m.1 <- data.matrix(dataset.1$ds)
   print(paste("Dimensions matrix A:", nrow(m.1), ncol(m.1)))

   if (normalize_after_selection == F) {  # Normalize input data here before selection
     if (input_normalization == "rank") {
         print("A before normalization:")
         print(m.1[1:5, ])

         max.n <- 10000
         for (i in 1:ncol(m.1)) m.1[,i] <- (max.n - 1) * (rank(m.1[,i]) - 1) /(nrow(m.1) - 1) + 1
 
         print("A after normalization:")
         print(m.1[1:5, ])
   
      } else if (input_normalization == "none") {   

      } else {
         stop(paste("ERROR: unknown input normalization:", input_normalization))  
      }
   }
   
   if (gene_subset == "gene-sets") {  # select relevant genes from gene sets

	max.G <- 0
	max.N <- 0
	for (gsdb in gene_sets_files) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		max.G <- max(max.G, max(GSDB$size.G))
		max.N <- max.N +  GSDB$N.gs
	}
	N.gs <- 0
	gs <- matrix("null", nrow=max.N, ncol=max.G)
	gs.names <- vector(length=max.N, mode="character")
	gs.descs <- vector(length=max.N, mode="character")
	size.G <- vector(length=max.N, mode="numeric")
	start <- 1
	for (gsdb in gene_sets_files) {  # Read all the gene sets from gene set files
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		N.gs <- GSDB$N.gs 
		gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
		gs.descs[start:(start + N.gs - 1)] <- GSDB$gs.desc
		size.G[start:(start + N.gs - 1)] <- GSDB$size.G
		gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
		start <- start + N.gs
	}
	N.gs <- max.N
	
	# Select desired gene sets
	
	locs <- match(gene_sets, gs.names)
        print(rbind(gene_sets, locs))
	N.gs <- sum(!is.na(locs))
	if(N.gs > 1) { 
           gs <- gs[locs,]
	} else { 
           gs <- t(as.matrix(gs[locs,]))   # Force vector to matrix if only one gene set specified
        }
	gs.names <- gs.names[locs]
	gs.descs <- gs.descs[locs]
	size.G <- size.G[locs]

        genes <- NULL
       	for (gs.i in 1:N.gs) {
   	   gene.set <- gs[gs.i, 1:size.G[gs.i]]
           genes <- c(genes, gene.set)
         }
        print(paste("Number of selected genes:", length(genes)))
        genes <- unique(genes)
        print(paste("Number of unique selected genes:", length(genes)))        
        m.2 <- m.1[genes,]
        print("Dimensions of selected input data:")
        print(dim(m.2))

   } else if (gene_subset == "all-genes") {     
      m.2 <- m.1
   }

   if (normalize_after_selection == T) {  # Normalize input data here after selection
   
      if (input_normalization == "rank") {
         print("A before normalization:")
         print(m.2[1:5, ])

         max.n <- 10000
         for (i in 1:ncol(m.2)) m.2[,i] <- (max.n - 1) * (rank(m.2[,i]) - 1) /(nrow(m.2) - 1) + 1

         print("A after normalization:")
         print(m.2[1:5, ])

       } else if (input_normalization == "none") {   

       } else {
         stop(paste("ERROR: unknown input normalization:", input_normalization))  
       }
    }
   
   # Perform Matrix Factorization to find components

   if (!is.null(preprojection_dataset)) {
      write.gct.2(gct.data.frame = m.2, descs = row.names(m.2), preprojection_dataset)
    }

   if (method == "NMF") {
      NMF.out <- NMF.div(V = m.2, k = number_of_comp, maxniter = 1000, seed = 123, stopconv = 40, stopfreq = 10)
      W <- NMF.out$W
      H <- NMF.out$H
      row.names(W) <- row.names(m.2)
      colnames(W) <- row.names(H) <- comp.names
      colnames(H) <- colnames(m.2)

     plot(seq(1, length(NMF.out$error.v)), NMF.out$error.v, xlab="time", ylab="Error [divergence]", pch=20, col="blue")

      
   } else if (method == "IMF") {

    # To be added

   } else if (method == "NMF_offset") {

     library(NMF)
     NMF.out <- nmf(m.2, number_of_comp, "offset", seed=123)
     W <- basis(NMF.out)
     H <- coef(NMF.out)
     row.names(W) <- row.names(m.2)
     colnames(W) <- row.names(H) <- comp.names
     colnames(H) <- colnames(m.2)
   }
   
   # end of matrix factorization


  ind <- order(row.names(W))
  W <- W[ind,]
  m.2 <- m.2[ind,]
   
   # Obtain H via W: Project original and additional dataset using non-negative solver

   print("W:")
   print(W[1:5, ])
   print(" m:")
   print(m.2[1:5, ])
   
   H_w <- matrix(0, nrow=number_of_comp, ncol= ncol(m.2), dimnames=list(row.names(H), colnames(H)))
   for (i in 1:ncol(H_w)) H_w[, i] <- nnls.fit(W, m.2[, i], wsqrt=1, eps=0, rank.tol=1e-07)

   print("H_w:")
   print(H_w[1:5, ])

   # Save W and H matrices

   write.gct.2(gct.data.frame = W, descs = row.names(W), filename = output_W_dataset)
   write.gct.2(gct.data.frame = H, descs = row.names(H), filename = output_H_dataset)
   write.gct.2(gct.data.frame = H_w, descs = row.names(H_w), filename = output_H_w_dataset)   

   # Plot sorted W, H  and A matrices

    mycol <- vector(length=512, mode = "numeric")   # Red/Blue "pinkogram" color map
   for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol <- rev(mycol)
   ncolors <- length(mycol)

   hc <- hclust(dist(t(W)), "complete")
   d1.W <- as.dendrogram(hc)
   hc2 <- hclust(dist(W), "complete")
   d2.W <- as.dendrogram(hc2)
   heatmap(W, Colv=d1.W, Rowv = d2.W,  scale="row", col=mycol, margins=c(15, 15), cexRow=0.10, cexCol=0.5, main="Sorted W Matrix",
             xlab = "Components", ylab= "Genes")

   hc <- hclust(dist(t(H)), "complete")
   d1.H <- as.dendrogram(hc)
   heatmap(H, Colv=d1.H, Rowv = d1.W,  scale="col", col=mycol, margins=c(15, 15), cexRow=0.10, cexCol=0.5, main="Sorted H Matrix",
             xlab = "Components", ylab= "Genes")

   hc <- hclust(dist(t(H_w)), "complete")
   d1.H <- as.dendrogram(hc)
   heatmap(H_w, Colv=d1.H, Rowv = d1.W,  scale="col", col=mycol, margins=c(15, 15), cexRow=0.10, cexCol=0.5, main="Sorted H_w Matrix",
             xlab = "Components", ylab= "Genes")

   dev.off()
  
 }

IC.v1 <-  function(x, y, n.grid=25) {  # Information Coefficient

    # For definitions of mutual information and the universal metric (NMI) see the 
    # definition of "Mutual Information" in wikipedia and Thomas and Cover's book

   x.set <- !is.na(x)
   y.set <- !is.na(y)
   overlap <- x.set & y.set

   x <- x[overlap] +  0.000000001*runif(length(overlap))
   y <- y[overlap] +  0.000000001*runif(length(overlap))

   if (length(x) > 2) {
   
      delta = c(bcv(x), bcv(y))
   
      rho <- cor(x, y)
      rho2 <- abs(rho)
      delta <- delta*(1 + (-0.75)*rho2)
      kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
      FXY <- kde2d.xy$z + .Machine$double.eps
      dx <- kde2d.xy$x[2] - kde2d.xy$x[1]
      dy <- kde2d.xy$y[2] - kde2d.xy$y[1]
      PXY <- FXY/(sum(FXY)*dx*dy)
      PX <- rowSums(PXY)*dy
      PY <- colSums(PXY)*dx
      HXY <- -sum(PXY * log(PXY))*dx*dy
      HX <- -sum(PX * log(PX))*dx
      HY <- -sum(PY * log(PY))*dy
      PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
      PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
      MI <- sum(PXY * log(PXY/(PX*PY)))*dx*dy
      IC <- sign(rho) * sqrt(1 - exp(- 2 * MI))
      if (is.na(IC)) IC <- 0
   } else {
      IC <- 0
   }

   return(IC)
}

   DISSECTOR_extract_signature.v1 <- function(
      #
      # Create a gene set from a signature TXT file
      #
      input_dataset,                    # Input signature dataset (TXT) in format c(file, gene_column, use_prefix)
      n_markers = 200,                  # Number of markers genes to select
      gene_set_name,                    # Name given to the extracted gene set
      up_and_down = T,                  # Choose n.markers genes from tghe top and bottom of the list (F = only top)
      reference_dataset = NULL,         # Select only markers genes that overlap the gene list of another dataset (GCT)                                        
      output_dataset)                   # Output dataset (GMT)
  {
       
   # Read input dataset
  
   tab <- read.table(input_dataset[[1]], header=T, sep="\t", skip=0, colClasses = "character")
   gene.names <- as.vector(tab[, input_dataset[[2]] ])
   if (input_dataset[[3]] == T) {
      for (i in 1:length(gene.names)) {
         gene.names[i] <- unlist(strsplit(as.character(gene.names[i]), " "))[1]
      }
    }
   gene.names <- unique(gene.names)
   print(gene.names[1:10])  # print top 10
   print(gene.names[seq(length(gene.names) - 10 + 1, length(gene.names))]) # print bottom 10
                         
   # Read reference dataset

   if (!is.null(reference_dataset)) {                                                                          
      dataset.1 <- MSIG.Gct2Frame(filename = reference_dataset)
      m.1 <- data.matrix(dataset.1$ds)
      print(paste("Dataset dimensions:", dim(m.1)))
      overlap <- intersect(gene.names, row.names(m.1))
      locs <- match(overlap, gene.names)
      locs <- locs[!is.na(locs)]
      print(paste("Overlap with reference dataset:", length(locs)))
      gene.names <- gene.names[locs]
     }
       
   # Select top n.markers (and bottom) n.markers gene markers  (signature)

   if (up_and_down == T) { # UP and DOWN genes
      marker.set <- c(seq(1, n_markers), seq(length(gene.names) - n_markers + 1, length(gene.names)))
      gene.names <- gene.names[marker.set]
      print(paste("Marker set size:", length(gene.names)))
      print(gene.names[1:10])
      print(gene.names[seq(length(gene.names) - 10 + 1, length(gene.names))])
    } else {  # Only UP genes
      marker.set <- seq(1, n_markers)
      gene.names <- gene.names[marker.set]
      print(length(gene,names))      
      print(gene.names[1:10])
    }

  # Save them in a GMT file
   
      row.header <- gene_set_name
      output.line <- paste(gene.names, collapse="\t")
      output.line <- paste(row.header, row.header, output.line, sep="\t")
      write(noquote(output.line), file = output_dataset, append = F, ncolumns = length(gene.names) + 2)

 }

                                                                          
Read.GeneSets.db <- function(
   gs.db,
   thres.min = 2,
   thres.max = 2000,
   gene.names = NULL)
  {

   temp <- readLines(gs.db)
   max.Ng <- length(temp)
   temp.size.G <- vector(length = max.Ng, mode = "numeric") 
   for (i in 1:max.Ng) {
      temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
   }
   max.size.G <- max(temp.size.G)      
   gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
   temp.names <- vector(length = max.Ng, mode = "character")
   temp.desc <- vector(length = max.Ng, mode = "character")
   gs.count <- 1
   for (i in 1:max.Ng) {
      gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
      gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
      gene.set.name <- gs.line[1] 
      gene.set.desc <- gs.line[2] 
      gene.set.tags <- vector(length = gene.set.size, mode = "character")
      for (j in 1:gene.set.size) {
         gene.set.tags[j] <- gs.line[j + 2]
      }
      if (is.null(gene.names)) {
         existing.set <- rep(TRUE, length(gene.set.tags))
      } else {
         existing.set <- is.element(gene.set.tags, gene.names)
      }
      set.size <- length(existing.set[existing.set == T])
      if ((set.size < thres.min) || (set.size > thres.max)) next
      temp.size.G[gs.count] <- set.size
      gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
      temp.names[gs.count] <- gene.set.name
      temp.desc[gs.count] <- gene.set.desc
      gs.count <- gs.count + 1
    }
   Ng <- gs.count - 1
   gs.names <- vector(length = Ng, mode = "character")
   gs.desc <- vector(length = Ng, mode = "character")
   size.G <- vector(length = Ng, mode = "numeric") 
   
   gs.names <- temp.names[1:Ng]
   gs.desc <- temp.desc[1:Ng]
   size.G <- temp.size.G[1:Ng]
   
   return(list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, size.G = size.G, max.N.gs = max.Ng))
 }


## ------------- Additional Functions from CNMF.v4 ------------------------


consensusNMF.2 <- function(input.ds, k.init, k.final, num.clusterings, maxniter, error.function, rseed=123456789, directory = "", stopconv = 40, stopfreq = 10, non.interactive.run = F, doc.string = "", ...) {

#
#  GenePattern Methodology for:
#
#  Metagenes and Molecular Pattern Discovery using Matrix Factorization
#  Jean-Philippe Brunet, Pablo Tamayo, Todd. R. Golub, and Jill P. Mesirov
# 
#  Author:  Pablo Tamayo (tamayo@genome.wi.mit.edu)
#
#  Based on the original matlab version written by Jean-Philippe Brunet (brunet@broad.mit.edu) and
#  with additional contributions from: Ted Liefeld (liefeld@broad.mit.edu)   
#  Date:  November 27, 2003
#
#  Last change March 3, 2005: modifications to make the output more readable.
#
#  Execute from an R console window with this command:
#  source("<this file>", echo = TRUE)
#  E.g. someoutput <- mynmf2(input.ds="c:\\nmf\\all_aml.res",k.init=2,k.final=5,num.clusterings=20,maxniter=500) 
#
#  For details on the method see:
#
#  Proc. Natl. Acad. Sci. USA 2004 101: 4164-4169
#  http://www.broad.mit.edu/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=89
#
#  Input parameters
#
#   input.ds
#                       input gene expression dataset in GCT or RES format
#   k.init
#                       initial value of k
#   k.final
#                       final value of k
#   num.clusterings
#                       number of NMF clusterings to build consensus matrix
#   maxniter
#                       maximum number of NMF iterations
#   error.function
#                       NMF error function: "divergence" of "euclidean"
#   rseed
#                       random number generator seed
#   directory
#                       file directory where to store the result files
#   stopconv
#                       how many no change checks are needed to stop NMF iterations (convergence)
#   stopfreq
#                       frequency (NMF iterations) of "no change" checks 
#   non.interactive.run 
#                       flag controlling if the plots are produced interatively (Rgui and saved) or only saved in files
#   doc.string
#                       prefix to be added to the output files
#
#  Output files are (prefix with "doc.string")
#
#   params.txt 
#                       run parameters and time of execution
#   membership.gct		
#			membership results for samples at all values of K
#   cophenetic.txt 
#			cophenetic values for each K
#   cophenetic.plot.jpeg
#			plot of cophenetic for each value of K		
#   consensus.k#.gct (for each value of K)
#			consensus matrix for k=#
#   consensus.plot.k#.jpeg (for each value of K)
#			plot of consensus matrix for k=#
#   graphs.k#.jpeg (for each value of K)

# save input parameters

filename <- paste(directory, doc.string, ".params.txt", sep="", collapse="")  

time.string <- as.character(as.POSIXlt(Sys.time(),"GMT"))
write(paste("Run of NMF on ", time.string), file=filename)

write(paste("input.ds =", input.ds, sep=" "), file=filename, append=T) 
write(paste("k.init = ", k.init, sep=" "), file=filename, append=T) 
write(paste("k.final =", k.final, sep=" "), file=filename, append=T) 
write(paste("num.clusterings =", num.clusterings, sep=" "), file=filename, append=T) 
write(paste("maxniter =", maxniter, sep=" "), file=filename, append=T) 
write(paste("error.function =", error.function, sep=" "), file=filename, append=T) 
write(paste("rseed =", rseed, sep=" "), file=filename, append=T) 
write(paste("directory =", directory, sep=" "), file=filename, append=T) 
write(paste("stopconv =", stopconv, sep=" "), file=filename, append=T) 
write(paste("stopfreq =", stopfreq, sep=" "), file=filename, append=T)
write(paste("non.interctive.run =", non.interactive.run, sep=" "), file=filename, append=T) 
write(paste("doc.string =", doc.string, sep=" "), file=filename, append=T) 


k.init<-as.integer(k.init)
k.final<-as.integer(k.final)
num.clusterings<-as.integer(num.clusterings)
n.iter<-as.integer(maxniter)
if (!is.na(rseed)){
     seed <- as.integer(rseed)
}


# library(mva)
# library(MASS)
# library(GenePattern)

D <- CNMF.read.dataset(input.ds)
A <- data.matrix(D)

# Threshold negative values to small quantity 

eps <- .Machine$double.eps
A[A < 0] <- eps



cols <- length(A[1,])
rows <- length(A[,1])

col.names <- names(D)

num.k <- k.final - k.init + 1

rho <- vector(mode = "numeric", length = num.k)
k.vector <- vector(mode = "numeric", length = num.k)

k.index <- 1

connect.matrix.ordered <- array(0, c(num.k, cols, cols))

filename <- paste(directory, doc.string, ".", "graphs.pdf", sep="", collapse="")
pdf(file=filename, width = 9, height = 11)

for (k in k.init:k.final) { 

   nf <- layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow=T), c(1, 1), c(1, 1, 1, 1), TRUE)
   assign <- matrix(0, nrow = num.clusterings, ncol = cols)

   for (i in 1:num.clusterings) {
	  
        print(paste("Computing clustering number=", i, " for k=", k, sep=""))

        if (error.function == "divergence"){
	    NMF.out <- NMF.div(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
	} else if (error.function == "euclidean"){
	    NMF.out <- NMF(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
	} else {
            stop(paste("Un-supported error function=", error.function, sep=""))
        }
        print(paste(NMF.out$t, " NMF iterations performed", sep=""))

        for (j in 1:cols) { # Find membership
            class <- order(NMF.out$H[,j], decreasing=T)
            assign[i, j] <- class[1]
        }

	if (i == 1) {  # Plot example for first clustering iteration
            H.saved <- NMF.out$H
            sub.string <- paste(doc.string, " k=", k, sep="")
            plot(1:NMF.out$t, NMF.out$error.v[1:NMF.out$t], pch = 20, cex = 1.5, col = 1, xlab="time", ylab="NMF error", sub=sub.string, main=paste("Convergence plot k=", k, " example", sep=""))


            if (rows < 1000) {
               W <- NMF.out$W
            } else {
               W <- NMF.out$W[sample(x = 1:rows, size = 1000),]
            }
            sub.string <- paste(doc.string, " k=", k, sep="")
            CNMF.matrix.abs.plot(W, sub = sub.string, log = F, main = "Example W matrix (orig. order)", ylab = "genes", xlab ="metasamples")
            CNMF.matrix.abs.plot(H.saved, sub = sub.string, log = F, main = "Example H matrix (orig. order)", ylab = "metagenes", xlab ="samples")
            CNMF.metagene.plot(H = H.saved, main = "Metagenes Example (orig. order)", sub = sub.string, xlab = "samples", ylab = "metagenes")

        }

        rm(NMF.out)

     }  ## end  for (i in 1:num.clusterings)

   
     # compute consensus matrix
     connect.matrix <- matrix(0, nrow = cols, ncol = cols)

     for (i in 1:num.clusterings) {
       for (j in 1:cols) {
          for (p in 1:cols) {
             if (j != p) {
                  if (assign[i, j] == assign[i, p]) {
                    connect.matrix[j, p] <- connect.matrix[j, p] + 1
                  } 
              } else {
                    connect.matrix[j, p] <- connect.matrix[j, p] + 1
              }
           }
       }
     }

     connect.matrix <- connect.matrix / num.clusterings

     dist.matrix <- 1 - connect.matrix
     dist.matrix <- as.dist(dist.matrix)
     HC <- hclust(dist.matrix, method="average")

     dist.coph <- cophenetic(HC)
     k.vector[k.index] <- k
     rho[k.index] <- cor(dist.matrix, dist.coph)
     rho[k.index] <- signif(rho[k.index], digits = 4)
   
#     connect.matrix.ordered <- matrix(0, nrow=cols, ncol = cols)

     for (i in 1:cols) {
        for (j in 1:cols) {
           connect.matrix.ordered[k.index, i, j] <- connect.matrix[HC$order[i], HC$order[j]]
         }
     }

     # compute consensus clustering membership

     membership <- cutree(HC, k = k)

     max.k <- max(membership)
     items.names.ordered <- col.names[HC$order]
     membership.ordered <- membership[HC$order]
     results <- data.frame(cbind(membership.ordered, items.names.ordered))

     if (k > k.init){
          all.membership <- cbind(all.membership, membership);
     } else {
          all.membership <- cbind(membership);
     }

     sub.string <- paste(doc.string, " k=", k, sep="")
     CNMF.matrix.abs.plot(connect.matrix.ordered[k.index,,], sub=sub.string, log = F, main = "Ordered Consensus Matrix", 
                          ylab = "samples", xlab ="samples")
     plot(HC, xlab="samples", cex = 0.75, labels = col.names, sub = sub.string, col = "blue", 
          main = paste("Ordered Linkage Tree. Coph=", rho[k.index]))

     matrixGct <- data.frame(connect.matrix.ordered[k.index,,])
     filename <- paste(directory, doc.string, ".", "matrix.k.",k, ".gct", sep="", collapse="")
     CNMF.write.gct.2(matrixGct, descs = "", filename) 

     resultsGct <- data.frame(membership.ordered)
     row.names(resultsGct) <- items.names.ordered
     filename <- paste(directory, doc.string, ".", "consensus.k.",k, ".gct", sep="", collapse="")
     CNMF.write.gct.2(resultsGct, descs = "", filename) 

     H.sorted <- H.saved[,HC$order]
     sub.string <- paste(doc.string, " k=", k, sep="")
     CNMF.matrix.abs.plot(H.sorted, sub = sub.string, log = F, main = "Example H matrix (ordered)", ylab = "metagenes", xlab ="samples")
     CNMF.metagene.plot(H = H.sorted, sub = sub.string, main = "Metagenes Example (ordered)", xlab = "samples", ylab = "metagenes")

     nf <- layout(matrix(c(1), 1, 1, byrow=T), c(1, 1), c(1, 1), TRUE)

     conlabel <- paste("Consensus k =", k, sep=" ", collapse="")

     sub.string <- paste("Consensus matrix k=", k, "; dataset= ", input.ds, sep="")
     CNMF.ConsPlot(connect.matrix.ordered[k.index,,], col.labels = membership.ordered, col.names = items.names.ordered, 
                   main = " ", sub=sub.string, xlab=" ", ylab=" ")

     k.index <- k.index + 1

} # end of loop over k


# Save consensus matrices in one file

  nf <- layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), 4, 4, byrow=T), c(1, 1, 1, 1), c(1, 1, 1, 1), TRUE)

  for (k in 1:num.k) { 
     CNMF.matrix.abs.plot(connect.matrix.ordered[k,,], log = F, main = paste("k=", k.vector[k]), 
                          sub = paste("Cophenetic coef.=", rho[k]), ylab = "samples", xlab ="samples")
  }
   
  y.range <- c(1 - 2*(1 - min(rho)), 1)
  plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, 
                          xlab = "k", ylab="Cophenetic correlation", type = "n")
  lines(k.vector, rho, type = "l", col = "black")
  points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")

# Write the membership matrix

resultsmembership <- data.frame(all.membership)
row.names(resultsmembership) <- col.names
colnames(resultsmembership) <- paste("k=", seq(k.init, k.final))

print("Membership:")

print(resultsmembership)

filename <- paste(directory, doc.string, ".", "membership", ".txt", sep="", collapse="")

col.names <- paste(colnames(resultsmembership), collapse = "\t")
col.names <- paste("SAMPLE", col.names, sep= "\t")
write(noquote(col.names), file = filename, append = F, ncolumns = length(col.names))
write.table(resultsmembership, file=filename, quote=F, col.names = F, row.names = T, append = T, sep="\t")

  nf <- layout(matrix(c(1), 1, 1, byrow=T), c(1), c(1), TRUE)
y.range <- c(1 - 2*(1 - min(rho)), 1)
plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, xlab = "k", ylab="Cophenetic correlation", type = "n")
lines(k.vector, rho, type = "l", col = "black")
points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")

filename <- paste(directory, doc.string, ".", "cophenetic.txt", sep="")

xx <- cbind(k.vector, rho)
write(noquote(c("k", "\t", "Cophenetic Coefficient")), file = filename, append = F, ncolumns = 1000)
write.table(xx, file = filename, append = T, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = F)

dev.off()
}


CNMF.read.dataset <- function(file) {
	result <- regexpr(paste(".gct","$",sep=""), tolower(file))
	if(result[[1]] != -1)
		return(CNMF.read.gct(file))
	result <- regexpr(paste(".res","$",sep=""), tolower(file))
	if(result[[1]] != -1)
		return(CNMF.read.res(file))
	stop("Input is not a res or gct file.")	
}

CNMF.matrix.abs.plot <- function(V, axes = F, log = F, norm = T, transpose = T, matrix.order = T, max.v = 1, min.v = 0, main = " ", sub = " ", xlab = " ", ylab = "  ") {
      rows <- length(V[,1])
      cols <- length(V[1,])
      if (log == T) {
         V <- log(V)
      }
      B <- matrix(0, nrow=rows, ncol=cols)
	for (i in 1:rows) {
           for (j in 1:cols) {
                if (matrix.order == T) {
                   k <- rows - i + 1
                } else {
                   k <- i
                }
                if (norm == T) {
                  if ((max.v == 1) && (min.v == 0)) {
                     max.val <- max(V)
                     min.val <- min(V)
                  } else {
		     	   max.val = max.v
                     min.val = min.v
                  }
               }
	     B[k, j] <-  max.val - V[i, j] + min.val
           }
      }
	if (transpose == T) {
	  B <- t(B)
        }
	if (norm == T) {
#            image(z = B, zlim = c(min.val, max.val), axes = axes, col = rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), main = main, sub = sub, xlab = xlab, ylab = ylab)
            image(z = B, zlim = c(min.val, max.val), axes = axes, col = rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75), main = main, sub = sub, xlab = xlab, ylab = ylab)           
      } else {
#            image(z = B, axes = axes, col = rainbow(100, s = 1, v = 0.6, start = 0.1, end = 0.9, gamma = 1), main = main, sub = sub, xlab = xlab, ylab = ylab)
            image(z = B, axes = axes, col = rainbow(100, s = 1, v = 0.6, start = 0.1, end = 0.9), main = main, sub = sub, xlab = xlab, ylab = ylab)         
      }
      return(list(B, max.val, min.val))
    }

CNMF.metagene.plot <- function(H, main = " ", sub = " ", xlab = "samples ", ylab = "amplitude") {
	k <- length(H[,1])
	S <- length(H[1,])
	index <- 1:S
	maxval <- max(H)
        minval <- min(H)
	plot(index, H[1,], xlim=c(1, S), ylim=c(minval, maxval), main = main, sub = sub, ylab = ylab, xlab = xlab, type="n")
	for (i in 1:k) {
	    lines(index, H[i,], type="l", col = i, lwd=2)
        }
}


CNMF.ConsPlot <- function(V, col.labels, col.names, main = " ", sub = " ", xlab=" ", ylab=" ") {

# Plots a heatmap plot of a consensus matrix

     cols <- length(V[1,])
     B <- matrix(0, nrow=cols, ncol=cols)
     max.val <- max(V)
     min.val <- min(V)
     for (i in 1:cols) {
         for (j in 1:cols) {
             k <- cols - i + 1
	     B[k, j] <-  max.val - V[i, j] + min.val
          }
     }

     col.names2 <- rev(col.names)
     col.labels2 <- rev(col.labels)
     D <- matrix(0, nrow=(cols + 1), ncol=(cols + 1))

     col.tag <- vector(length=cols, mode="numeric")
     current.tag <- 0
     col.tag[1] <- current.tag
     for (i in 2:cols) {
        if (col.labels[i] != col.labels[i - 1]) {
             current.tag <- 1 - current.tag
        }
        col.tag[i] <- current.tag
     }
     col.tag2 <- rev(col.tag)
     D[(cols + 1), 2:(cols + 1)] <- ifelse(col.tag %% 2 == 0, 1.02, 1.01)
     D[1:cols, 1] <- ifelse(col.tag2 %% 2 == 0, 1.02, 1.01)
     D[(cols + 1), 1] <- 1.03
     D[1:cols, 2:(cols + 1)] <- B[1:cols, 1:cols]

#     col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), "#BBBBBB", "#333333", "#FFFFFF")
     col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75), "#BBBBBB", "#333333", "#FFFFFF")     
     image(1:(cols + 1), 1:(cols + 1), t(D), col = col.map, axes=FALSE, main=main, sub=sub, xlab= xlab, ylab=ylab)
     for (i in 1:cols) {
         col.names[i]  <- paste("      ", substr(col.names[i], 1, 12), sep="")
         col.names2[i] <- paste(substr(col.names2[i], 1, 12), "     ", sep="")
     }

     axis(2, at=1:cols, labels=col.names2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.50, font.axis=1, line=-1)
     axis(2, at=1:cols, labels=col.labels2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.65, font.axis=1, line=-1)

     axis(3, at=2:(cols + 1), labels=col.names, adj= 1, tick=FALSE, las = 3, cex.axis=0.50, font.axis=1, line=-1)
     axis(3, at=2:(cols + 1), labels=as.character(col.labels), adj = 1, tick=FALSE, las = 1, cex.axis=0.65, font.axis=1, line=-1)

     return()
   }

CNMF.read.res <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in RES format and converts it into an R data frame
#
   header.cont <- readLines(filename, n = 1)
   temp <- unlist(strsplit(header.cont, "\t"))
   colst <- length(temp)
   header.labels <- temp[seq(3, colst, 2)]
   ds <- read.delim(filename, header=F, row.names = 2, sep="\t", skip=3, blank.lines.skip=T, comment.char="", as.is=T)
   colst <- length(ds[1,])
   cols <- (colst - 1)/2
   rows <- length(ds[,1])
   A <- matrix(nrow=rows - 1, ncol=cols)
   A <- ds[1:rows, seq(2, colst, 2)]
   table1 <- data.frame(A)
   names(table1) <- header.labels
   return(table1)
}

CNMF.read.gct <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
   ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T)
   ds <- ds[-1]
   return(ds)
}

CNMF.write.gct.2 <- function(gct.data.frame, descs = "", filename) 
{
    f <- file(filename, "w")
    cat("#1.2", "\n", file = f, append = TRUE, sep = "")
    cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
    cat("Name", "\t", file = f, append = TRUE, sep = "")
    cat("Description", file = f, append = TRUE, sep = "")

    colnames <- colnames(gct.data.frame)
    cat("\t", colnames[1], file = f, append = TRUE, sep = "")

    if (length(colnames) > 1) {
       for (j in 2:length(colnames)) {
           cat("\t", colnames[j], file = f, append = TRUE, sep = "")
       }
     }
    cat("\n", file = f, append = TRUE, sep = "\t")

    oldWarn <- options(warn = -1)
    m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
    m[, 1] <- row.names(gct.data.frame)
    if (length(descs) > 1) {
        m[, 2] <- descs
    } else {
        m[, 2] <- row.names(gct.data.frame)
    }
    index <- 3
    for (i in 1:dim(gct.data.frame)[2]) {
        m[, index] <- gct.data.frame[, i]
        index <- index + 1
    }
    write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
    close(f)
    options(warn = 0)

}

NMF.div <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {

        N <- length(V[,1])
        M <- length(V[1,])
        set.seed(seed)
        W <- matrix(runif(N*k), nrow = N, ncol = k)  # Initialize W and H with random numbers
        H <- matrix(runif(k*M), nrow = k, ncol = M)
        VP <- matrix(nrow = N, ncol = M)
        error.v <- vector(mode = "numeric", length = maxniter)
        new.membership <- vector(mode = "numeric", length = M)
        old.membership <- vector(mode = "numeric", length = M)
        no.change.count <- 0
        eps <- .Machine$double.eps
        for (t in 1:maxniter) {
                VP = W %*% H
                W.t <- t(W)
                H <- H * (W.t %*% (V/VP)) + eps
                norm <- apply(W, MARGIN=2, FUN=sum)
                for (i in 1:k) {
                    H[i,] <- H[i,]/norm[i]
                }
                VP = W %*% H
                H.t <- t(H)
                W <- W * ((V/(VP + eps)) %*% H.t) + eps
                norm <- apply(H, MARGIN=1, FUN=sum)
                for (i in 1:k) {
                    W[,i] <- W[,i]/norm[i]
                }
               error.v[t] <- sum(V * log((V + eps)/(VP + eps)) - V + VP)/(M * N)
               if (t %% stopfreq == 0) {

                    for (j in 1:M) {
                        class <- order(H[,j], decreasing=T)
                        new.membership[j] <- class[1]
                     }
                     if (sum(new.membership == old.membership) == M) {
                        no.change.count <- no.change.count + 1
                     } else {
                        no.change.count <- 0
                     }
                     if (no.change.count == stopconv) break
                     old.membership <- new.membership
               }
        }
        return(list(W = W, H = H, t = t, error.v = error.v))
}

MSIG.Preprocess.Dataset <- function(
   input.ds, 
   output.ds,
   thres = NULL, 
   ceil = NULL, 
   shift = NULL,
   fold = NULL, 
   delta = NULL, 
   normalization = NULL,
   cntrl.genes = NULL) {

   print(c("Running MSIG.Preprocess.Dataset... on:", input.ds))
   print(c("output file:", output.ds))
   print(c("normalization =", normalization))
   
# Read dataset

   dataset <- MSIG.Gct2Frame(filename = input.ds)
   m <- data.matrix(dataset$ds)
   gs.names <- dataset$row.names
   gs.descs <- dataset$descs
   sample.names <- dataset$names

# threshold, ceiling and shift

   if (!is.null(thres)) {
     m[m < thres] <- thres
   }
   if (!is.null(ceil)) {
      m[m > ceil] <- ceil
   }
   if (!is.null(shift)) {
      m <- m + shift
   }

   # identify and save control genes

   if (!is.null(cntrl.genes)) {
      gene.names2 <- intersect(cntrl.genes, gs.names)
      locs <- match(gene.names2, gs.names, nomatch=0)
      msig.cntrl <- m[locs, ]
      msig.cntrl.genes <- gs.names[locs]
      msig.cntrl.descs <- gs.descs[locs]
      m <- m[-locs, ]
      gs.names <- gs.names[-locs]
      gs.descs <- gs.descs[-locs]
    }

   # variation filter

   if ((!is.null(fold)) && (!is.null(delta))) {
      temp <- MSIG.VarFilter(V = m, fold = fold, delta = delta, gene.names = gs.names, gene.descs = gs.descs) 
      m <- temp$V
      gs.names <- temp$new.gene.names
      gs.descs <- temp$new.gene.descs
      dim(m) 
   }

   # restore control genes

   if (!is.null(cntrl.genes)) {
      m <- rbind(m, msig.cntrl)
      gs.names <- c(gs.names, msig.cntrl.genes)
      gs.descs <- c(gs.descs, msig.cntrl.descs)
    }

# normalization

   if (!is.null(normalization)) {
      if (normalization == 1) {
         m <- MSIG.NormalizeCols.Rank(m)
      } else if (normalization == 2) {
         m <- MSIG.NormalizeCols.Rank(m)/length(m[,1])
      } else if (normalization == 3) {
         m <- GSEA.NormalizeCols(m) + 3
         m <- GSEA.Threshold(m, 0.001, 100000) 
      } else if (normalization == 4) {
         m <- MSIG.NormalizeCols.Rank(m)/length(m[,1])
      } else if (normalization == 5) {
         m <- MSIG.NormalizeCols.Rescale(m)
      } else if (normalization == 6) {
         cols <- length(m[1,])
         for (j in 1:cols) {  # column rank normalization from 0 to N - 1
            m[,j] <- rank(m[,j], ties.method = "average") - 1
         }
         m <- 10000*m/(length(m[,1]) - 1)
      } else if (normalization == 7) {
         m <- ((100*MSIG.NormalizeCols.Rank(m))%/%length(m[,1]) + 1)
      } else if (normalization == 8) { 
          row.mean <- apply(m, MARGIN=1, FUN=mean)
          for (i in 1:length(m[,1])) {
             m[i,] <- m[i,] / row.mean[i]
          }
      }
   }
   
   V <- data.frame(m)
   names(V) <- sample.names
   row.names(V) <- gs.names
   write.gct(gct.data.frame = V, descs = gs.descs, filename = output.ds)  

 }

GSEA.Threshold <- function(V, thres, ceil) { 
#
# Threshold and ceiling pre-processing for gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        V[V < thres] <- thres
        V[V > ceil] <- ceil
        return(V)
}

GSEA.VarFilter <- function(V, fold, delta, gene.names = "") { 
#
# Variation filter pre-processing for gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        cols <- length(V[1,])
        rows <- length(V[,1])
        row.max <- apply(V, MARGIN=1, FUN=max)
               row.min <- apply(V, MARGIN=1, FUN=min)
        flag <- array(dim=rows)
        flag <- (row.max /row.min >= fold) & (row.max - row.min >= delta)
        size <- sum(flag)
        B <- matrix(0, nrow = size, ncol = cols)
        j <- 1
        if (length(gene.names) == 1) {
           for (i in 1:rows) {
              if (flag[i]) {
                 B[j,] <- V[i,]
                 j <- j + 1
               }
           }
        return(B)
        } else {
            new.list <- vector(mode = "character", length = size)
            for (i in 1:rows) {
              if (flag[i]) {
                 B[j,] <- V[i,]
                 new.list[j] <- gene.names[i]
                 j <- j + 1
              }
            }
        return(list(V = B, new.list = new.list))
        }
}

MSIG.VarFilter <- function(V, fold, delta, gene.names = "", gene.descs = "") { 

# Variation filter pre-processing for gene expression matrix

        cols <- length(V[1,])
        rows <- length(V[,1])
        row.max <- apply(V, MARGIN=1, FUN=max)
        row.min <- apply(V, MARGIN=1, FUN=min)
        flag <- array(dim=rows)
        flag <- (row.max /row.min >= fold) & (row.max - row.min >= delta)
        size <- sum(flag)
        B <- matrix(0, nrow = size, ncol = cols)
        j <- 1
        if (length(gene.names) == 1) {
           for (i in 1:rows) {
              if (flag[i]) {
                 B[j,] <- V[i,]
                 j <- j + 1
               }
           }
        return(B)
        } else {
            new.gene.names <- vector(mode = "character", length = size)
            new.gene.descs <- vector(mode = "character", length = size)
            for (i in 1:rows) {
              if (flag[i]) {
                 B[j,] <- V[i,]
                 new.gene.names[j] <- gene.names[i]
                 new.gene.descs[j] <- gene.descs[i]
                 j <- j + 1
              }
            }
        return(list(V = B, new.gene.names = new.gene.names, new.gene.descs = new.gene.descs, locations = flag))
        }
}

GSEA.NormalizeRows <- function(V) { 
#
# Stardardize rows of a gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        row.mean <- apply(V, MARGIN=1, FUN=mean)
        row.sd <- apply(V, MARGIN=1, FUN=sd)

        row.n <- length(V[,1])
        for (i in 1:row.n) {
             if (row.sd[i] == 0) {
                  V[i,] <- 0
           } else {
              V[i,] <- (V[i,] - row.mean[i])/row.sd[i]
           }
        }
        return(V)
}

GSEA.NormalizeCols <- function(V) { 
#
# Stardardize columns of a gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        col.mean <- apply(V, MARGIN=2, FUN=mean)
               col.sd <- apply(V, MARGIN=2, FUN=sd)
        col.n <- length(V[1,])
        for (i in 1:col.n) {
             if (col.sd[i] == 0) {
                  V[i,] <- 0
           } else {
              V[,i] <- (V[,i] - col.mean[i])/col.sd[i]
           }
        }
        return(V)
}

GSEA.NormalizeCols.Rank <- function(V) { 
#
      cols <- length(V[1,])
      rows <- length(V[,1])
      for (j in 1:cols) {  # column rank normalization
         V[,j] <- rank(V[,j], ties.method = "average")
      }

      return(V)
}


MSIG.NormalizeCols.Rank <- function(V) { 

      cols <- length(V[1,])
      rows <- length(V[,1])
      for (j in 1:cols) {  # column rank normalization
         V[,j] <- rank(V[,j], ties.method = "average")
      }

      return(V)
}

MSIG.NormalizeCols.Rescale <- function(V) { 

      epsilon <- 0.00001
      cols <- length(V[1,])
      for (j in 1:cols) {  # column rank normalization
         max.v <- max(V[,j])
         min.v <- min(V[,j])
         V[,j] <- (V[,j] - min.v + epsilon)/(max.v - min.v)
      }

      return(V)
    }

MSIG.Gct2Frame <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

   ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, na.strings = "")
   descs <- ds[,1]
   ds <- ds[-1]
   row.names <- row.names(ds)
   names <- names(ds)
   return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}


write.gct <- function(gct.data.frame, descs = "", filename) 
{
    f <- file(filename, "w")
    cat("#1.2", "\n", file = f, append = TRUE, sep = "")
    cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
    cat("Name", "\t", file = f, append = TRUE, sep = "")
    cat("Description", file = f, append = TRUE, sep = "")

    names <- names(gct.data.frame)
    cat("\t", names[1], file = f, append = TRUE, sep = "")

    if (length(names) > 1) {
       for (j in 2:length(names)) {
           cat("\t", names[j], file = f, append = TRUE, sep = "")
       }
     }
    cat("\n", file = f, append = TRUE, sep = "\t")

    oldWarn <- options(warn = -1)
    m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
    m[, 1] <- row.names(gct.data.frame)
    if (length(descs) > 1) {
        m[, 2] <- descs
    } else {
        m[, 2] <- row.names(gct.data.frame)
    }
    index <- 3
    for (i in 1:dim(gct.data.frame)[2]) {
        m[, index] <- gct.data.frame[, i]
        index <- index + 1
    }
    write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
    close(f)
    options(warn = 0)

  }

cluster.and.resort <- function(V, rows, cols) {
        dist.matrix <- dist(t(V))
        HC <- hclust(dist.matrix, method="complete")
        V <- V[, HC$order]
        cols <- cols[HC$order]
        dist.matrix <- dist(V)
        HC <- hclust(dist.matrix, method="complete")
        V <- V[HC$order, ]
        rows <- rows[HC$order]
        return(list(V = V, rows = rows, cols = cols))
      }

HeatMapPlot <- function(
V, 
row.names = NULL,
col.names = NULL,
main = " ", 
sub = " ", 
xlab=" ", 
ylab=" ",
row.norm = TRUE,
cmap.type = 1)   # 1 = red/bluecologram, 2 = scale of violets
{
       n.rows <- length(V[,1])
       n.cols <- length(V[1,])

       if (row.norm == TRUE) {
          row.mean <- apply(V, MARGIN=1, FUN=mean)
          row.sd <- apply(V, MARGIN=1, FUN=sd)
          row.n <- length(V[,1])
          for (i in 1:n.rows) {
	     if (row.sd[i] == 0) {
    	         V[i,] <- 0
             } else {
	         V[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
             }
             V[i,] <- ifelse(V[i,] < -6, -6, V[i,])
             V[i,] <- ifelse(V[i,] > 6, 6, V[i,])
          }
        }

       if (cmap.type == 1) { 
           red.blue.palette <- colorRampPalette(c("red", "white", "blue"), space = "rgb")
           mycol <- rev(red.blue.palette(20))
        } else if (cmap.type == 2) {  # range of violet
           violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
           mycol <- rev(violet.palette(20))
        }
        ncolors <- length(mycol) - 2

        heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
        heatm[1:n.rows,] <- V[seq(n.rows, 1, -1),]
        maxv <- max(V)
        minv <- min(V)
        rangev <- maxv - minv
#        windows(width=14, height=9)
        par(mar = c(6, 12, 3, 3))
        image(1:n.cols, 1:n.rows, t(heatm), col=mycol, axes=FALSE, main=main, sub = sub, xlab= xlab, ylab=ylab)
        if (!is.null(row.names)) {
            numC <- nchar(row.names)
            
            size.row.char <- 35/(ifelse(n.rows > 50, 50, n.rows) + 12)
            for (i in 1:n.rows) {
               row.names[i] <- substr(row.names[i], 1, 35)
            }
            axis(2, at=1:n.rows, labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=2, line=-1)
        }
        if (!is.null(col.names)) {
            numC <- nchar(col.names)
            size.col.char <- 20/(n.cols + 15)
            for (i in 1:n.cols) {
               col.names[i] <- substr(col.names[i], 1, 35)
            }
           axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
        }

# version with sorted row and cols

#        dist.matrix <- dist(t(heatm))
#        HC <- hclust(dist.matrix, method="complete")
#        heatm <- heatm[, HC$order]
#        col.names <- col.names[HC$order]
#        dist.matrix <- dist(heatm)
#        HC <- hclust(dist.matrix, method="complete")
#        heatm <- heatm[HC$order, ]
#        row.names <- row.names[HC$order]
#        windows(width=14, height=9)
#        par(mar = c(6, 12, 3, 3))
#        image(1:n.cols, 1:n.rows, t(heatm), col=mycol, axes=FALSE, main=main, sub = sub, xlab= xlab, ylab=ylab)
#        if (!is.null(row.names)) {
#            axis(2, at=1:n.rows, labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=2, line=-1)
#        }
#        if (!is.null(col.names)) {
#            axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
#        }
       
	return()
}

NMF <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {
        N <- length(V[,1])
        M <- length(V[1,])
        set.seed(seed)
        W <- matrix(runif(N*k), nrow = N, ncol = k)  # Initialize W and H with random numbers
        H <- matrix(runif(k*M), nrow = k, ncol = M)
        VP <- matrix(nrow = N, ncol = M)
        error.v <- vector(mode = "numeric", length = maxniter)
        new.membership <- vector(mode = "numeric", length = M)
        old.membership <- vector(mode = "numeric", length = M)
        eps <- .Machine$double.eps
        for (t in 1:maxniter) {
              VP = W %*% H
              H <- H * (crossprod(W, V)/crossprod(W, VP)) + eps
              VP = W %*% H
              H.t <- t(H)
              W <- W * (V %*% H.t)/(VP %*% H.t) + eps
              error.v[t] <- sqrt(sum((V - VP)^2))/(N * M)
               if (t %% stopfreq == 0) {
                    for (j in 1:M) {
                        class <- order(H[,j], decreasing=T)
                        new.membership[j] <- class[1]
                     }
                     if (sum(new.membership == old.membership) == M) {
                        no.change.count <- no.change.count + 1
                     } else {
                        no.change.count <- 0
                     }
                     if (no.change.count == stopconv) break
                     old.membership <- new.membership
               }
        }
        return(list(W = W, H = H, t = t, error.v = error.v))
}

nnls.fit <- function(x,y,wsqrt=1,eps=0,rank.tol=1e-07) {
  ## Purpose: Nonnegative Least Squares (similar to the S-Plus function
  ## with the same name) with the help of the R-library quadprog
  ## ------------------------------------------------------------------------
  ## Attention:
  ## - weights are square roots of usual weights
  ## - the constraint is coefficient>=eps
  ## ------------------------------------------------------------------------
  ## Author: Marcel Wolbers, July 99
  ##
  ##========================================================================
  require ("quadprog")
  m <- NCOL(x)
  if (length(eps)==1) eps <- rep(eps,m)
  x <- x * wsqrt
  y <- y * wsqrt
#  sometimes a rescaling of x and y helps (if solve.QP.compact fails otherwise)
  xscale <- apply(abs(x),2,mean)
  yscale <- mean(abs(y))
  x <- t(t(x)/xscale)
  y <- y/yscale
  Rinv <- backsolve(qr.R(qr(x)),diag(m))
  cf <- solve.QP.compact(Dmat=Rinv,dvec=t(x)%*%y,Amat=rbind(rep(1,m)),
                   Aind=rbind(rep(1,m),1:m),bvec=eps*xscale/yscale,
                         factorized=TRUE)$sol
  cf <- cf*yscale/xscale  #scale back
  cf
}
                                              
write.gct.2 <- function(gct.data.frame, descs = "", filename) 
{
    f <- file(filename, "w")
    cat("#1.2", "\n", file = f, append = TRUE, sep = "")
    cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
    cat("Name", "\t", file = f, append = TRUE, sep = "")
    cat("Description", file = f, append = TRUE, sep = "")

    colnames <- colnames(gct.data.frame)
    cat("\t", colnames[1], file = f, append = TRUE, sep = "")

    if (length(colnames) > 1) {
       for (j in 2:length(colnames)) {
           cat("\t", colnames[j], file = f, append = TRUE, sep = "")
       }
     }
    cat("\n", file = f, append = TRUE, sep = "\t")

    oldWarn <- options(warn = -1)
    m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
    m[, 1] <- row.names(gct.data.frame)
    if (length(descs) > 1) {
        m[, 2] <- descs
    } else {
        m[, 2] <- row.names(gct.data.frame)
    }
    index <- 3
    for (i in 1:dim(gct.data.frame)[2]) {
        m[, index] <- gct.data.frame[, i]
        index <- index + 1
    }
    write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
    close(f)
    options(warn = 0)

}

