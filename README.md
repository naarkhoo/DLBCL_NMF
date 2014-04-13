DLBCL_NMF
=========

Introduction
This pipeline calculates clusters based on a consensus non-negative matrix factorization (NMF) clustering method , . This pipeline has the following features:

1) Convert input data set to a non-negitive matrix by column rank normalization.
2) Classify samples into consensus clusters.
3) Determine differentially expressed marker genes for each subtype.

Input:
mRNAseq of normalized RSEM/RPKM value with log2 transformed was as the input RNAseq data for the clustering.
RSEM is used to estimate gene and transcript abundances and these values are normalized to a fixed upper quaritile value of 1000 for gene and 300 for transcript level estimates.
RPKM for a given GeneX is calculated by: (raw read counts * 10^9) / (total reads * length of GeneX). Total reads is the lane yield after removing poor quality reads and the length of GeneX is defined as the median length of all transcripts associated with GeneX.



CNMF method:
Non-negative matrix factorization (NMF) is an unsupervised learning algorithm that has been shown to identify molecular patterns when applied to gene expression data , . Rather than separating gene clusters based on distance computation, NMF detects contextdependent patterns of gene expression in complex biological systems.

Coephentic Correlation Coefficient:
We use the cophenetic correlation coefficient to determine the cluster that yields the most robust clustering. The cophenetic correlation coefficient is computed based on the consensus matrix of the CNMF clustering, and measures how reliably the same samples are assigned to the same cluster across many iterations of the clustering lgorithm with random initializations. The cophenetic correlation coefficient lies between 0 and 1, with higher values indicating more stable cluster assignments. We select the number of clusters k based on the largest observed correlation coefficient for all tested values of k.


References:
[1] Brunet, J.P., Tamayo, P., Golub, T.R. & Mesirov, J.P., Metagenes and molecular pattern discovery using matrix factorization, Proc Natl Acad Sci U S A 12(101):4164-9 (2004)
[2] Broad Genepattern: NMFConsensus
[3] Rousseeuw, P.J., Silhouettes: A graphical aid to the interpretation and validation of cluster analysis., J. Comput. Appl. Math. 20:53-65 (1987)
