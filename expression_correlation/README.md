## Introduction

In these folders we calculate the correlation of the expression profiles between
eIF4A3 and different sets of genes. Expression profiles were obtained from TCGA
samples for different cancer types. Expression profiles were obtained from RNA
sequencing data and pre-processed using the TCGAbiolinks R package as described
in:

	Kumar, M., Papaleo, E.  A pan-cancer assessment of alterations of the kinase
	domain of ULK1, an upstream regulator of autophagy. Sci Rep 10, 14874 (2020). 
	https://doi.org/10.1038/s41598-020-71527-4

We have selected two lists of target genes, one that is just the list of TFEB
targets and one that is the list of genes identified in the gene enrichment
analysis of the RNASeq experiment (see the rnaseq_analysis folder of this
repository). We perform this calculation only for cancer types and
corresponding subtypes for which eIF4A3 is deregulated (abs(logFC > 0.5)). When
plotting and clustering, we further restrict it to those cancer types for which
eIF4A3 is upregulated (logFC > 0.5) hence excluding KICH from the analysis.

The same procedure is repeated for both gene lists, separately. Also, we
distinguish the cases in which we work with samples from cancer types and cancer
subtypes, using slightly different procedures. 

For cancer types (deregulated_eIF4A3_per_type directory) it works as follows:

	* we calculate, for each of the two lists, the Spearman correlation coefficient
	of the expression profiles between eIF4A3 and each gene of the
	genes of the two lists, separately, for tumor primary samples only. We also
	calculate the Spearman's rho statistic as per the cor.test function in R. We
	save all our correlation values and p-values as matrices in csv files.

	* afterwards, we load these matrices and do as follow. We:
		* remove the KICH type (see above)
		* remove types with less than 10 samples
		* remove genes with 50% or more missing data points
		* save the resulting heatmap of correlation values
		* we also filter the correlation values heatmap for the respective
		p-values obtained in the correlation analysis by removing all elements
		with p-value > 0.05, and saving that matrix
		* we cluster the full correlation matrix (see the paper for details) both
		for genes and for cancer type and save the clustering results, both as a
		matrix and as a figure

For cancer subtypes (deregulated_eIF4A3_per_subtype directory), we perform a
similar analysis, albeit with one extra final step. In this case, since we have
less cancer samples per subtype, we try and remove from the matrix before
clustering those gene and cancer subtypes that have very few correlation values
that are deemed to be significant (p-value < 0.05) by the statistical test. This
means that we remove genes and subtypes that have less than 20% of significant
correlations over the total from the matrix, and the cluster it or filter it
similarly as what done with the cancer types.

## Running

just enter one of the two directories above and run the two scripts in this
order:

	Rscript corr_genes_TPonly.R 
	Rscript plot_and_cluster.R

## Output

The final figures that were included in the paper were
heatmap_tfeb_full_clustering_euclidean.pdf for the per type analysis and 
heatmap_tfeb_0.05_0.20_clustering_euclidean.pdf for the subtype analysis. These
include the clustering of the full matrices of all correlation values, with or
without 20% filtering (see above) for cancer subtypes and types respectively
