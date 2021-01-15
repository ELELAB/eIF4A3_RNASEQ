## Introduction

This directory contains summarizing and plotting of the log of fold change
values from differential expression analysis (DEA) for eIF4A3 and different
cancer types and subtypes, in the respective scripts, estimating the difference
in expression between primary tumor and normal samples.

The DEA data used in this directory was originally generated in the context
of this publication:

	Kumar, M., Papaleo, E.  A pan-cancer assessment of alterations of the kinase
	domain of ULK1, an upstream regulator of autophagy. Sci Rep 10, 14874 (2020). 
	https://doi.org/10.1038/s41598-020-71527-4

Briefly, we 

* retained only those cancer types for which at least five normal samples when
we compared tumor primary (TP) tissue samples with normal tissue  (NT) samples.
For the breast cancer tumor type (BRCA) we removed samples belonging to male
patients. We retained datasets with at least five cancer subtype samples when we
carried out analyses at the subtype level. Subtype groups have been defined in
accordance with the guidelines provided by the PanCancerAtlas consortium.

* performed a differential gene expression analysis (DEA) comparing tumor versus
normal samples using the limma-voom pipeline, based on an empirical Bayes
procedure, as we implemented in TCGAbiolinks . We corrected for the Tissue
Source Site (TSS) as a batch effect incorporating it directly in the design
matrix for DEA. we have used an FDR  (adjusted p-value) cutoff of 0.05 and a
logFC cut-off higher than 0.5  (up-regulated genes) or lower than -0.5
(down-regulated genes).

* selected the eIF4A3 gene from the dataset and written the
resulting fold-change values per gene type and subtype (see run_per_type.R and
run_per_subtype.R scripts). We have finally visualized the results using the
plot.R script.

## Running

Run the run_per_type.R and run_per_subtype.R scripts first:

	Rscript run_per_type.R
	Rscript run_per_subtype.R

then run the plotting script:

	Rscript plot.R

## Output

The output figure deregulation_eIF4A3_per_type.pdf was included in the paper
