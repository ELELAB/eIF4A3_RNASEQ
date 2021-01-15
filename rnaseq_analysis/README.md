## Introduction

Here we start from an already-performed differential expression analysis (DEA)
as stored in the control-VS-case1.DEseq2_Method.GeneDiffExp_edit.csv and
control-VS-case2.DEseq2_Method.GeneDiffExp_edit.csv files, that correspond
respectively to DEA analysis performed on RNA sequencing data of two cell lines
in which eIF4A3 was depleted by means of two different siRNA (see paper). Data
pre-processing and DEA was performed by BGI (see paper for details).

We have first saved the original Excel files in the csv format to make them
easier to work on.

In the run.R script we, separately for up and down-regulated genes:

* Keep those genes that are considered to be up or downregulated respect to
control in either case, from which we get a list of genes

* use enrichR to perform gene enrichment analysis on such list, on the 
Biological processes and Cellular components gene ontology (libraries
GO_Biological_Process_2018 and GO_Cellular_Components_2018)

* sort the resulting enriched GO terms by P-value and consider the top 10 for
plotting

* generate a Venn diagram detailing the overlap of upregulated genes only for
both cases - this was not done for downregulated genes

We plot both the first 10 categories as described above, and a more
comprehensive heatmap with genes that
* were found to be significantly upregulated
* feature GO terms that were found to be significantly enriched in the
dataset 
* feature Go terms that were related to the biological phenomenon under
study. These are:
	lysosome (GO:0005764)
	lytic vacuole (GO:0000323)
	lytic vacuole membrane (GO:0098852)
	lysosomal membrane (GO:0005765)
	autophagosome (GO:0005776)
	autophagosome membrane (GO:0000421)

we save all data as csv files as well.

## Running

just run the run.R script that performs all the required steps:

	Rscript run.r

## Output

the final plots that were used in the paper were categories_up_first_10_GO_cell_comp.pdf
and heatmap_union_functionally_related_gos.pdf
