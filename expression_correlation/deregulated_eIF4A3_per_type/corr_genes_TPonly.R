#!/usr/bin/env R

# corr_genes.R - calculate Spearman correlation between gene expression patterns
# Copyright (C) 2020, Matteo Tiberti, Elena Papaleo
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

library('plyr')
library('testit')
library('TCGAbiolinks')

# get correlation values. Arguments:
#   cancer_type: cancer type
#   data_location: location where we can find the data
#   reference_gene: HGNC of the gene with which correlations will be calculated
#   genes of interest: HGNC of the genes for which correlations will be calculated respect to 


get_correlation = function(cancer_type, data_location, reference_gene, genes_of_interest) {

  cancer_type = cancer_type[1]
  
  print(sprintf("calculating correlation for %s", cancer_type))

  # load gene expression matrix for specific subtype
  dataFilt = get(load(sprintf("%s/%s/filtering/%s_dataFilt.rda",
                        data_location,
                        cancer_type,
                        cancer_type)))
  
  # find tumor primary samples in the dataset
  TP_samples = TCGAquery_SampleTypes(colnames(dataFilt), c("TP"))
  
  # keep only TP samples
  dataFilt = dataFilt[, TP_samples]

  # check if our genes of interest are in there
  if (! reference_gene %in% rownames(dataFilt)) {
    print("reference gene not in the dataset!")
    return(FALSE)
  }
  
  # notify user if some of our genes of interest are 
  # not available, but proceed
  for (g in genes_of_interest) {
    if (! g %in% rownames(dataFilt)) {
      print(paste(g, 'not in the dataset!'))
    }
  }
  
  # keep only those genes in the list that are actually available
  filtered_genes_of_interest = intersect(rownames(dataFilt), genes_of_interest)
  
  # check if any gene of interest is left
  if (length(filtered_genes_of_interest) < 1) {
    printf("no genes of interest left for tumour type %s", project)
    return(FALSE)
  }
          
  # create matrix with just reference gene
  mat_reference_gene = subset(dataFilt, rownames(dataFilt) == reference_gene)
  
  # create matrix with just all the genes of interest
  mat_genes_of_interest = subset(dataFilt, rownames(dataFilt) %in% filtered_genes_of_interest)
    
  # prepare output matrices for this type
  # these matrices will have one row (cancer type) and one column per gene
  # of interest
  this_corr_mat = matrix(data=rep(0, dim(mat_genes_of_interest)[[1]]), nrow=1)
  this_pv_mat = matrix(data=rep(0, dim(mat_genes_of_interest)[[1]]), nrow=1)
  colnames(this_corr_mat) = rownames(mat_genes_of_interest)
  colnames(this_pv_mat) = rownames(mat_genes_of_interest)
  
  # for each gene of interest...
  for (i in 1:dim(mat_genes_of_interest)[1]) {
    # calculate correlation and p-value and store them
    t = cor.test(mat_reference_gene[1,], mat_genes_of_interest[i,], method = 'spearman')
    this_corr_mat[[i]] = t$estimate
    this_pv_mat[[i]] = t$p.value
  }
  
  # add type names, subtype names and number of samples
  # turn into a single-row dataframes and return them
  corr_mats_df = cbind(data.frame(cancer.type=c(cancer_type),
                                  nsamples=c(dim(dataFilt)[2]),
                                  as.data.frame(this_corr_mat)))
  pv_mats_df   = cbind(data.frame(cancer.type=c(cancer_type),
                                  nsamples=c(dim(dataFilt)[2]),
                                  as.data.frame(this_pv_mat)))

  return(list(corr_mats_df, pv_mats_df))
}

get_correlation_matrix = function(cancer_types, data_location, reference_gene, genes_of_interest) {

  # calculate correlation values for each (type, subtype) combination
  # in the cancer_types dataframe
  full_dat = apply(cancer_types, 1, get_correlation, 
    data_location=data_location,
    reference_gene=reference_gene, 
    genes_of_interest=genes_of_interest)

  # get all correlation/pvalues together 
  all_corrs = lapply(full_dat, `[[`, 1)
  all_pvals = lapply(full_dat, `[[`, 2)

  # build final dataframes
  correlations = rbind.fill(all_corrs)
  pvalues = rbind.fill(all_pvals)
  
  return(list(correlations, pvalues))
}

##############################################################################
##############################################################################

# main script

# we exclude CHOL because it has no subtypes defined
excluded_types = c("CHOL")

data_location="../../../../TCGA_data/Expression/cancer_types/"

# load cancer types and filter them for log of fold-change in cancer subtypes
# to find in which subtypes eIF4A3 is deregulated
logFC_co = 0.5
types_DEA = read.table('../../deregulation_eIF4A3/eIF4A3_DEA_per_type.tsv', 
                            sep='\t', 
                            header=TRUE, 
                            stringsAsFactors=FALSE)

# filter out types with low |logFC| and excluded types
cancer_types = data.frame(cancer.type=types_DEA[ (abs(types_DEA$logFC) > logFC_co) & (! types_DEA$cancer.type %in% excluded_types), c("cancer.type")])

# gene against all correlations will be calculated against
reference_gene = 'EIF4A3'

# load list of genes to calculate correlation
genes_of_interest = sort(unique(read.table('../../tfeb_genes.tsv', header=TRUE, stringsAsFactors=FALSE)$genes))

# calculate and save correlation for TFEB genes
data = get_correlation_matrix(cancer_types,
                               data_location,
                               reference_gene,
                               genes_of_interest)

write.table(data[[1]], file="correlation_tfeb.tsv", sep='\t', quote=FALSE)
write.table(data[[2]], file="p-values_tfeb.tsv",    sep='\t', quote=FALSE)

print("the following genes were never in the dataset:")
print(setdiff(genes_of_interest, colnames(data[[1]])))

# calculate and save correlation for genes from enrichment analysis
genes_of_interest = sort(unique(read.table('../../rnaseq_analysis/genes_union_functionally_related_gos.tsv', header=TRUE, stringsAsFactors=FALSE)$genes))

data = get_correlation_matrix(cancer_types,
                               data_location,
                               reference_gene,
                               genes_of_interest)

write.table(data[[1]], file="correlation_enriched.tsv", sep='\t', quote=FALSE)
write.table(data[[2]], file="p-values_enriched.tsv",    sep='\t', quote=FALSE)

print("the following genes were never in the dataset:")
print(setdiff(genes_of_interest, colnames(data[[1]])))

