#!/usr/bin/env R

# corr_genes.R - calculate Spearman correlation between gene expression patterns
# Copyright (C) 2019, Matteo Tiberti <matteo.tiberti@gmail.com>
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

#set working directory
#setwd("/data/user/shared_projects/tcga_eIF4A3/correlations_study/")
#to enquire the required library

library('plyr')
library('R.utils')
library('TCGAutils')
library('TCGAbiolinks')
library('testit')

# round all numeric variables
# x: data frame 
# digits: number of digits to round

get_subtypes_barcodes = function(project, sample_type) {
  
  subtypes_table = PanCancerAtlas_subtypes()
  
  # get subtypes for this project
  subtypes_table = subtypes_table[subtypes_table$cancer.type == project,]
  
  # get available subtypes
  subtypes = unique(subtypes_table$Subtype_Selected)
  print(paste0("identified subtypes for ", project, ": ", subtypes))
  
  # prepare output list
  subtype_barcodes = vector(mode="list", length=length(subtypes))
  names(subtype_barcodes) = subtypes
  
  # find length of stored barcodes - we are enforcing same length for each project
  # this will also raise an error if that cancer type isn't in the subtype table
  barcodes_length = unique(as.vector(sapply(subtypes_table$pan.samplesID, nchar)))
  assert("degenerate barcodes", {
    length(barcodes_length) == 1
  })
  
  # if barcodes are complete just store them
  if (barcodes_length[[1]] == 28) {
    print(paste0("project ", project, " all barcodes are complete. They will not be looked up."))
    for (s in 1:length(subtypes)) {
      subtype_barcodes[[s]] = subtypes_table[subtypes_table$Subtype_Selected == subtypes[s],]$pan.samplesID
    }
  } else {
    # otherwise we need to look up the 
    print(paste0("project ", project, " Barcodes are NOT complete, they will be looked up."))
    tcga_project <- paste0("TCGA-", project)
    for (s in 1:length(subtypes)) {
      barcodes = subtypes_table[subtypes_table$Subtype_Selected == subtypes[s],]$pan.samplesID
      # try to fetch barcodes - here the tryCatch clause is mandatory to handle those cases
      # in wich GDCquery doesn't find any matching barcodes and raises an error. Here
      # we just ignore the error and leave NULL as the barcodes for the corresponding
      # subtype; this will be cleaned up later
      tryCatch({
        stage.query.exp <- GDCquery(project = tcga_project,
                                    data.category = "Transcriptome Profiling",
                                    data.type = "Gene Expression Quantification",
                                    workflow.type = "HTSeq - Counts",
                                    barcode = barcodes,
                                    sample.type = sample_type,
                                    legacy = FALSE)
        subtype_barcodes[[s]] = getResults(stage.query.exp, cols="cases")
        print(paste0("subtype ", subtypes[s], " ", length(subtype_barcodes[[s]]), " found, was ", length(barcodes)))
      },
      error=function(cond) {
        print(paste0("subtype ", subtypes[s], " no samples found"))        
      })
    }
  }
  return(subtype_barcodes)
}


# get correlation values. Arguments:
#   project: TCGA project acronym (BRCA, LUSC, UVM ...)
#   sample_types: same as sample.type argument in GDCquery
#   reference_gene: HGCN of the gene with which correlations will be calculated
#   genes of interes: genes for which correlations will be calculated respect to 

get_gene_matrix = function(dataFilt, project, sample_type, with_normal=FALSE) {
  
  # get short letter code for sample type, after lowercasing to go around
  # TCGA biolinks lowercase descriptions
  this_sampleTypes = as.data.frame(sampleTypes)
  sampleTypes$Definition = tolower(sampleTypes$Definition)
  sample_type_slc = sampleTypes[sampleTypes$Definition == tolower(sample_type),]$Short.Letter.Code
  
  # get barcodes per subtype
  subtype_barcodes = get_subtypes_barcodes(project, sample_type)
  
  # remove subtypes with no barcodes (NULL)
  subtype_barcodes = compact(subtype_barcodes)
  
  # prepare output
  subtype_mats = vector(mode="list", length=length(subtype_barcodes))
  names(subtype_mats) = names(subtype_barcodes)
  
  # get NT samples for this tumor type if necessary
  if (with_normal) {
      NT = TCGAquery_SampleTypes(colnames(dataFilt), "NT")
      print(paste0("NT is ", length(NT) , " samples"))
  }
  
  # for evey subtype, 
  for (st in 1:length(subtype_barcodes)) {
    # filter the barcodes keeping tumor samples only (TP in most cases)
    tumor_samples = TCGAquery_SampleTypes(subtype_barcodes[[st]], sample_type_slc)
    print(subtype_barcodes[[st]])
    print(length(tumor_samples))
    print(tumor_samples)
    # get the corresponding samples that are present in the data
      this_subtype_mat = dataFilt[, colnames(dataFilt) %in% tumor_samples ]
    if (with_normal) {
      # add NT samples if required
      subtype_mats[[st]] = cbind(dataFilt[, colnames(dataFilt) %in% NT ], this_subtype_mat)
    } else {
      subtype_mats[[st]] = this_subtype_mat
    }
  }
  return(subtype_mats)
}

get_correlation = function(project, reference_gene, genes_of_interest, sample_type, with_normal=FALSE,
  data_location="../../../../TCGA_data/Expression") {
  
  print(paste0("calculating correlation for ", project))

  # load gene expression matrix
  dataFilt = get(load(sprintf("%s/%s/filtering/%s_dataFilt.rda", data_location, project, project)))
    
  # check if our genes of interest are in there
  if (! reference_gene %in% rownames(dataFilt)) {
    print("reference gene not in the dataset!")
    return(FALSE)
  }
  
  for (g in genes_of_interest) {
    if (! g %in% rownames(dataFilt)) {
      print(paste(g, 'not in the dataset!'))
    }
  }
  
  # keep only those genes in the list that are actually available
  filtered_genes_of_interest = intersect(rownames(dataFilt), genes_of_interest)
  
  if (length(filtered_genes_of_interest) < 1) {
    printf("no genes of interest left for tumour type %s", project)
    return(FALSE)
  }
  
  subtype_mats = get_gene_matrix(dataFilt, project, sample_type, with_normal=with_normal)
  
  # prepare output
  corr_mats = vector(mode="list", length=length(subtype_mats))
  pv_mats = vector(mode="list", length=length(subtype_mats))
  samples_register = vector(length=length(subtype_mats), mode='numeric')
  
  # for each subtype...
  for (st in 1:length(subtype_mats)) {
    
    # create matrix with just reference gene
    mat_reference_gene = subset(subtype_mats[[st]], rownames(subtype_mats[[st]]) == reference_gene)
    
    # create matrix with just all the genes of interest
    mat_genes_of_interest = subset(subtype_mats[[st]], rownames(subtype_mats[[st]]) %in% filtered_genes_of_interest)
    
    # take note of the number of samples
    samples_register[st] = dim(subtype_mats[[st]])[2]
    
    # prepare output matrices for this subtype
    this_corr_mat = matrix(data=rep(0, dim(mat_genes_of_interest)[[1]]), nrow=1)
    this_pv_mat = matrix(data=rep(0, dim(mat_genes_of_interest)[[1]]), nrow=1)
    colnames(this_corr_mat) = rownames(mat_genes_of_interest)
    colnames(this_pv_mat) = rownames(mat_genes_of_interest)
    
    print(dim(mat_genes_of_interest))
    # for each gene of interest...
    for (i in 1:dim(mat_genes_of_interest)[1]) {
      # calculate correlation and p-value, save them separately
      if (dim(mat_reference_gene)[2] > 1) {
        t = cor.test(mat_reference_gene[1,], mat_genes_of_interest[i,], method = 'spearman')
        this_corr_mat[[i]] = t$estimate
        this_pv_mat[[i]] = t$p.value
      } else {
        this_corr_mat[[i]] = NA
        this_pv_mat[[i]] = NA
      }
    }
    corr_mats[[st]] = this_corr_mat
    pv_mats[[st]] = this_pv_mat
  }
  
  
  # join subtype matrices together in a single dataframe
  corr_mats_df = as.data.frame(rbind.fill.matrix(corr_mats))
  pv_mats_df   = as.data.frame(rbind.fill.matrix(pv_mats))
  
  # add type names, subtype names and number of samples
  corr_mats_df = cbind(data.frame(type=rep(project, dim(corr_mats_df)[1]),
                                  subtypes=names(subtype_mats),
                                  nsamples=samples_register), corr_mats_df)
  pv_mats_df   = cbind(data.frame(type=rep(project, dim(pv_mats_df)[1]), 
                                  subtypes=names(subtype_mats),
                                  nsamples=samples_register), pv_mats_df)
  
  return(list(corr_mats_df, pv_mats_df))
}

get_correlation_matrix = function(cancer_types, reference_gene, genes_of_interest, sample_type, with_normal=FALSE) {
   
  # prepare output correlations dataframe
  correlations = data.frame(matrix(ncol = 3 + length(genes_of_interest), nrow=0), stringsAsFactors = FALSE)
  colnames(correlations) = c("type", "subtypes", "nsamples", genes_of_interest)
  correlations$type = as.character(correlations$type)
  correlations$subtypes = as.character(correlations$subtypes)

  # prepare output p-values dataframe
  pvalues = data.frame(matrix(ncol = 3 + length(genes_of_interest), nrow=0), stringsAsFactors = FALSE)
  colnames(pvalues) = c("type", "subtypes", "nsamples", genes_of_interest)
  pvalues$type = as.character(pvalues$type)
  pvalues$subtypes = as.character(pvalues$subtypes)


  # for each tumour types...
  for (tt in cancer_types) {
    
    # get correlation values
    corr = get_correlation(tt, reference_gene, genes_of_interest, sample_type, with_normal=with_normal)
    if (identical(corr, FALSE)) {
      # if some problem occurred, skip to the next iteration
      # this would better be done with tryCatch()
      next  
    }
    # build up output dataframe
    correlations = rbind.fill(correlations, corr[[1]])
    pvalues = rbind.fill(pvalues, corr[[2]])
  }

  return(list(correlations, pvalues))

}

##############################################################################
##############################################################################

# main script

# load cancer types
cancer_types = read.table('../../cancer_types.tsv', header=TRUE, stringsAsFactors=FALSE)
cancer_types = unique(cancer_types[! is.na(cancer_types$cancer.subtype), ]$cancer.type)

# sample type
sample_type = c('Primary solid Tumor')

# gene against all correlations will be calculated against
reference_gene = 'EIF4A3'

# load list of genes to calculate correlation
genes_of_interest = sort(unique(read.table('../../tfeb_genes.tsv', header=TRUE, stringsAsFactors=FALSE)$genes))

data = get_correlation_matrix(cancer_types, 
                        reference_gene, 
                        genes_of_interest, 
                        sample_type,
                        with_normal=FALSE)

write.table(data[[1]], file="correlation_tfeb.tsv", sep='\t', quote=FALSE)
write.table(data[[2]], file="p-values_tfeb.tsv",   sep='\t', quote=FALSE)


genes_of_interest = sort(unique(read.table('../../rnaseq_analysis/genes_union_all_gos.tsv', header=TRUE, stringsAsFactors=FALSE)$genes))

data = get_correlation_matrix(cancer_types, 
                        reference_gene, 
                        genes_of_interest, 
                        sample_type,
                        with_normal=FALSE)

write.table(data[[1]], file="correlation_enriched.tsv", sep='\t', quote=FALSE)
write.table(data[[2]], file="p-values_enriched.tsv",   sep='\t', quote=FALSE)





