#!/usr/bin/env R

# plot_and_cluster.R - plot and cluster correlation data
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

library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(factoextra)

analyze_correlation = function(name, alpha) {

    # color scale 
    colf=colorRamp2(c(-1,       0,      1   ),
                    c("blue", "white", "red"))

    corrs_fname = sprintf("correlation_%s.tsv", name)
    pvals_fname = sprintf("p-values_%s.tsv", name)

    # read data
    corrs = read.table(corrs_fname, stringsAsFactors=FALSE, sep='\t', quote="")
    pvals = read.table(pvals_fname, stringsAsFactors=FALSE, sep='\t', quote="")

    rownames(corrs) = corrs$cancer.type
    rownames(pvals) = pvals$cancer.type

    # remove genes that are NA in at least
    # half of the values (for clustering) as well as  
    # types that have less than 10 samples

    orig_rownames = rownames(corrs)
    orig_colnames = colnames(corrs)

    print(sprintf("original matrix %s: %d, %d", name, dim(corrs)[[1]], dim(corrs)[[2]]))
    corrs = corrs[, colSums(! is.na(corrs)) > 4 ]
    corrs = corrs[corrs$nsamples >= 10, ]
    pvals = pvals[, colSums(! is.na(pvals)) > 4 ]
    pvals = pvals[pvals$nsamples >= 10, ]
    print("lost rows")
    print(setdiff(orig_rownames, rownames(corrs)))
    print("lost cols")
    print(setdiff(orig_colnames, colnames(corrs)))

    #############################

    # generate matrices
    corrs_m = as.matrix(corrs[,3:dim(corrs)[2]])
    pvals_m = as.matrix(pvals[,3:dim(pvals)[2]])

    # plot
    ha = columnAnnotation(samples = anno_text(corrs$nsamples))
    pdf_fname = sprintf("heatmap_%s_full.pdf", name)
    pdf(pdf_fname, width=10, height=12)
    hm = Heatmap(t(corrs_m),
            col=colf,
            na_col="grey",
            cluster_columns=F,
            cluster_rows=F,
            column_names_side = "top",
            top_annotation=ha,
            row_names_gp = gpar(fontsize = 8))

    draw(hm)
    dev.off()

    # euclidean cluster and plot
    pdf_fname = sprintf("heatmap_%s_full_clustering_euclidean.pdf", name)
    pdf(pdf_fname, width=10, height=12)
    hm = Heatmap(t(corrs_m), 
            col=colf,
            cluster_columns=T,
            cluster_rows=T,
            column_names_side = "top",
            row_names_gp = gpar(fontsize = 8),
            column_dend_height = unit(4, "cm"),
            row_dend_width = unit(4, "cm"),
            column_dend_side="top",
            row_dend_side = "right")
    draw(hm)
    dev.off()
    
    tsv_fname = sprintf("heatmap_%s_full_clustering_euclidean.tsv", name)
    write.table(t(corrs_m)[row_order(hm), column_order(hm)], file=tsv_fname, sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
    
    corrs_m[pvals_m > alpha] = NA
    
    corrs_m = corrs_m[column_order(hm), row_order(hm)] # because I'm plotting the transposed matrix
    
    pdf_fname = sprintf("heatmap_%s_clustering_euclidean_filt.pdf", name)
    pdf(pdf_fname, width=10, height=12)
    hm = Heatmap(t(corrs_m),
                 na_col="grey", 
                 col=colf,
                 cluster_columns=F,
                 cluster_rows=F,
                 column_names_side = "top",
                 #column_labels=subtypes,
                 #column_split=corrs$type,
                 #top_annotation=ha,
                 row_names_gp = gpar(fontsize = 8))
    draw(hm)
    dev.off()
    
    
    #############################

    # regenerate matrices
    corrs_m = as.matrix(corrs[,3:dim(corrs)[2]])
    pvals_m = as.matrix(pvals[,3:dim(pvals)[2]])

    # filter out elements by p-value
    corrs_m[pvals_m > alpha] = NA

    # create vectors with fraction of significant values 
    cols_full = colSums( (! is.na(corrs_m))/dim(corrs_m)[1] )
    rows_full = rowSums( (! is.na(corrs_m))/dim(corrs_m)[2] )
    
    # calculate empirical cumulative distribution functions 
    pdf_fname = sprintf("rows_ecdf_%s.pdf", name)
    title = sprintf("ECDF for rows, %s gene set", name)
    pdf(pdf_fname)
    plot(ecdf(as.vector(rows_full)), main=title)
    dev.off()

    pdf_fname = sprintf("cols_ecdf_%s.pdf", name)
    title = sprintf("ECDF for columns, %s gene set", name)
    pdf(pdf_fname)
    plot(ecdf(as.vector(cols_full)), main=title)
    dev.off()

    # plot filtered matrix
    pdf_fname = sprintf("heatmap_%s_full_filt.pdf", name)
    pdf(pdf_fname, width=10, height=12)
    hm = Heatmap(t(corrs_m),
            col=colf,
            na_col="grey",
            cluster_columns=F,
            cluster_rows=F,
            column_names_side = "top",
            top_annotation=ha,
            row_names_gp = gpar(fontsize = 8))
    draw(hm)
    dev.off()

    tsv_name = sprintf("heatmap_%s_full_filt.tsv", name)
    write.table(t(corrs_m)[row_order(hm), column_order(hm)], file=tsv_fname, sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
    print(sprintf("full matrix %s, %d genes, %d types", name, dim(corrs_m)[[2]], dim(corrs_m)[[1]]))
    print(sprintf("full matrix %s filt, %d cells remain (%.1f %%)", name, sum( ! is.na(corrs_m)), 
                                                               sum( ! is.na(corrs_m)) / (dim(corrs_m)[[1]]*dim(corrs_m)[[2]])*100))
}
    #############################

# run analysis
analyze_correlation("tfeb", 0.05)

analyze_correlation("enriched", 0.05)
