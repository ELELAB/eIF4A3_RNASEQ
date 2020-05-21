
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(factoextra)

excluded_subtypes = c("BRCA.Normal",
                   "GBM_LGG.NA",
                   "KIRC.NA",
                   "THCA.NA",
                   "UCEC.NA")

corrs = read.table("correlation_enriched.tsv", stringsAsFactors=FALSE, sep='\t', quote="")
pvals = read.table("p-values_enriched.tsv", stringsAsFactors=FALSE, sep='\t', quote="")
rownames(corrs) = paste(corrs$type, corrs$subtypes, sep="_")
rownames(pvals) = paste(pvals$type, pvals$subtypes, sep="_")

# generate mixed type-subtype names

# remove genes that are NA everywhere
#corrs = corrs[, colSums(! is.na(corrs)) > 0 ]
corrs = corrs[! corrs$subtypes %in% excluded_subtypes, ]
corrs = corrs[corrs$nsamples >= 10, ]
#pvals = pvals[, colSums(! is.na(pvals)) > 0 ]
pvals = pvals[! pvals$subtypes %in% excluded_subtypes, ]
pvals = pvals[pvals$nsamples >= 10,]

#############################

al = 0.05

corrs_m = as.matrix(corrs[,4:dim(corrs)[2]])
pvals_m = as.matrix(pvals[,4:dim(pvals)[2]])

corrs_m[pvals_m > al] = NA

cols_full = colSums( (! is.na(corrs_m))/dim(corrs_m)[1] )
rows_full = rowSums( (! is.na(corrs_m))/dim(corrs_m)[2] )

png('rows_ecdf_enriched.png')
plot(ecdf(rows_full))
dev.off()

png('cols_ecdf_enriched.png')
plot(ecdf(cols_full))
dev.off()

#############################

corrs_m = as.matrix(corrs[,4:dim(corrs)[2]])
pvals_m = as.matrix(pvals[,4:dim(pvals)[2]])

corrs_m[pvals_m > al] = NA

rows_not_NA = rowSums(! is.na(corrs_m)) > 0
cols_not_NA = colSums(! is.na(corrs_m)) > 0

cols_full = colSums( (! is.na(corrs_m))/dim(corrs_m)[1] ) > 0.2
rows_full = rowSums( (! is.na(corrs_m))/dim(corrs_m)[2] ) > 0.2

rows_idxs = rows_not_NA & rows_full
cols_idxs = cols_not_NA & cols_full

corrs_m2 = as.matrix(corrs[,4:dim(corrs)[2]])
corrs_m2 = corrs_m2[rows_idxs, cols_idxs]

pcs = prcomp(corrs_m2, center=TRUE, scale=TRUE)

png('eigval_subtypes_0.2_enriched.png')
fviz_eig(pcs)
dev.off()

png('proj12_subtypes_0.2_enriched.png')
fviz_pca_ind(pcs, repel=TRUE, axes=c(1,2), labelsize=3)
dev.off()

png('proj13_subtypes_0.2_enriched.png')
fviz_pca_ind(pcs, repel=TRUE, axes=c(1,3), labelsize=3)
dev.off()

colf=colorRamp2(c(-1,       0,      1   ), 
                c("blue", "white", "red"))

#ha = columnAnnotation(samples = anno_text(nsamples))
png("heatmap_correlations_0.05_0.2_TPonly_cluster_enriched_euclidean.png", width=1000, height=800)
hm = Heatmap(t(corrs_m2), 
        col=colf,
        cluster_columns=T,
        cluster_rows=T,
        #clustering_distance_columns = "spearman",
        #clustering_distance_rows = "spearman",
        column_names_side = "top",
        #column_labels=subtypes,
        #column_split=corrs$type,
        #top_annotation=ha,
        row_names_gp = gpar(fontsize = 8),
        column_dend_height = unit(4, "cm"),
        row_dend_width = unit(4, "cm"),
        column_dend_side="top",
        row_dend_side = "right")
draw(hm)
dev.off()

genes_0.05_0.2_TPonly_cluster_euclidean = rownames(t(corrs_m2))[row_order(hm)]

corrs_m = corrs_m[rows_idxs, cols_idxs]
corrs_m = corrs_m[column_order(hm), row_order(hm)] # because I'm plotting the transposed matrix

png("heatmap_correlations_0.05_0.2_TPonly_cluster_enriched_euclidean_filt.png", width=1000, height=800)
Heatmap(t(corrs_m),
        na_col="grey", 
        col=colf,
        cluster_columns=F,
        cluster_rows=F,
        column_names_side = "top",
        #column_labels=subtypes,
        #column_split=corrs$type,
        #top_annotation=ha,
        row_names_gp = gpar(fontsize = 8))
dev.off()

png("heatmap_correlations_0.05_0.2_TPonly_cluster_enriched_pearson.png", width=1000, height=800)
hm = Heatmap(t(corrs_m2), 
        col=colf,
        na_col="grey", 
        cluster_columns=T,
        cluster_rows=T,
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",
        column_names_side = "top",
             #column_labels=subtypes,
             #column_split=corrs$type,
             #top_annotation=ha,
        row_names_gp = gpar(fontsize = 8),
        column_dend_height = unit(4, "cm"),
        row_dend_width = unit(4, "cm"),
        column_dend_side="top",
        row_dend_side = "right")
draw(hm)
dev.off()

genes_0.05_0.2_TPonly_cluster_pearson = rownames(t(corrs_m2))[row_order(hm)]

corrs_m = corrs_m[column_order(hm), row_order(hm)] # because I'm plotting the transposed matrix

png("heatmap_correlations_0.05_0.2_TPonly_cluster_enriched_pearson_filt.png", width=1000, height=800)
Heatmap(t(corrs_m), 
        col=colf,
        na_col="grey", 
        cluster_columns=F,
        cluster_rows=F,
        column_names_side = "top",
        #column_labels=subtypes,
        #column_split=corrs$type,
        #top_annotation=ha,
        row_names_gp = gpar(fontsize = 8))
dev.off()

png("heatmap_correlations_0.05_0.2_TPonly_cluster_enriched_spearman.png", width=1000, height=800)
hm = Heatmap(t(corrs_m2), 
        col=colf,
        na_col="grey", 
        cluster_columns=T,
        cluster_rows=T,
        clustering_distance_columns = "spearman",
        clustering_distance_rows = "spearman",
        column_names_side = "top",
             #column_labels=subtypes,
             #column_split=corrs$type,
             #top_annotation=ha,
        row_names_gp = gpar(fontsize = 8),
        column_dend_height = unit(4, "cm"),
        row_dend_width = unit(4, "cm"),
        column_dend_side="top",
        row_dend_side = "right")
draw(hm)
dev.off()

corrs_m = corrs_m[column_order(hm), row_order(hm)] # because I'm plotting the transposed matrix

png("heatmap_correlations_0.05_0.2_TPonly_cluster_enriched_spearman_filt.png", width=1000, height=800)
Heatmap(t(corrs_m), 
        col=colf,
        na_col="grey", 
        cluster_columns=F,
        cluster_rows=F,
        column_names_side = "top",
        #column_labels=subtypes,
        #column_split=corrs$type,
        #top_annotation=ha,
        row_names_gp = gpar(fontsize = 8))
dev.off()

write.table(data.frame(euclidean=genes_0.05_0.2_TPonly_cluster_euclidean,
                     pearson=genes_0.05_0.2_TPonly_cluster_pearson),
                     'GO_genes_clustering_0.2.csv',
                     row.names=FALSE,
                     quote=FALSE,
                     sep='\t')



#############################
#############################
#############################
#############################




corrs_m = as.matrix(corrs[,4:dim(corrs)[2]])
pvals_m = as.matrix(pvals[,4:dim(pvals)[2]])

corrs_m[pvals_m > al] = NA

rows_not_NA = rowSums(! is.na(corrs_m)) > 0
cols_not_NA = colSums(! is.na(corrs_m)) > 0

cols_full = colSums( (! is.na(corrs_m))/dim(corrs_m)[1] ) > 0.15
rows_full = rowSums( (! is.na(corrs_m))/dim(corrs_m)[2] ) > 0.15

rows_idxs = rows_not_NA & rows_full
cols_idxs = cols_not_NA & cols_full

corrs_m2 = as.matrix(corrs[,4:dim(corrs)[2]])
corrs_m2 = corrs_m2[rows_idxs, cols_idxs]

pcs = prcomp(corrs_m2, center=TRUE, scale=TRUE)

png('eigval_subtypes_0.15_enriched.png')
fviz_eig(pcs)
dev.off()

png('proj12_subtypes_0.15_enriched.png')
fviz_pca_ind(pcs, repel=TRUE, axes=c(1,2), labelsize=3)
dev.off()

png('proj13_subtypes_0.15_enriched.png')
fviz_pca_ind(pcs, repel=TRUE, axes=c(1,3), labelsize=3)
dev.off()

colf=colorRamp2(c(-1,       0,      1   ), 
                c("blue", "white", "red"))

png("heatmap_correlations_0.05_0.15_TPonly_cluster_enriched_euclidean.png", width=1000, height=800)
hm = Heatmap(t(corrs_m2), 
             col=colf,
             cluster_columns=T,
             cluster_rows=T,
             #clustering_distance_columns = "spearman",
             #clustering_distance_rows = "spearman",
             column_names_side = "top",
             #column_labels=subtypes,
             #column_split=corrs$type,
             #top_annotation=ha,
             row_names_gp = gpar(fontsize = 8),
             column_dend_height = unit(4, "cm"),
             row_dend_width = unit(4, "cm"),
             column_dend_side="top",
             row_dend_side = "right")
draw(hm)
dev.off()

genes_0.05_0.15_TPonly_cluster_euclidean = rownames(t(corrs_m2))[row_order(hm)]

corrs_m = corrs_m[rows_idxs, cols_idxs]
corrs_m = corrs_m[column_order(hm), row_order(hm)] # because I'm plotting the transposed matrix

png("heatmap_correlations_0.05_0.15_TPonly_cluster_enriched_euclidean_filt.png", width=1000, height=800)
Heatmap(t(corrs_m),
        na_col="grey", 
        col=colf,
        cluster_columns=F,
        cluster_rows=F,
        column_names_side = "top",
        #column_labels=subtypes,
        #column_split=corrs$type,
        #top_annotation=ha,
        row_names_gp = gpar(fontsize = 8))
dev.off()

png("heatmap_correlations_0.05_0.15_TPonly_cluster_enriched_pearson.png", width=1000, height=800)
hm = Heatmap(t(corrs_m2), 
             col=colf,
             na_col="grey", 
             cluster_columns=T,
             cluster_rows=T,
             clustering_distance_columns = "pearson",
             clustering_distance_rows = "pearson",
             column_names_side = "top",
             #column_labels=subtypes,
             #column_split=corrs$type,
             #top_annotation=ha,
             row_names_gp = gpar(fontsize = 8),
             column_dend_height = unit(4, "cm"),
             row_dend_width = unit(4, "cm"),
             column_dend_side="top",
             row_dend_side = "right")
draw(hm)
dev.off()

genes_0.05_0.15_TPonly_cluster_pearson = rownames(t(corrs_m2))[row_order(hm)]

corrs_m = corrs_m[column_order(hm), row_order(hm)] # because I'm plotting the transposed matrix

png("heatmap_correlations_0.05_0.15_TPonly_cluster_enriched_pearson_filt.png", width=1000, height=800)
Heatmap(t(corrs_m), 
        col=colf,
        na_col="grey", 
        cluster_columns=F,
        cluster_rows=F,
        column_names_side = "top",
        #column_labels=subtypes,
        #column_split=corrs$type,
        #top_annotation=ha,
        row_names_gp = gpar(fontsize = 8))
dev.off()

png("heatmap_correlations_0.05_0.15_TPonly_cluster_enriched_spearman.png", width=1000, height=800)
hm = Heatmap(t(corrs_m2), 
             col=colf,
             na_col="grey", 
             cluster_columns=T,
             cluster_rows=T,
             clustering_distance_columns = "spearman",
             clustering_distance_rows = "spearman",
             column_names_side = "top",
             #column_labels=subtypes,
             #column_split=corrs$type,
             #top_annotation=ha,
             row_names_gp = gpar(fontsize = 8),
             column_dend_height = unit(4, "cm"),
             row_dend_width = unit(4, "cm"),
             column_dend_side="top",
             row_dend_side = "right")
draw(hm)
dev.off()

corrs_m = corrs_m[column_order(hm), row_order(hm)] # because I'm plotting the transposed matrix

png("heatmap_correlations_0.05_0.15_TPonly_cluster_enriched_spearman_filt.png", width=1000, height=800)
Heatmap(t(corrs_m), 
        col=colf,
        na_col="grey", 
        cluster_columns=F,
        cluster_rows=F,
        column_names_side = "top",
        #column_labels=subtypes,
        #column_split=corrs$type,
        #top_annotation=ha,
        row_names_gp = gpar(fontsize = 8))
dev.off()

write.table(data.frame(euclidean=genes_0.05_0.15_TPonly_cluster_euclidean,
                     pearson=genes_0.05_0.15_TPonly_cluster_pearson),
          'GO_genes_clustering_0.15.csv',
          row.names=FALSE,
          quote=FALSE,
          sep='\t')

