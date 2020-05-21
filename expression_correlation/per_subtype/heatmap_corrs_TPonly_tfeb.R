
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

excluded_subtypes = c("BRCA.Normal",
                   "GBM_LGG.NA",
                   "KIRC.NA",
                   "THCA.NA",
                   "UCEC.NA")

corrs = read.table("correlation_tfeb.tsv", stringsAsFactors=FALSE, sep='\t', quote="")
pvals = read.table("p-values_tfeb.tsv", stringsAsFactors=FALSE, sep='\t', quote="")

# remove genes that are NA everywhere
corrs = corrs[, colSums(! is.na(corrs)) > 0 ]
corrs = corrs[! corrs$subtypes %in% excluded_subtypes, ]
pvals = pvals[, colSums(! is.na(pvals)) > 0 ]
pvals = pvals[! pvals$subtypes %in% excluded_subtypes, ]

#################################

corrs_m = as.matrix(corrs[,4:dim(corrs)[2]])
pvals_m = as.matrix(pvals[,4:dim(pvals)[2]])

colf=colorRamp2(c(-1,       0,      1   ), 
                c("blue", "white", "red"))

ha = columnAnnotation(samples = anno_text(corrs$nsamples))
png("heatmap_correlations_full_TPonly_tfeb.png", width=1000, height=800)
Heatmap(t(corrs_m), 
        col=colf,
        na_col="grey", 
        cluster_columns=F,
        cluster_rows=F,
        column_names_side = "top",
        column_labels=corrs$subtypes,
        column_split=corrs$type,
        top_annotation=ha,
        row_names_gp = gpar(fontsize = 8))
dev.off()

values = as.data.frame(as.vector(corrs_m))
colnames(values) = c('v1')

png("distribution_correlations_full_TPonly_tfeb.png")
ggplot(values, aes(v1)) + 
    geom_histogram(aes(y=..density..), na.rm=TRUE) +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()


#############################

corrs_m = as.matrix(corrs[,4:dim(corrs)[2]])
pvals_m = as.matrix(pvals[,4:dim(pvals)[2]])

alpha=0.05

corrs_m[pvals_m > alpha] = NA

min_nsamples = corrs$nsamples >= 10

row_selection = which(rowSums(! is.na(corrs_m)) > 0 & min_nsamples, )
corrs_m = corrs_m[row_selection,]
corrs = corrs[row_selection,]

colf=colorRamp2(c(-1,       0,      1   ), 
                c("blue", "white", "red"))

ha = columnAnnotation(samples = anno_text(corrs$nsamples))
png("heatmap_correlations_0.05_TPonly_tfeb.png", width=1000, height=800)
Heatmap(t(corrs_m), 
        col=colf,
        na_col="grey", 
        cluster_columns=F,
        cluster_rows=F,
        column_names_side = "top",
        column_labels=corrs$subtypes,
        column_split=corrs$type,
        top_annotation=ha,
        row_names_gp = gpar(fontsize = 8))
dev.off()

values = as.data.frame(as.vector(corrs_m))
colnames(values) = c('v1')

png("distribution_correlations_0.05_TPonly_tfeb.png")
ggplot(values, aes(v1)) + 
    geom_histogram(aes(y=..density..),na.rm=TRUE) +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()

########################################

corrs = read.table("correlation_tfeb.tsv", stringsAsFactors=FALSE, sep='\t', quote="")
corrs = corrs[, colSums(! is.na(corrs)) > 0 ]
corrs = corrs[! corrs$subtypes %in% excluded_subtypes, ]

corrs_m = as.matrix(corrs[,4:dim(corrs)[2]])
corrs_m = corrs_m[row_selection,]
corrs =   corrs[row_selection,]

colf=colorRamp2(c(-1,       0,      1   ), 
                c("blue", "white", "red"))

ha = columnAnnotation(samples = anno_text(corrs$nsamples))
png("heatmap_correlations_full_TPonly_not_filtered_tfeb.png", width=1000, height=800)
Heatmap(t(corrs_m), 
        col=colf,
        na_col="grey", 
        cluster_columns=F,
        cluster_rows=F,
        column_names_side = "top",
        column_labels=corrs$subtypes,
        column_split=corrs$type,
        top_annotation=ha,
        row_names_gp = gpar(fontsize = 8))
dev.off()

values = as.data.frame(as.vector(corrs_m))
colnames(values) = c('v1')

png("distribution_correlations_full_TPonly_not_filtered_tfeb.png")
ggplot(values, aes(v1)) + 
    geom_histogram(aes(y=..density..), na.rm=TRUE) +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()


