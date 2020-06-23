# run.R - analyze RNAseq differential expression analysis data
# Copyright (C) 2020 Matteo Tiberti, Elena Papaleo

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(viridis)
library(stringr)
library(enrichR)
library(tidyr)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(scales)
library(VennDiagram)
library(RColorBrewer)

plot_heatmap <- function(enrichment_table, 
                    fc_matrix, 
                    tfeb_genes, 
                    int_gos, 
                    tfeb_color, 
                    output_file,
                    output_gene_file,
                    width,
                    height) {

    # here we are plotting a ComplexHeatmap with two parts (heatmap and an extra annotation).
    # The heatmap shows the fold-change of the identified genes and the annotation
    # shows their association with a certain GO term of interest as identified by the
    # enrichment analysis.

    # we first start working on (ii). here we transform the enrichment dataframe to a dataframe in which 
    # each line is a single gene from the FC matrix, each column is a GO term of interest
    # and the entries contain Y or N if the gene is associated with the GO term of interest
    # according to the enrichment analysis

    all_gos = data.frame()

    for (g in 1:length(int_gos)) {
        go_genes = str_split(enrichment_table[ enrichment_table["Term"] == int_gos[[g]], ]$Genes, ';')[[1]]
        df_col = ifelse(rownames(fc_matrix) %in% go_genes, 'Y', 'N')
        df = data.frame(df_col)
        colnames(df) = c(g)
        if (g ==1 ) {
            all_gos = df
        }
        else {
            all_gos = cbind(all_gos, 
                            df)
        }
    }
    colnames(all_gos) = int_gos
    rownames(all_gos) = rownames(fc_matrix)

    # here we start preparing (i) as well.
    # we filter out those lines (genes) which have no associated GO term among those
    # of interest from the FC matrix
    any_go = apply(all_gos, 1, function(x)(any(x=='Y')))
    fc_matrix = fc_matrix[any_go,]

    # write gene list as tsv
    fc_df = data.frame(genes=rownames(fc_matrix))
    write.table(fc_df, output_gene_file, sep='\t', quote=FALSE, row.names=FALSE)

    # likewise, we select in the all_gos matrix only those genes that are left in the
    # FC dataframe
    all_gos_any = all_gos[rownames(all_gos) %in% rownames(fc_matrix), ]

    # prepare colors for the annotation
    ha_col = c("Y" = 'red', "N" = 'white')
    ha_cols = list(ha_col)[rep(1,length(int_gos))]
    names(ha_cols) = colnames(all_gos)

    # prepare the annotation
    ha = rowAnnotation(df = all_gos_any,
                       col = ha_cols,
                       show_legend=F)

    # color the names of the TFEB genes of a different color, specified by
    # tfeb_color. Note that default is black (= don't color them)
    tfeb_color = ifelse(rownames(fc_matrix) %in% tfeb_genes, tfeb_color, 'black')

    # set up legend
    hm_legend_parms = list(title="log(FC)",
                           legend_height=unit(7, "cm"),
                           grid_width=unit(2, "cm"))

    # plot heatmap with annotation
    pdf(output_file, width=width, height=height)
    hm = Heatmap(fc_matrix,
                 cluster_columns=F,
                 cluster_rows=F,
                 row_names_gp = gpar(col = tfeb_color, fontsize=12),
                 col=viridis(1000),
                 right_annotation=ha,
                 heatmap_legend_param=hm_legend_parms)
    draw(hm)
    dev.off()
}

# load list of TFEB-controlled genes
tfeb_genes = sort(unique(read.table('../tfeb_genes.tsv', header=TRUE, stringsAsFactors=FALSE)$genes))

# load data
cases_cols = c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "character", "numeric", "character")
case1 = read.table('control-VS-case1.DEseq2_Method.GeneDiffExp_edit.csv', sep=';', header=T, colClasses=cases_cols)
case2 = read.table('control-VS-case2.DEseq2_Method.GeneDiffExp_edit.csv', sep=';', header=T, colClasses=cases_cols)

# identify upregulated genes
case1_up   = case1[ case1$classification == 'Up' ,]
case2_up   = case2[ case2$classification == 'Up' ,]

# plot venn diagram
myCol = brewer.pal(4, "Pastel2")

venn.diagram(list(c1=case1_up$Symbol, c2=case2_up$Symbol),
             filename="venn_upregulated.svg", 
             category.names=list("Case 1","Case 2"),
             fontfamily = "sans",
             cat.fontfamily = "sans",
             lty = 'blank',
             imagetype='svg',
             fill = myCol[1:2])

venn.diagram(list(c1=case1_up$Symbol, c2=case2_up$Symbol),
             filename="venn_upregulated.tiff", 
             category.names=list("Case 1","Case 2"),
             fontfamily = "sans",
             cat.fontfamily = "sans",
             lty = 'blank',
             imagetype='tiff',
             fill = myCol[1:2])


# list of genes that are upregulated in either experiments
union_genes = unique(c(case1_up$Symbol, case2_up$Symbol))

# run enrichR on such gene list
enrichr_dbs = c("GO_Biological_Process_2018", "GO_Cellular_Component_2018")
enrichr_dbs_labels = c("GO_biol_proc", "GO_cell_comp")

up_enrichment = enrichr(union_genes, enrichr_dbs)

for (i in range(1, length(enrichr_dbs))) {

    this_enrichment = up_enrichment[[i]]
    write.table(this_enrichment, sprintf("enrichment_up_%s.tsv", enrichr_dbs_labels[i]), sep='\t', row.names=FALSE, quote=FALSE)
    
    this_enrichment = this_enrichment[ this_enrichment$P.value <= 0.05, ]
    
    # sort by P-value
    sorted_enrichment = this_enrichment[order(this_enrichment$P.value, decreasing=TRUE),]
    
    # number of genes per category
    gene_count = lapply(sorted_enrichment$Genes, function(x) length(strsplit(x, ';' )[[1]]))
    sorted_enrichment$gene.count = gene_count
    
    # unlist the unlistable to make it compatible with ggplot
    sorted_enrichment$P.value = unlist(sorted_enrichment$P.value)
    sorted_enrichment$Term = unlist(sorted_enrichment$Term)
    sorted_enrichment$gene.count = unlist(sorted_enrichment$gene.count)
    
    sorted_enrichment$Term = factor(sorted_enrichment$Term, levels=sorted_enrichment$Term)
    
    # do the plotting of GO terms, number of genes and P-pvalues
    #pdf(height=3000, width=1000)
    plot = ggplot(data=sorted_enrichment, aes_string(x="P.value", y="Term", size="gene.count")) + 
    geom_point() + 
    scale_size(range=c(0, 6)) + 
    ylab(NULL) +
    scale_x_continuous( trans = "log10",
                        breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x)))
    ggsave(sprintf("categories_up_%s.pdf", enrichr_dbs_labels[i]), plot=plot, width=10, height=dim(sorted_enrichment)[1]/8)
    
    df_l = dim(sorted_enrichment)[1]
    plot = ggplot(data=sorted_enrichment[(df_l-9):df_l,], aes_string(x="P.value", y="Term", size="gene.count")) + 
    geom_point() + 
    scale_size(range=c(0, 6)) + 
    ylab(NULL) +
    scale_x_continuous( trans = "log10", 
                        breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x)))
    ggsave(sprintf("categories_up_first_10_%s.pdf", enrichr_dbs_labels[i]), plot=plot, width=10, height=3.5)
}

this_enrichment = up_enrichment[[2]]
this_enrichment = this_enrichment[ this_enrichment$P.value <= 0.05, ]

# prepare logFC matrix to be plotted
# get from both dataframes those genes that are in the union list
case1_overlap = case1[ case1$Symbol %in% union_genes, c("Symbol", "logFC")]
case2_overlap = case2[ case2$Symbol %in% union_genes, c("Symbol", "logFC")]
    
# join the two dataframes into one 
common_up_fc = full_join(case1_overlap, case2_overlap, by="Symbol" )

# set row names as gene names 
rownames(common_up_fc) = common_up_fc$Symbol

# get only the log fold change columns and prepare the matrix
common_up_fc = common_up_fc[,c('logFC.x', 'logFC.y')]
common_up_fc_matrix = as.matrix(common_up_fc)
colnames(common_up_fc_matrix) = c("case 1", "case 2")

# order matrix by per-gene logFC average
matrix_order = order(rowMeans(common_up_fc_matrix))
common_up_fc_matrix = common_up_fc_matrix[ matrix_order, ]


# do the plotting, using all GO terms
int_gos = unique(this_enrichment$Term)

plot_heatmap(this_enrichment, 
                    common_up_fc_matrix, 
                    tfeb_genes, 
                    int_gos, 
                    "red", 
                    "heatmap_union_all_gos.pdf",
                    "genes_union_all_gos.tsv",
                    width=8, height=35)

plot_heatmap(this_enrichment, 
             common_up_fc_matrix, 
             tfeb_genes, 
             int_gos, 
             "black", 
             "heatmap_union_all_gos_black.pdf",
             "genes_union_all_gos_black.tsv",
             width=8, height=35)


# do the plotting, using selected GO terms
int_gos = c("lysosome (GO:0005764)", 
            "lytic vacuole (GO:0000323)",
            "lytic vacuole membrane (GO:0098852)",
            "lysosomal membrane (GO:0005765)", 
            "autophagosome (GO:0005776)",
            "autophagosome membrane (GO:0000421)")

plot_heatmap(this_enrichment, 
                    common_up_fc_matrix, 
                    tfeb_genes, 
                    int_gos, 
                    "black", 
                    "heatmap_union_functionally_related_gos.pdf",
                    "genes_union_functionally_related_gos.tsv",
                    width=8, height=20)

cases_cols = c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "character", "numeric", "character")
case1 = read.table('control-VS-case1.DEseq2_Method.GeneDiffExp_edit.csv', sep=';', header=T, colClasses=cases_cols)
case2 = read.table('control-VS-case2.DEseq2_Method.GeneDiffExp_edit.csv', sep=';', header=T, colClasses=cases_cols)

# perform enrichment analysis on downregulated genes as well

# identify downregulated genes
case1_down   = case1[ case1$classification == 'Down' ,]
case2_down   = case2[ case2$classification == 'Down' ,]

# list of genes that are downregulated in either experiments
union_genes = unique(c(case1_down$Symbol, case2_down$Symbol))

# run enrichR on such gene list
down_enrichment = enrichr(union_genes, enrichr_dbs)

for (i in range(1, length(enrichr_dbs))) {
    
    this_enrichment = down_enrichment[[i]]
    write.table(this_enrichment, sprintf("enrichment_down_%s.tsv", enrichr_dbs_labels[i]), sep='\t', row.names=FALSE, quote=FALSE)
    
    this_enrichment = this_enrichment[ this_enrichment$P.value <= 0.05, ]
    
    # sort by P-value
    sorted_enrichment = this_enrichment[order(this_enrichment$P.value, decreasing=TRUE),]
    
    # number of genes per category
    gene_count = lapply(sorted_enrichment$Genes, function(x) length(strsplit(x, ';' )[[1]]))
    sorted_enrichment$gene.count = gene_count
    
    # unlist the unlistable to make it compatible with ggplot
    sorted_enrichment$P.value = unlist(sorted_enrichment$P.value)
    sorted_enrichment$Term = unlist(sorted_enrichment$Term)
    sorted_enrichment$gene.count = unlist(sorted_enrichment$gene.count)
    
    sorted_enrichment$Term = factor(sorted_enrichment$Term, levels=sorted_enrichment$Term)
    
    # do the plotting of GO terms, number of genes and P-pvalues
    plot = ggplot(data=sorted_enrichment, aes_string(x="P.value", y="Term", size="gene.count")) + 
        geom_point() + 
        scale_size(range=c(0, 6)) + 
        ylab(NULL) +
        scale_x_continuous( trans = "log10", 
                            breaks = trans_breaks("log10", function(x) 10^x),
                            labels = trans_format("log10", math_format(10^.x)))
    ggsave(sprintf("categories_down_%s.pdf", enrichr_dbs_labels[i]), plot, width=10, height=dim(sorted_enrichment)[1]/10, limitsize=FALSE)

    df_l = dim(sorted_enrichment)[1]
    plot = ggplot(data=sorted_enrichment[(df_l-9):df_l,], aes_string(x="P.value", y="Term", size="gene.count")) + 
        geom_point() + 
        scale_size(range=c(0, 6)) + 
        ylab(NULL) +
        scale_x_continuous( trans = "log10", 
                            breaks = trans_breaks("log10", function(x) 10^x),
                            labels = trans_format("log10", math_format(10^.x)))
    ggsave(sprintf("categories_down_first_10_%s.pdf", enrichr_dbs_labels[i]), plot, width=10, height=3.5)
}
