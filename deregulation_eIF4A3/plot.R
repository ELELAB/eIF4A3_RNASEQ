# plot.R - plot differential expression analysis data
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

library(ggplot2)
library(RColorBrewer)

# read subtype data
data = read.table("eIF4A3_DEA_per_subtype.tsv", sep='\t', header=TRUE, stringsAsFactors=FALSE)

# generate plot labels
data$labels = paste(data$cancer.type, data$cancer.subtype, sep=', ')

# prepare colors
n.types = length(unique(data$cancer.type))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
colors=getPalette(n.types)

# plot
png('deregulation_eIF4A3_per_subtype.png')
ggplot(data=data, aes(x=labels, y=logFC, fill=cancer.type)) +
geom_bar(stat="identity") +
geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
geom_hline(yintercept=-0.5, linetype="dashed", color = "red") +
labs(y="logFC", x="Cancer subtype", fill="Cancer type") +
theme_linedraw() + 
theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
scale_fill_manual(values=colors) +
coord_flip()
dev.off()

# read type data
data.type = read.table("eIF4A3_DEA_per_type.tsv", sep='\t', header=TRUE, stringsAsFactors=FALSE)

# keep only types for which we have subtypes
data.type = data.type[ data.type$cancer.type %in% data$cancer.type, ]

# plot
png('deregulation_eIF4A3_per_type.png')
ggplot(data=data.type, aes(x=cancer.type, y=logFC, fill=cancer.type)) +
geom_bar(stat="identity", show.legend=FALSE) +
geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
geom_hline(yintercept=-0.5, linetype="dashed", color = "red") +
labs(y="logFC", x="Cancer type") +
theme_linedraw() + 
theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
scale_fill_manual(values=colors) +
coord_flip()
dev.off()
