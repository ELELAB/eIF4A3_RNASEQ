# run.R - analyze differential expression analysis data
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

# load cancer types
cancer_types = read.table('../cancer_types.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE)

# filter out types without subtype definition
cancer_types = cancer_types[! is.na(cancer_types$cancer.subtype), ]
print(cancer_types)

# directory where data is stored
data_repository = "../../../TCGA_data/DEA/cancer_subtypes/"

# gene of interest
gene = 'EIF4A3'

# prepare output data frame
data = data.frame(cancer.type = character(),
				  cancer.subtype = character(),
				  logFC = numeric(),
				  AveExpr = numeric(),
				  t = numeric(),
				  P.Value = numeric(),
				  adj.P.Val = numeric(),
				  B = numeric())

main_types = unique(cancer_types$cancer.type)

# for each cancer type
for (t in main_types) {

	subtypes = cancer_types[cancer_types$cancer.type == t,]

	for (st in subtypes$cancer.subtype) {	
		load(paste0(data_repository, t, "/", t, "_", st, "_dataDEGs.rda"))
                print(dim(dataDEGs))

		# get data relative to eIF4A3
		this_data = dataDEGs[rownames(dataDEGs) == gene, ]

		# if the gene is actually found add it to the data frame, keeping
		# information on the cancer type
		if (dim(this_data)[[1]] == 1) {
			this_data$cancer.type = t
			this_data$cancer.subtype = st
			data = rbind(data, this_data) 
		}
		else if (dim(this_data)[[1]] > 1) {
			print(paste("gene", gene, "has more than one entry for subtype", st))
		} 
		else {
			print(paste("gene", gene, "not found for subtype", st))
		}
	}
}

# write output
write.table(data, "eIF4A3_DEA_per_subtype.tsv", sep='\t', quote=FALSE, row.names=FALSE)
