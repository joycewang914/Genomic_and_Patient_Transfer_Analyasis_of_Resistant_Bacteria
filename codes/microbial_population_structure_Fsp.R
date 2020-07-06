# Calculate SNV-based gene flow between NFs using an adaptation of Wright F statistic (from Donker et al, Microbial Genomics 2017).
# e.g. if the genetic distance of bacterial isolates found within the same nursing facility is equal to that of between two facilities, the two populations are homogenous and Fsp = 0; if two populations are completely different (i.e. no gene flow) then Fsp = 1.

# Only use samples included in pruned phylogenetic trees (see FigS6).

# For Fsp analysis, only facilities with > 5 isolates were included.


library(ape)
source("Fsp_calculation.R")
path_facils = as.character(1:6)

# Read in SNV matrix
mrsa = read.dna("2019-07-19_Filtered_Pathways_MRSA", format = "interleaved")
vrefm = read.dna("2019-07-19_Filtered_Pathways_VREfm", format = "interleaved")
vrefc = read.dna("2019-07-19_Filtered_Pathways_VREfc", format = "interleaved")
ecol = read.dna("2019-07-19_Filtered_Pathways_CipREc", format = "interleaved")
path_aro_clusters = list("MRSA" = mrsa, "VREfm" = vrefm, "VREfc" = vrefc, "CipREc" = ecol)
path_facils = as.character(1:6)

# Read in tree where putative intra-facility transmission isolates have been removed
org_tree = readRDS("all_vs_pruned_trees.RDS")

# Fsp calculation
inter_NH_Fsp_list = list()
pw_facil_mat = list()

for (org in names(path_aro_clusters)){

  dna = path_aro_clusters[[org]]
  dna = dna[,apply(as.character(dna), 2, FUN = function(x){
    sum(x %in% c("-", "n", "N")) == 0
  })]
  
  temp_mat = dna[grep(paste(org_tree[[org]][[2]]$tip.label, collapse = "|"), rownames(dna)),]

  temp_fsp = facility_fsp(temp_mat, path_facils)
  inter_NH_Fsp_list[[org]] = temp_fsp
  
  # Count how many facility pairs have at least 5 isolates of each species
  facil_count_table = table(sapply(rownames(temp_mat), 
                                   FUN = function(x){substr(strsplit(x, "-")[[1]][2], 1, 1)}))
  facil_counts = outer(facil_count_table, facil_count_table,
                       Vectorize(FUN = function(y,z){
                         y >= 5 & z >= 5
                       })) 
  
  pw_facil_mat[[org]] = facil_counts
}
saveRDS(file = "pruned_inter_NH_Fsp_list.RDS", object = inter_NH_Fsp_list)
saveRDS(file = "ARO_facil_counts.RDS", object = pw_facil_mat)
