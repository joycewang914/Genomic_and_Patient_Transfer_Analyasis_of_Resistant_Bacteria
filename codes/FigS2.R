# Fig S2 showing the proportion of sequence types (ST) within each species. 
# Include all isolates 

library(wesanderson)
library(ape)

# Load DNA data (all isolates in one single cluster)
mrsa_dna = read.dna("single_cluster/2019_04_10_17_09_08_MRSA_USA_100_1_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta",  format = "fasta")
ecol_dna = read.dna("single_cluster/2019_04_10_17_09_08_EC958_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta",  format = "fasta")
vrefc_dna = read.dna("single_cluster/2019_04_10_11_04_35_Efaecalis_V583_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta",  format = "fasta")
vrefm_dna = read.dna("single_cluster/2019_03_19_18_03_26_3399-4058-0-RVRE_final_ordered_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta",  format = "fasta")
mdro_single_clusters = list("MRSA" = mrsa_dna, "VREfm" = vrefm_dna, "VREfc" = vrefc_dna, "CipREc" = ecol_dna)

# Load DNA data
mrsa = read.dna("2019-07-19_Filtered_Pathways_MRSA", format = "interleaved")
vrefm = read.dna("2019-07-19_Filtered_Pathways_VREfm", format = "interleaved")
vrefc = read.dna("2019-07-19_Filtered_Pathways_VREfc", format = "interleaved")
ecol = read.dna("2019-07-19_Filtered_Pathways_CipREc", format = "interleaved")
path_aro_clusters = list("MRSA" = mrsa, "VREfm" = vrefm, "VREfc" = vrefc, "CipREc" = ecol)

# Load MLST data
mlst_ecoli = read.table("CipREc_MLST.txt"); 
mlst_saureus = read.table("MRSA_MLST.txt"); 
mlst_efc = read.csv("VREfc_MLST.csv")
rownames(mlst_efc) = mlst_efc$X
mlst_efm = read.table("VREfm_MLST.txt"); 
mlst_efm = mlst_efm[!rownames(mlst_efm) %in% c("702-3020-0-RVRE_S18_L001", "546-4011-0-RVRE_S58_L001"), ]
mlst_clusters = list("MRSA" = mlst_saureus, "VREfm" = mlst_efm, "VREfc" = mlst_efc, "CipREc" = mlst_ecoli)

path_facils = as.character(1:6)

# Focus on isolates from Mody CID 2018:
# All isolates
path_single_cluster_mlsts = sapply(names(mlst_clusters), FUN = function(o){print(o)
  mdro_pt = sapply(rownames(mdro_single_clusters[[o]]), FUN = function(x){strsplit(x, "-")[[1]][2]})
  path_mdro_pt = mdro_pt[substr(mdro_pt, 1, 1) %in% path_facils]
  mlst_pt = rownames(mlst_clusters[[o]])[sapply(rownames(mlst_clusters[[o]]), 
                                                FUN = function(x){substr(strsplit(x, "-")[[1]][2],1, 1) %in% path_facils})]
  
  path_mlst_pt = sapply(mlst_pt, FUN = function(x){strsplit(x, "-")[[1]][2]})
  common_pt = intersect(mdro_pt, path_mlst_pt); print(length(common_pt))
  mlst = sapply(common_pt, FUN = function(p){
    as.character(mlst_clusters[[o]][grep(paste0("-", p, "-"), rownames(mlst_clusters[[o]])), "ST"])})
  
  mlst
})

# Analyzed isolates (belonging to most prevalent ST)
path_mlsts = sapply(names(mlst_clusters), FUN = function(o){print(o)
  mdro_pt = sapply(rownames(path_aro_clusters[[o]]), FUN = function(x){strsplit(x, "-")[[1]][2]})
  path_mdro_pt = mdro_pt[substr(mdro_pt, 1, 1) %in% path_facils]
  mlst_pt = rownames(mlst_clusters[[o]])[sapply(rownames(mlst_clusters[[o]]), 
                                                FUN = function(x){substr(strsplit(x, "-")[[1]][2],1, 1) %in% path_facils})]
  
  path_mlst_pt = sapply(mlst_pt, FUN = function(x){strsplit(x, "-")[[1]][2]})
  common_pt = intersect(mdro_pt, path_mlst_pt); print(length(common_pt))
  mlst = sapply(common_pt, FUN = function(p){
    as.character(mlst_clusters[[o]][grep(paste0("-", p, "-"), rownames(mlst_clusters[[o]])), "ST"])})
  
  mlst
})


# Plot MLST barplot
par(mfrow = c(2, 2))
col = rev(wes_palette("Zissou1", 10, type = "continuous"))
for (o in names(path_single_cluster_mlsts)){
  
  all_samples = path_single_cluster_mlsts[[o]]

  if (o %in%  c("VREfm")){all_samples_wrt_common_st = all_samples} else{
    analyzed_samples = path_mlsts[[o]]
    # analyzed_samples = gsub("*", "", analyzed_samples, fixed = TRUE)
    most_common_st = names(which.max(table(analyzed_samples)))
    all_samples_wrt_common_st = sapply(names(all_samples), FUN = function(pt){
      if (pt %in% names(analyzed_samples)){
        if (length(unique(analyzed_samples)) == 1){leg = paste0("ST", most_common_st)}else{
          leg = paste0("ST", most_common_st, " and ST", most_common_st, "-like")}
      }else{
        if (o != "MRSA"){
          leg = paste0("Others")} else{
            most_common_not_analyzed_st = names(which.max(table(all_samples[!names(all_samples) %in% names(analyzed_samples)])))
            leg = paste0("ST", most_common_not_analyzed_st, " and others")
          }
      }
    })
  }
  mlst_table = rev(sort(table(all_samples_wrt_common_st)))/sum(table(all_samples_wrt_common_st))
  mlst_col = col[1:length(mlst_table)]
  
  barplot(mlst_table, col = mlst_col, ylim = c(0, 1), ylab = "", cex.axis = 1, 
          main = o, yaxt = "n", cex.lab = 1, cex.main = 1, cex.names = 1, las = 2, border = FALSE)
  title(ylab = "Proportion (%)", cex.lab = 1, line = 5)
  
  axis(side=2, font=1, at = c(0, 0.25, 0.5, 0.75, 1),
       labels = c("0", "25", "50", "75", "100"), las =2, cex.lab = 1, cex.axis = 1)
}

