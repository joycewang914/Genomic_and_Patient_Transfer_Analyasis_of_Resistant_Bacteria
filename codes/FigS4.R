# Plot FigS4 to visualize the homogeneity of antibiotic-resitant lineages in southeast Michigan over time (2010 - 2015). Blue indicates isolates from an earlier study (Mody 2015 JAMA IM) and red indicates isolates from current analysis (Mody 2018 CID).

library(gdata)
library(ape)
library(pheatmap)

# load DNA data from current analysis
mrsa = read.dna("../../../../CDC_grant/Year1/Analysis/2019_06_13_Pathways_genomic_analysis/2019-07-19_Filtered_Pathways_MRSA", format = "interleaved")
vrefm = read.dna("../../../../CDC_grant/Year1/Analysis/2019_06_13_Pathways_genomic_analysis/2019-07-19_Filtered_Pathways_VREfm", format = "interleaved")
vrefc = read.dna("../../../../CDC_grant/Year1/Analysis/2019_06_13_Pathways_genomic_analysis/2019-07-19_Filtered_Pathways_VREfc", format = "interleaved")
ecol = read.dna("../../../../CDC_grant/Year1/Analysis/2019_06_13_Pathways_genomic_analysis/2019-07-19_Filtered_Pathways_CipREc", format = "interleaved")
mdro_clusters = list("MRSA" = mrsa, "VREfm" = vrefm, "VREfc" = vrefc, "CipREc" = ecol)

# load DNA data from both current and previous studies
# Load DNA data (all isolates in one single cluster)
mrsa_dna = read.dna("single_cluster/2019_04_10_17_09_08_MRSA_USA_100_1_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta",  format = "fasta")
ecol_dna = read.dna("single_cluster/2019_04_10_17_09_08_EC958_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta",  format = "fasta")
vrefc_dna = read.dna("single_cluster/2019_04_10_11_04_35_Efaecalis_V583_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta",  format = "fasta")
vrefm_dna = read.dna("single_cluster/2019_03_19_18_03_26_3399-4058-0-RVRE_final_ordered_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta",  format = "fasta")
mdro_single_clusters = list("MRSA" = mrsa_dna, "VREfm" = vrefm_dna, "VREfc" = vrefc_dna, "CipREc" = ecol_dna)

# load MLST data (remove isolates from same patient)
mlst_ecoli = read.table("CipREc_MLST.txt"); 
mlst_ecoli = mlst_ecoli[!rownames(mlst_ecoli) %in% c("222-K008-30-GX-Ecoli-TIP_S27_L001",
                                                     "3038-F035-60-OX2-Ecoli-TIP",
                                                     "3099-F035-90-OX2-Ecoli-TIP_S14_L001",
                                                     "3872-L020-1-G-TIPmvre_S35_L001",
                                                     "442-K008-60-RX-Ecoli-TIP_S41_L001",
                                                     "466-F011-1-GX1-Ecoli-TIP",
                                                     "519-F001-30-RX1-Ecoli-TIP_S14_L001", 
                                                     "829-F012-60-GX1-Ecoli-TIP_S9_L001",
                                                     "651-F012-30-FX-Ecoli-TIP_S42_L001",
                                                     "721-F001-60-WX1-Ecoli-TIP_S1_L001" ), ]

mlst_saureus = read.table("MRSA_MLST.txt"); 
mlst_saureus = mlst_saureus[grep("34-A001-1-F-MRSA-TIP_S6_L001", rownames(mlst_saureus), invert = T), ]

mlst_efc = read.csv("VREfc_MLST.csv")
rownames(mlst_efc) = mlst_efc$X
mlst_efc = mlst_efc[!rownames(mlst_efc) %in% c("3366-F050-14-G-VRE-TIP_S7_L001", 
                                               "3447-F050-30-G-VRE-TIP_S15_L001",
                                               "3448-F050-30-R-VRE-TIP_S19_L001", 
                                               "513-E013-14-G-VRE-TIP_S55_L001",
                                               "682-E013-DC-G-VRE-TIP_S49_L001"), ] 
mlst_efm = read.table("VREfm_MLST.txt"); 
mlst_efm = mlst_efm[!rownames(mlst_efm) %in% c("702-3020-0-RVRE_S18_L001", 
                                               "546-4011-0-RVRE_S58_L001",
                                               "3362-E077-30-R-TIPmvre_S22_L001", 
                                               "3363-A015-14-R-TIPmvre_S29_L001"), ]

mlst_clusters = list("MRSA" = mlst_saureus, "VREfm" = mlst_efm, "VREfc" = mlst_efc, "CipREc" = mlst_ecoli)

path_facils = as.character(1:6) # current study
TIP = LETTERS[1:12] # previous study

# Single cluster

# Some STs are similar to the ones we analyzed - use Pathways assignment to help find the equivalent STs in TIP collection
single_cluster_mlsts = sapply(names(mlst_clusters), FUN = function(o){print(o)
  
  # Pathways patients only
  mdro_pt = sapply(rownames(mdro_clusters[[o]]), FUN = function(x){strsplit(x, "-")[[1]][2]})
  
  # Pathways MLST
  mdro_pt_mlst = sapply(mdro_pt, FUN = function(x){
    as.character(mlst_clusters[[o]][grep(paste0("-", x, "-"), rownames(mlst_clusters[[o]])), "ST"])
  })
  
  # Patients in Pathways and TIP with prevalent ST as per Pathways data (but remove ND and Novel)
  tip_pts = unlist(sapply(rownames(mlst_clusters[[o]][(mlst_clusters[[o]]$ST %in% unique(mdro_pt_mlst)) &
                            !(mlst_clusters[[o]]$ST %in% c("ND", "Novel")),]), FUN = function(x){
    strsplit(x, "-")[[1]][2][substr(strsplit(x, "-")[[1]][2], 1, 1) %in% TIP]
  }))
  
  path_tip_pts = c(mdro_pt, tip_pts)
  
  # Make sure we use samples that we have both DNA and MLST data
  common_pts = path_tip_pts[sapply(path_tip_pts, FUN = function(x){
    length(grep(paste0("-", x, "-"), rownames(mdro_single_clusters[[o]]))) > 0
  })]
  
  mlst = sapply(common_pts, FUN = function(p){
    as.character(mlst_clusters[[o]][grep(paste0("-", p, "-"), rownames(mlst_clusters[[o]])), "ST"])
    })

  prev_mlst = names(which.max(sort(table(mlst))))
  
  list(pts = common_pts, st = prev_mlst)
})


# Make heatmaps showing only the most prevalent ST and label by pathways versus TIP
library(pheatmap)
for (o in names(mdro_clusters)){print(o)
  st_pts_full = single_cluster_mlsts[,o][["pts"]]

  st_pts_ind = sapply(names(st_pts_full), FUN = function(x){
    paste(strsplit(x, "-")[[1]][1:3], collapse = "-")
  })
  
  dna_dist = dist.dna(mdro_single_clusters[[o]][grep(paste(st_pts_ind, collapse = "|"), 
          rownames(mdro_single_clusters[[o]])),], pairwise.deletion = FALSE, model = "N", as.matrix = T)
  
  study_ind = sapply(rownames(dna_dist), FUN = function(x){
    ifelse(substr(strsplit(x, "-")[[1]][2], 1, 1) %in% path_facils, "Current", "Past")
  })

  annotdf <- data.frame(row.names = names(study_ind), 
                        Study = study_ind)
  
  newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotdf$Study))))
  mycolors <- newCols(length(unique(annotdf$Study)))
  names(mycolors) <- sort(unique(annotdf$Study))
  mycolors <- list(Study = mycolors)
  
  st = single_cluster_mlsts[,o]$st

  pheatmap(dna_dist,
           scale="none",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           border_color=NA,
           main=paste(ifelse(o %in% "VREfm", paste0(o, "(clade A)"), paste0(o, " (ST ", st, " and like)"))), # VREfc mostly belongs to ST 6 - edited figure afterwards
           annotation_row = annotdf,
           annotation_col = annotdf,
           annotation_names_col = F,
           annotation_names_row = F,
           cellwidth = 2,
           cellheight = 2,
           annotation_colors = mycolors,legend = TRUE,show_rownames = F, 
           show_colnames = F,fontsize = 10
           )

}