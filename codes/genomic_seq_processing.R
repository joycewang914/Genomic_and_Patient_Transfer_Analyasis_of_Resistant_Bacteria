# This script lists steps for generating final single-nucleotide variant matrices used for phylogenetic and population structure analyses.

# Working witih 4 common antibiotic-resistant organisms (AROs) in Michigan nursing homes:
# MRSA (mrsa): methicillin-resistant Staphylococcus aureus
# Escherichia coli (ecol): ciprofloxacin-resistant E. coli
# Vancomycin-resistant Enterococcus faecalis (vrefc): Vancomycin-resistant Enterococcus faecalis
# Vancomycin-resistant Enterococcus faecium (vrefm): Vancomycin-resistant Enterococcus faecium

library(ape)

# Read in gubbin-filtered SNV files
## Load fasta files

mrsa_dna = read.dna("2019_03_19_14_07_22_803-A001-120-N-MRSA-TIP__ordered_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta", format = "fasta")
ecol_dna = read.dna("2019_03_19_10_08_45_3425-5015-0-R2X_ordered_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta",  format = "fasta")
vrefc_dna = read.dna("2019_03_19_14_29_34_7605-4117-0-HVRE_ordered_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta",  format = "fasta")
vrefm_dna = read.dna("2019_03_19_18_03_26_3399-4058-0-RVRE_final_ordered_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites.fasta",  format = "fasta")
mdro_dnas = list("mrsa" = mrsa_dna, "ecol" = ecol_dna, "vrefc" = vrefc_dna, "vrefm" = vrefm_dna)

# Remove E. coli isolates on a long branch (genetically divergent from the rest of collection)
mdro_dnas$ecol = mdro_dnas$ecol[grep("6044|6048|6057", rownames(mdro_dnas$ecol), invert = TRUE),]

# Only analyzing samples from a recent collection (Mody et al, CID 2018)
# Facilities denoted as "path_facils" (Pathways study facilities)
path_facils = 1:6
path_dnas = lapply(mdro_dnas, FUN = function(x){
  x[sapply(rownames(x), FUN = function(y){substr(strsplit(y, "-")[[1]][2], 1, 1) %in% path_facils}),]
})

## Read in SNV positions
mrsa_snv = read.table("2019_03_19_14_07_22_803-A001-120-N-MRSA-TIP__ordered_genome_aln_w_alt_allele_unmapped_gubbins.summary_of_snp_distribution.vcf", header = TRUE, skip = 3)
vrefm_snv = read.table("2019_03_19_18_03_26_3399-4058-0-RVRE_final_ordered_genome_aln_w_alt_allele_unmapped_gubbins.summary_of_snp_distribution.vcf", header = TRUE, skip = 3)
vrefc_snv = read.table("2019_03_19_14_29_34_7605-4117-0-HVRE_ordered_genome_aln_w_alt_allele_unmapped_gubbins.summary_of_snp_distribution.vcf", header = TRUE, skip = 3)
ecol_snv = read.table("2019_03_19_10_08_45_3425-5015-0-R2X_ordered_genome_aln_w_alt_allele_unmapped_gubbins.summary_of_snp_distribution.vcf", header = TRUE, skip = 3)

snvs = list("mrsa" = mrsa_snv, "vrefm" = vrefm_snv, "vrefc" = vrefc_snv, "ecol" = ecol_snv)

# Read in RAST file (genome annotation) to screen for regions associated with horizontal gene transfer
mrsa_rast = read.table("803-A001-120-N-MRSA.gff", skip = 1, sep = "\t", quote = "")
vre_fm_rast = read.table("3399-4058-0-RVRE.gff", skip = 1, sep = "\t", quote = "")
vre_fc_rast = read.table("7605-4117-0-HVRE.gff", skip = 1, sep = "\t", quote = "")
ecol_rast = read.table("3425-5015-0-R2X.gff", skip = 1, sep = "\t", quote = "")
rasts =  list("mrsa" = mrsa_rast, "vrefm" = vre_fm_rast, "vrefc" = vre_fc_rast, "ecol" = ecol_rast)

HGT_patterns = c("mobile","insertion","transposon","transposase","phage","conjugative")

rast_annot_list = list()
for (org in names(rasts)){print(org)
  
  rast_annot = rasts[[org]]
  rast_annot_ind = grep(paste(HGT_patterns, collapse = "|"), as.character(rast_annot$V9), ignore.case = TRUE)
  HGT_pos = unique(unlist(apply(rast_annot[rast_annot_ind, ], 1, FUN = function(x){x[4]:x[5]})))
  rast_annot_list[[org]] = HGT_pos
}

# Remove SNVs within HGT regions from fasta file

filter_path_dna = list()
for (org in names(rast_annot_list)){print(org)
  dna = path_dnas[[org]]
  snv = snvs[[org]]
  snv = snv[,2:ncol(snv)]
  rast_annot = rast_annot_list[[org]]
  
  snv_filter = which(snv$POS %in% intersect(snv$POS, rast_annot))
  print(paste0(length(snv_filter), " positions filtered for HGT-associated annotations"))
  dna_filter = dna[, -snv_filter]
  
  filter_path_dna[[org]] = dna_filter
}


write.dna(filter_path_dna$mrsa, file = "2019-07-19_Filtered_Pathways_MRSA", format = "interleaved")
write.dna(filter_path_dna$vrefm, file = "2019-07-19_Filtered_Pathways_VREfm", format = "interleaved")
write.dna(filter_path_dna$vrefc, file = "2019-07-19_Filtered_Pathways_VREfc", format = "interleaved")
write.dna(filter_path_dna$ecol, file = "2019-07-19_Filtered_Pathways_CipREc", format = "interleaved")