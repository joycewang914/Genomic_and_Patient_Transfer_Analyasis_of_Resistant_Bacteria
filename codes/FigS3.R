# Fig S3 showing the distribution and frequency of pair-wise SNV distance within each antibiotic-resistant bacterial species. 

library(wesanderson)
library(ape)

mrsa = read.dna("2019-07-19_Filtered_Pathways_MRSA", format = "interleaved")
vrefm = read.dna("2019-07-19_Filtered_Pathways_VREfm", format = "interleaved")
vrefc = read.dna("2019-07-19_Filtered_Pathways_VREfc", format = "interleaved")
ecol = read.dna("2019-07-19_Filtered_Pathways_CipREc", format = "interleaved")
path_aro_clusters = list("MRSA" = mrsa, "VREfm" = vrefm, "VREfc" = vrefc, "CipREc" = ecol)
path_facils = as.character(1:6)

aro_dna_dist = list()
for (o in names(path_aro_clusters)){
  
  dna_mat = path_aro_clusters[[o]][,apply(as.character(as.matrix(path_aro_clusters[[o]])), 2, FUN = function(x){
    sum(x %in% c("-", "n")) == 0
  })]
  dna_dist = dist.dna(dna_mat, as.matrix = TRUE)*ncol(dna_mat)
  
  aro_dna_dist[[o]] = dna_dist
}

cols = structure(wes_palette(name = "Darjeeling1", 
                             length(aro_dna_dist) + 1, type = "continuous")[1:4],
                 names = names(aro_dna_dist))
cols[cols %in% "#F2AD00"] = "lightgrey"

# Plot histograms
par(mfrow = c(2, 2))
for (o in names(aro_dna_dist)){print(o)
  
  max_snv = round(ceiling(max(aro_dna_dist[[o]]))*1.1, digits = -1)
  hist(aro_dna_dist[[o]][upper.tri(aro_dna_dist[[o]])], 
       xlim = c(0, max_snv), xaxt = "n", 
       breaks = seq(0, max_snv, 5), col = cols[o], 
       main = o, xlab = "Pairwise SNV distance", lty = "blank", ylab = "Isolate pair frequency")
  axis(side = 1, seq(0, max_snv, 5), labels = seq(0, max_snv, 5), 
       cex.axis = 0.75, tick = TRUE, lwd.ticks = 0.05)
}