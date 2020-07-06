# Plot FigS1 visualizing time of isolation of each antibiotic-resistant organism.

library(ape)
library(wesanderson)

dna_mrsa = read.dna("2019-07-19_Filtered_Pathways_MRSA", format = "interleaved")
dna_vrefm = read.dna("2019-07-19_Filtered_Pathways_VREfm", format = "interleaved")
dna_vrefc = read.dna("2019-07-19_Filtered_Pathways_VREfc", format = "interleaved")
dna_ecol = read.dna("2019-07-19_Filtered_Pathways_CipREc", format = "interleaved")
path_aro_clusters = list("MRSA" = dna_mrsa, 
                         "VREfm" = dna_vrefm, 
                         "VREfc" = dna_vrefc, 
                         "CipREc" = dna_ecol)
path_facils = as.character(1:6)

cols = structure(wes_palette(name = "Darjeeling1", n = length(path_aro_clusters)), 
                 names = names(path_aro_clusters)) 
cols[cols %in% "#F2AD00"] = "darkgrey"

par(mfrow = c(2, 2))
for (org in names(path_aro_clusters)){print(org)
  
  temp_mat = structure(rep(0, 8), names = c("0", "14","30", "60", "90","120", "150" ,"180"))
  visits = sapply(rownames(path_aro_clusters[[org]]), FUN = function(x){
    as.numeric(as.character(strsplit(x, "-")[[1]][3]))
  }) 
  col = cols[org]
  temp_mat[names(table(visits))] = table(visits)
  bp = barplot(temp_mat, col = col, main = org, border = FALSE, xlab = "Day since enrollment",
               ylab = "Frequency", xaxt = "n", plot = TRUE)
  text(x = bp,y = -max(temp_mat)/7, names(temp_mat), xpd = TRUE, srt = 90)
  
}

