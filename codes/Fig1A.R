# Plot Fig 1A showing the prevalence of antibiotic-resistant organisms (AROs) in regional nursing facilities due to endemic spread of epidemic lineages. 
# load dna matrix data (MGE removed)

# Fig 1B and 1C are a sample of Supplemental Figure S6. Please refer to FigS6.R for details.

library(ape)
library(wesanderson)

mrsa = read.dna("2019-07-19_Filtered_Pathways_MRSA", format = "interleaved")
vrefm = read.dna("2019-07-19_Filtered_Pathways_VREfm", format = "interleaved")
vrefc = read.dna("2019-07-19_Filtered_Pathways_VREfc", format = "interleaved")
ecol = read.dna("2019-07-19_Filtered_Pathways_CipREc", format = "interleaved")
path_aro_clusters = list("MRSA" = mrsa, "VREfm" = vrefm, "VREfc" = vrefc, "CipREc" = ecol)

path_facils = as.character(1:6)

total_enrol_pt_per_facil = structure(c(133, 82, 169, 137, 55, 76), names = 1:6)
  
path_aro_cluster_info = sapply(names(path_aro_clusters), FUN = function(x){
  t(sapply(rownames(path_aro_clusters[[x]]), FUN = function(s){
    info = strsplit(as.character(s), "-")[[1]];
    id = info[1]
    pt = info[2]
    fac = substr(pt, 1, 1)
    visit = info[3]
    site = substr(info[4], 1, 1)
    cbind(paste(id, pt, visit, site, sep = "-"), id, pt, visit, site, fac)
  }))
})

# ARO distribution by NH
aro_prev_by_fac = sapply(path_aro_clusters, FUN = function(y){
  pt = sapply(rownames(y), FUN = function(x){strsplit(x, "-")[[1]][2]})
  fac = substr(pt, 1, 1)
  table(fac)/total_enrol_pt_per_facil
})


# Plot barplot showing prevalence by facility
nh_cols = wes_palette(name = "Darjeeling1", length(path_aro_clusters))
nh_cols[nh_cols %in% "#F2AD00"] = "lightgrey"
par(mar = c(5.1, 8, 4.1, 4.1))
barplot(t(aro_prev_by_fac), 
        col = nh_cols, 
        border="white", 
        font.axis=1, 
        cex.lab = 2, 
        cex.axis = 2,
        cex.names = 2,
        beside= T, 
        xlab="Nursing Facility", 
        font.lab=1, 
        ylab = "Percentage of\ncolonized patients", 
        yaxt = "n",
        ylim = c(0, 0.3))
axis(side = 2, font = 1, at = seq(0, .3, 0.15),
     labels = seq(0, .3*100, 15), las = 2, cex.axis = 2, cex.lab = 2)
legend("top", legend = colnames(aro_prev_by_fac), 
       fill = nh_cols, col = nh_cols, xpd = TRUE, inset = -0.2, text.font = 2,
       cex = 2, bty = "n", horiz = TRUE, border = NA, text.width = 3, x.intersp = 0.25)
