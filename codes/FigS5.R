# # Plot Fig S6 using violin plots to show the distribution of single nucleotide variant (SNV) distance among pairs of nearest neighbour isolates between and within nursing facility and collected at admission and during follow-up visits using all isolates and only admission isolates and follow-up isolates that did not have a closely related neighbour within the same facility.

library(hues)
library(vioplot)
library(scales)

# Read in all + pruned trees (generated from FigS6.R)
org_tree = readRDS("all_vs_pruned_trees.RDS")

# Read in DNA data
dna_mrsa = read.dna("2019-07-19_Filtered_Pathways_MRSA", format = "interleaved")
dna_vrefm = read.dna("2019-07-19_Filtered_Pathways_VREfm", format = "interleaved")
dna_vrefc = read.dna("2019-07-19_Filtered_Pathways_VREfc", format = "interleaved")
dna_ecol = read.dna("2019-07-19_Filtered_Pathways_CipREc", format = "interleaved")
path_dna = list("MRSA" = dna_mrsa, "VREfm" = dna_vrefm, "VREfc" = dna_vrefc, "CipREc" = dna_ecol)


# Make vioplots
hist_cols = structure(iwanthue(6), 
                      names = c("Nearest admission isolates between facilities",
                                "Nearest admission isolates within facility",
                                "Nearest follow-up isolates between facilities",
                                "Nearest follow-up isolates within facilities",
                                "Nearest admission isolates to follow-up isolates within facility",
                                "Nearest follow-up isolates to admission isolates between facilities"))

m = matrix(c(1,3,2, 4, 5, 7, 6, 8, 9, 9), nrow = 2, ncol = 5, widths = c(0.5, 0.5, 0.5, 0.5, 0.5))
layout(m)
# Plot violin plots
for (org in names(org_tree)){
  
  tree = org_tree[[org]]
  
  temp_dna = path_dna[[org]]
  total_dna_dist = dist.dna(temp_dna, pairwise.deletion = FALSE, model = "N", as.matrix = T)
  
  
  for (n in names(tree)){print(n)
    
    ind = grep(paste(tree[[n]]$tip.label, collapse = "|"), rownames(temp_dna))
    dna_dist = total_dna_dist[ind, ind]
    
    admission_isolates = rownames(temp_dna)[ind][sapply(rownames(temp_dna)[ind], FUN = function(x){
      (strsplit(x, "-")[[1]][3]) %in% "0"
    })]
    admission_isolates_facs = structure(sapply(admission_isolates, FUN = function(x){
      substr(strsplit(x, "-")[[1]][2], 1, 1)
    }), names = admission_isolates)
    
    acq_isolates = rownames(temp_dna)[ind][sapply(rownames(temp_dna)[ind], FUN = function(x){
      (strsplit(x, "-")[[1]][3]) != "0"
    })]
    acq_isolates_facs = structure(sapply(acq_isolates, FUN = function(x){
      substr(strsplit(x, "-")[[1]][2], 1, 1)
    }), names = acq_isolates)
    
    
    closest_admit_between = sapply(names(admission_isolates_facs), FUN = function(x){
      fac1 = admission_isolates_facs[x]
      min(dna_dist[x, names(admission_isolates_facs[!admission_isolates_facs %in% fac1])])
    })
    
    closest_admit_within = sapply(names(admission_isolates_facs), FUN = function(x){
      fac1 = admission_isolates_facs[x]
      min(dna_dist[x, names(admission_isolates_facs[admission_isolates_facs %in% fac1 & 
                                                      !names(admission_isolates_facs) %in% x])])
    })
    
    closest_acq_between = sapply(names(acq_isolates_facs), FUN = function(x){
      fac1 = acq_isolates_facs[x]
      min(dna_dist[x, names(acq_isolates_facs[!acq_isolates_facs %in% fac1])])
    })
    
    closest_acq_within = sapply(names(acq_isolates_facs), FUN = function(x){
      fac1 = acq_isolates_facs[x]
      if (sum(acq_isolates_facs %in% fac1 & !names(acq_isolates_facs) %in% x) == 0){NA}else{
        min(dna_dist[x, names(acq_isolates_facs[acq_isolates_facs %in% fac1 & 
                                                  !names(acq_isolates_facs) %in% x])])}
    })
    
    closest_acq_admit_between = sapply(names(acq_isolates_facs), FUN = function(x){
      fac1 = acq_isolates_facs[x]
      min(dna_dist[x, names(admission_isolates_facs[!admission_isolates_facs %in% fac1])])
    })
    
    closest_acq_admit_within = sapply(names(acq_isolates_facs), FUN = function(x){
      fac1 = acq_isolates_facs[x]
      min(dna_dist[x, names(admission_isolates_facs[admission_isolates_facs %in% fac1])])
    })
    
    
    
    pw_dist = list(
      "Nearest admission isolates\nbetween facilities\n" = closest_admit_between,
      "Nearest admission isolates\nwithin facility\n" = closest_admit_within,
      "Nearest follow-up isolates\nbetween facilities\n" = closest_acq_between,
      "Nearest follow-up isolates\nwithin facilities\n" = closest_acq_within,
      "Nearest admission to follow-up\nisolates between facilities\n" = closest_acq_admit_between,
      "Nearest admission to follow-up\nisolates within facility\n" = closest_acq_admit_within)
    
    if (n %in% names(tree)[1]){plot_title = paste0(org, "\n", n)}else{plot_title = n}
    
    max_yval = max(total_dna_dist)
    vioplot(pw_dist[[1]], 
            pw_dist[[2]],
            pw_dist[[3]],
            pw_dist[[4]], 
            pw_dist[[5]],
            pw_dist[[6]],border = hist_cols, ylim = c(0, max_yval),
            col = alpha(hist_cols, alpha = 0.25),
            xlab = "Pairwise SNV distance", ylab = "",
            main = plot_title,cex.main = 1.5, names = FALSE, horizontal = T)
    
    for (i in 1:length(pw_dist)){
      temp_type = names(pw_dist)[[i]]
      stripchart(pw_dist[[i]], at = i, add = T, 
                 bg = hist_cols[i], pch = 21, method = "jitter", vertical = F)
    }
    
  }
}
plot(0, 0, type = "n", ann = F, axes = F)
legend("right", legend = rev(c("Nearest admission isolates\nbetween facilities\n",
                               "Nearest admission isolates\nwithin facility\n",
                               "Nearest follow-up isolates\nbetween facilities\n",
                               "Nearest follow-up isolates\nwithin facilities\n",
                               "Nearest admission to follow-up\nisolates between facilities\n",
                               "Nearest admission to follow-up\nisolates within facility\n")),
       fill = rev(hist_cols),col = rev(hist_cols), xpd = TRUE, inset = -0.15,
       border = FALSE, bty = "n",  
       cex = 1.5)
