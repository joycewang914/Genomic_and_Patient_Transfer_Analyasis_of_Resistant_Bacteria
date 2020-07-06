# Plot Fig3 
# A) correlation between patient transfer pattern and ARO population structure between nursing facilities
# B) correlation between geographical distance and ARO population structure between nursing facilities

library(wesanderson)

path_fac_dist = readRDS("NF_geocode.RDS")
inter_NH_Fsp_list = readRDS("pruned_inter_NH_Fsp_list.RDS")
all_pt_transfer_dist = readRDS("patient_transfer_KL_distance.RDS")
pw_facil_mat = readRDS("ARO_facil_counts.RDS") # at least 5 isolates in each facility
geo_mat_names = outer(gsub("NF ", "", rownames(all_pt_transfer_dist)), gsub("NF ", "", colnames(all_pt_transfer_dist)), Vectorize(FUN = function(x, y){paste0(x, '_', y)}))

cols = structure(wes_palette(name = "Darjeeling1", n = length(inter_NH_Fsp_list)), names = names(inter_NH_Fsp_list)) 
cols[cols %in% "#F2AD00"] = "darkgrey"

# A) 
par(mar = c(5, 5, 5, 8), xpd = F)
plot(0:1, 0:1,
     ylab = "Isolate similarity\nbetween NFs", 
     xlab = "Patient sharing between NFs", col = "transparent",
     ylim = c(0,1), xlim = c(0, 1.1))

for (org in names(inter_NH_Fsp_list)){
  
  Fsp_mat = inter_NH_Fsp_list[[org]]
  facil_counts = pw_facil_mat[[org]]
  
  final_df = data.frame("geodist" = path_fac_dist[upper.tri(path_fac_dist)],
                        "pt_transfer" = all_pt_transfer_dist[upper.tri(all_pt_transfer_dist)],
                        "fsp" = Fsp_mat[upper.tri(Fsp_mat)],
                        "pair" = geo_mat_names[upper.tri(geo_mat_names)])
  final_df = final_df[facil_counts[upper.tri(facil_counts)] > 0, ]
  
  pt_model = lm(final_df$fsp ~ final_df$pt_transfer)
  
  lines(x = final_df$pt_transfer[grep("1_", final_df$pair, invert = TRUE)], 
        y = final_df$fsp[grep("1_", final_df$pair, invert = T)], type = "p",  pch = 16, col = cols[org])
  lines(x = final_df$pt_transfer[grep("1_", final_df$pair)], 
        y = final_df$fsp[grep("1_", final_df$pair)], type = "p",  pch = 21, col = cols[org])
  
  abline(pt_model, col = cols[org], lwd = 3)
  
  cor_dat = cor.test(final_df$fsp, final_df$pt_transfer, method = "spearman")
  text(1.3, coef(summary(pt_model))[2,1]*1+coef(summary(pt_model))[1,1]*1.1,
       paste0(org, " rho = ", round(cor_dat$estimate, 3), " (p = ", round(cor_dat$p.value, 3), ")"), xpd = TRUE, 
       col = cols[org], cex = .75)  
}

text(0, 1.1, pos = 2, labels = expression(bold("A")), cex = 1.5, xpd = TRUE)

# B
par(mar = c(5, 5, 5, 8), xpd = F)
plot(0:1, 0:1,
     ylab = "Isolate similarity\nbetween NFs", 
     xlab = "Geographical distance\nbetween NFs (kilometers)", col = "transparent",
     ylim = c(0,1), xlim = c(0, 40))

for (org in names(inter_NH_Fsp_list)){
  
  Fsp_mat = inter_NH_Fsp_list[[org]]
  facil_counts = pw_facil_mat[[org]]
  
  final_df = data.frame("geodist" = path_fac_dist[upper.tri(path_fac_dist)],
                        "pt_transfer" = all_pt_transfer_dist[upper.tri(all_pt_transfer_dist)],
                        "fsp" = Fsp_mat[upper.tri(Fsp_mat)],
                        "pair" = geo_mat_names[upper.tri(geo_mat_names)])
  final_df = final_df[facil_counts[upper.tri(facil_counts)] > 0, ]
  
  geodist_model = lm(final_df$fsp ~ final_df$geodist)
  
  lines(x = final_df$geodist[grep("1_", final_df$pair, invert = TRUE)], 
        y = final_df$fsp[grep("1_", final_df$pair, invert = T)], type = "p",  pch = 16, col = cols[org])
  lines(x = final_df$geodist[grep("1_", final_df$pair)], 
        y = final_df$fsp[grep("1_", final_df$pair)], type = "p",  pch = 21, col = cols[org])
  
  abline(geodist_model, col = cols[org], lwd = 3)
  
  cor_dat = cor.test(final_df$fsp, final_df$geodist, method = "spearman")
  text(52, coef(summary(geodist_model))[2,1]*40+coef(summary(geodist_model))[1,1]*1.1,
       paste0(org, " rho = ", round(cor_dat$estimate, 3), " (p = ", round(cor_dat$p.value, 3), ")"), xpd = TRUE, 
       col = cols[org], cex = .75)  
}

text(0, 1.1, pos = 2, labels = expression(bold("B")), cex = 1.5, xpd = TRUE)

