# Plot Figure S10 showing the correlation between the popoulation structure of each antibiotic-resistant lineage and patient sharing pattern between nursing facilities.


library(wesanderson)

inter_NH_Fsp_list = readRDS("pruned_inter_NH_Fsp_list.RDS")
all_pt_transfer_dist = readRDS("patient_transfer_KL_distance.RDS")
pw_facil_mat = readRDS("ARO_facil_counts.RDS") # at least 5 isolates in each facility
geo_mat_names = outer(gsub("NF ", "", rownames(all_pt_transfer_dist)), gsub("NF ", "", colnames(all_pt_transfer_dist)), Vectorize(FUN = function(x, y){paste0(x, '_', y)}))

# Fsp vs patient transfer
cols = structure(wes_palette(name = "Darjeeling1", n = length(inter_NH_Fsp_list)), names = names(inter_NH_Fsp_list)) 
cols[cols %in% "#F2AD00"] = "darkgrey"

par(mfrow = c(2, 2))

for (org in names(inter_NH_Fsp_list)){
    temp_fsp = inter_NH_Fsp_list[[org]]

    final_df = data.frame("pt_transfer" = all_pt_transfer_dist[upper.tri(all_pt_transfer_dist)],
                          "fsp" = temp_fsp[upper.tri(temp_fsp)],
                          "pair" = geo_mat_names[upper.tri(geo_mat_names)])
    final_df = final_df[pw_facil_mat[[org]][upper.tri(pw_facil_mat[[org]])] > 0, ]
    
    pt_model = lm(final_df$fsp ~ final_df$pt_transfer)
    
    plot(final_df$pt_transfer, final_df$fsp, pch = 16, 
         ylab = "Isolate similarity between NFs\n(low:similar; high:different)", 
         xlab = "Patient sharing between NFs\n(low:similar; high:different)", col= cols[org])
    abline(pt_model, col = cols[org], lwd = 3)
    text(final_df$fsp ~ final_df$pt_transfer, labels=pair,data=final_df, cex=0.75, font=2)
    
    cor_dat = cor.test(final_df$fsp, final_df$pt_transfer, method = "spearman")
    title(paste0(org, "\nSpearman rho = ", round(cor_dat$estimate, 3), "\n(p = ", round(cor_dat$p.value, 3), ")"))
  }