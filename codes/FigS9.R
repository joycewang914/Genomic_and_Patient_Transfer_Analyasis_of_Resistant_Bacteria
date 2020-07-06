# Plot Fig S9 showing correlation between bacterial population structure and patient transfer or distance between nursing facilities. 

library(entropy)
source("spearman_correlation.R")

# As patient-level hospital transfer information is linked to sensitive clinical data, here we provide pre-processed 1) patient transfer data (see Fig2.R for details) and 2) distance between nursing facilities for downstream analysis.

# Patient transfer KL distance
pt_transfer_dist = readRDS("patient_transfer_KL_distance.RDS")

# Nursing facility distance data
path_fac_dist = readRDS("NF_geocode.RDS")

# Create a matrix containing facility pair names (useful for selecting which pairs to include for analysis later)
geo_mat_names =  outer(rownames(path_fac_dist), colnames(path_fac_dist), 
                       Vectorize(FUN = function(x, y){paste0(x, '_', y)}))



geo_df = data.frame('geodist' = path_fac_dist[upper.tri(path_fac_dist)],
                    'pt_transfer' = pt_transfer_dist[upper.tri(pt_transfer_dist)],
                    'pair' = geo_mat_names[upper.tri(geo_mat_names)])

# Make a correlation plot showing the correlation between geographical distance and patient transfer between nursing facilities
model = lm(geo_df$pt_transfer ~ geo_df$geodist)
plot(geo_df$geodist, geo_df$pt_transfer, col="lightblue", 
     pch = 19, 
     xlab = "Geographical distance between NFs (kilometers)", 
     ylab = "Difference in patient transfer\npattern between NF pairs",
     panel.first = abline(model, col = "pink", lwd = 3))
text(pt_transfer ~ geodist, labels=pair,data=geo_df, cex = 0.75, font=2)

cor_dat = round(spearman_rho_ci(geo_df$pt_transfer, geo_df$geodist, df = geo_df), 2)
p_val = ifelse(cor_dat[2] == 0, "< 0.001", paste0("= ", formatC(cor_dat[2], digits = 2, format = "e")))
rho_ci = paste0(cor_dat[1])
title(main = paste0("Spearman rho = ", rho_ci, "\n",
                    "p ", p_val), line = 0.5, cex.main=1)