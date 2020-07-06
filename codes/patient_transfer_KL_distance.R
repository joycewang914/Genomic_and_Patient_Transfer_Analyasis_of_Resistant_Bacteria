# Calculate Kullback-Leibler divergence of patient transfer pattern between nursing facilities
# e.g. if two nursing homes admit their patients from the same set of acute-care hospitals at equal frequencies, then dissimilarity score = 0; if completely different hospitals then dissimilarity score = 1

library(entropy)

all_pt_mat = readRDS("patient_transfer_hospital_NF.RDS")

# Get patient transfer pattern
all_pt_facil_mat = t(all_pt_mat[grep("Hospital", rownames(all_pt_mat)), grep("NF", colnames(all_pt_mat))])
all_pt_pseudo_facil_mat = sapply(colnames(all_pt_facil_mat), FUN = function(x){
  all_pt_facil_mat[,x] + (colSums(all_pt_facil_mat)/10)[x]
}) # add 1/10 of discharged from each hospital to each nursing home 

all_pt_transfer_dist = matrix(outer(rownames(all_pt_pseudo_facil_mat), rownames(all_pt_pseudo_facil_mat),
                                    Vectorize(FUN = function(x,y){
                                      KL.plugin(all_pt_pseudo_facil_mat[x,]/sum(all_pt_pseudo_facil_mat[x,]),
                                                all_pt_pseudo_facil_mat[y,]/sum(all_pt_pseudo_facil_mat[y,]))
                                    })),
                              ncol = nrow(all_pt_pseudo_facil_mat), 
                              dimnames = list(rownames(all_pt_pseudo_facil_mat), 
                                              rownames(all_pt_pseudo_facil_mat)))

# saveRDS(object = all_pt_transfer_dist, file = "patient_transfer_KL_distance.RDS")