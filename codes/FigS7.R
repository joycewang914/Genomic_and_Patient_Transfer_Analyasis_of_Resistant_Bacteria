# Plot Fig S7 of permutation result to test for clustering of nursing facilities on each phylogeny.

org_tree = readRDS("all_vs_pruned_trees.RDS")

for (org in names(org_tree)){
  
  tree = org_tree[[org]][[2]]
  sts = subtrees(tree)
  nf_tree_facil = sapply(tree$tip.label, FUN = function(x){substr(strsplit(x, "-")[[1]][2], 1, 1)})
  
  perm_num = 1000
  
  # Empirical data

  temp_mat =  vector(mode = "list", length = perm_num + 1)
  names(temp_mat) = c("empirical", 1:perm_num)
  intra_tn_facil_mat = structure(rep(0, length(path_facils)), names = path_facils)
  
  for (f1 in path_facils){
    
    subtree_counts = sum(unlist(lapply(sts, FUN = function(y){sum(nf_tree_facil[y$tip.label] %in% f1) > 1})))
    intra_tn_facil_mat[f1] = subtree_counts
  }
  
  temp_mat[["empirical"]] = intra_tn_facil_mat
  
  # Randomized data
  set.seed(20191221)

  for (i in 1:perm_num){
    
    intra_tn_facil_mat = structure(rep(0, length(path_facils)), names = path_facils)
    
    rand_tree_facil = structure(sample(nf_tree_facil, replace = FALSE), names = names(nf_tree_facil))
    
    for (f1 in path_facils){
      
      rand_subtree_counts = sum(unlist(lapply(sts, FUN = function(y){
        sum(rand_tree_facil[y$tip.label] %in% f1) > 1 })))
      
      intra_tn_facil_mat[f1] = rand_subtree_counts
    }
    
    temp_mat[[i + 1]] = intra_tn_facil_mat
  }
  
  temp_mat_p = structure(rep(0.5, length(path_facils)), names = path_facils)
  
  for (p in path_facils){
    
    temp_mat_p[p] = max(1/perm_num, sum(unlist(lapply(temp_mat[!names(temp_mat) %in% "empirical"], 
                                                      FUN = function(x){
                                                        temp_mat[["empirical"]][p] < x[p]}))/perm_num))
    
  }
  
  for (p in path_facils){
    rand_data = unlist(lapply(temp_mat[2:(perm_num+1)], FUN = function(x){x[p]}))
    hist(rand_data, 50, xlab = "# subtrees", ylab = "Frequency", xlim= c(0, max(rand_data)+1),
         main = paste0("NF ",p, " (N = ", table(nf_tree_facil)[p], ")\nperm p = ", temp_mat_p[p]),
         cex.main = 2, cex.lab = 2, cex.axis = 2)
    
    abline(v =temp_mat$empirical[p], col = "red")
    
  }
  
}
