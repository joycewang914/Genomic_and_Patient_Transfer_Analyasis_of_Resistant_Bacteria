# Plot Fig S6 showing the source of each patient isolate by nursing facility. 

# Prune genetically similar isolates found within the same NF during follow-up visits to rule out the possiblity of within-facility transmission. 

# Reroot each phylgoenetic tree at midpoint

# Add sequence type information by each isolate.

library(ape)
library(phytools)
library(hues)


# Load RaXML trees
tree_mrsa = read.tree("RAxML_bipartitions.2019_03_19_14_07_22_803-A001-120-N-MRSA-TIP__ordered_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites_raxML")
tree_vrefm = read.tree("RAxML_bipartitions.2019_03_19_18_03_26_3399-4058-0-RVRE_final_ordered_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites_raxML")
tree_vrefc = read.tree("RAxML_bipartitions.2019_03_19_14_29_34_7605-4117-0-HVRE_ordered_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites_raxML")
tree_ecol = read.tree("RAxML_bipartitions.2019_03_19_10_08_45_3425-5015-0-R2X_ordered_genome_aln_w_alt_allele_unmapped_gubbins.filtered_polymorphic_sites_raxML")
mdro_trees = list("MRSA" = tree_mrsa, "VREfm" = tree_vrefm, "VREfc" = tree_vrefc, "CipREc" = tree_ecol)
path_trees = lapply(names(path_dna), FUN = function(p){
  isolates = rownames(path_dna[[p]])
  keep.tip(mdro_trees[[p]], isolates)
})
names(path_trees) = names(path_dna)
rooted_path_trees = lapply(path_trees, midpoint.root)
remove(path_trees)

# Load MLST data
mlst_ecoli = read.table("CipREc_MLST.txt"); 
mlst_saureus = read.table("MRSA_MLST.txt"); 
mlst_efc = read.csv("VREfc_MLST.csv")
rownames(mlst_efc) = mlst_efc$X
mlst_efm = read.table("VREfm_MLST.txt"); 
mlst_efm = mlst_efm[!rownames(mlst_efm) %in% c("702-3020-0-RVRE_S18_L001", "546-4011-0-RVRE_S58_L001"), ]
mlst_clusters = list("MRSA" = mlst_saureus, "VREfm" = mlst_efm, "VREfc" = mlst_efc, "CipREc" = mlst_ecoli)

path_tree_mlst = list()
for (m in names(path_dna)){
  
  mlst = mlst_clusters[[m]]
  path_tree = rooted_path_trees[[m]]
  
  path_mlst = mlst[sapply(rownames(mlst), FUN = function(x){substr(strsplit(x, "-")[[1]][2],1,1) %in% path_facils}), ]
  
  overlap_ind = intersect(substr(rownames(path_mlst), 1, 8), substr(path_tree$tip.label, 1, 8))
  
  path_tree_mlst[[m]] = list("mlst" = path_mlst[grep(paste(overlap_ind, collapse = "|"), rownames(path_mlst)),], 
                             "tree" = keep.tip(path_tree, grep(paste(overlap_ind, collapse = "|"), path_tree$tip.label)))
  
}

# Plot trees - by NF
nf_cols = structure(iwanthue(length(path_facils)), names = path_facils)
nf_cols[nf_cols %in% "#4C3E45"] = "orange"

org_tree = list()
for (org in names(rooted_path_trees)){print(org)
  par(mar = c(5,5,5,10))
  tree = rooted_path_trees[[org]]
  tree_visit = sapply(tree$tip.label, FUN = function(x){strsplit(x, "-")[[1]][3]})
  
  # Isolates coloured by NF
  # isolates collected at the time of NF admission are shown as solid
  # circles. 
  # Follow-up isolates that are genetically distinct from other isolates within the same
  # NF are shown as circles with a black border, otherwise indicated as a solid black         circle and pruned from subsequent analysis. 
  # Sequence type (ST) of each isolate is shown on the right. 
  first_sts = subtrees(tree)
  
  pruned_tree_tips = sapply(1:length(first_sts), FUN = function(x){print(x)
    
    temp_tree = first_sts[[x]]; 
    
    temp_tree_mat = t(sapply(temp_tree$tip.label, FUN = function(x){
      cbind(substr(strsplit(x, "-")[[1]][2], 1, 1), strsplit(x, "-")[[1]][3])}))
    
    if (length(unique(temp_tree_mat[,1])) > 1){
      keep_tree = temp_tree   
      } else if (
        sum(temp_tree_mat[,2] != 0) == 0){
        keep_tree = temp_tree;
        } else if (
          sum(temp_tree_mat[,2] == 0) > 0){
           drop_tip = temp_tree_mat[temp_tree_mat[,2] != 0,,drop = FALSE]
           keep_tree = drop.tip(temp_tree, rownames(drop_tip))
        } else if (
          length(unique(temp_tree_mat[,2][temp_tree_mat[,2] != 0])) > 0){
          
          # This chunk will not run because path_visit file is not provided on github.
          # Will keep both tips for now.
          
          # pt_visit_date = sapply(rownames(temp_tree_mat),FUN = function(z){
          #   pt = strsplit(z, "-")[[1]][2]
          #   visit = strsplit(z, "-")[[1]][3]
          #   date = path_visit[path_visit$studyID %in% pt & path_visit$visit %in% visit, "date"]
          # })
          # earlier_pt_date = unique(c(names(which.min(pt_visit_date)), names(pt_visit_date[grep("-0$", names(pt_visit_date))])))
          
          earlier_pt_date = ""
          pt_visit_date = ""
          keep_tree = drop.tip(temp_tree, setdiff(names(pt_visit_date), earlier_pt_date))
          }
    
    diff = setdiff(temp_tree$tip.label, keep_tree$tip.label)
  })    
  
  pruned_tree_tips = unique(unlist(pruned_tree_tips))
  
  # NF trees
  nf_tree_facil = sapply(tree$tip.label, FUN = function(x){substr(strsplit(x, "-")[[1]][2], 1, 1)})
  nf_tree_facil_cols = sapply(nf_tree_facil, FUN = function(x){nf_cols[x]})
  
  nf_tree_facil_cols[grep(paste(pruned_tree_tips, collapse = "|"), names(nf_tree_facil_cols))] = "black"
  
  tip_shape = rep("transparent", length(nf_tree_facil_cols))
  tip_shape[!tree_visit %in% "0"] = "black"
  
  plot(tree, no.margin = FALSE, show.tip.label = FALSE,main = org, use.edge.length = T)
  tiplabels(pch = 21, col = tip_shape, cex = ifelse(length(tree$tip.label) > 100, 0.25, 0.5),
            bg = nf_tree_facil_cols)
  legend("left", legend  = names(nf_cols),
         fill = nf_cols, title = "NF", border = FALSE, bty = "n", cex = 1, inset = -.1, xpd = TRUE)
  add.scale.bar(0, -.5, cex = 0.7, font = 2, col = "black", lwd = 2, xpd = TRUE)

  # Add MLST 
  temp_mlst = path_tree_mlst[[org]]$mlst
  if (org %in% c("MRSA","CipREc")){
    temp_mlst$new_ST = paste0(names(which.max(table(temp_mlst[,"ST"]))), 
                              " and ", names(which.max(table(temp_mlst[,"ST"]))), "-like")}else{
      temp_mlst$new_ST = temp_mlst$ST
                              }
  temp_mlst = as.data.frame(temp_mlst)
  mlst_cols = structure(rainbow(length(as.character(unique(temp_mlst[,"new_ST"])))), 
                        names = as.character(unique(temp_mlst[,"new_ST"])))
  mlst_col = sapply(tree$tip.label, FUN = function(x){
    ST = temp_mlst[grep(substr(x, 1, 8), rownames(temp_mlst)),"new_ST"]; 
    mlst_cols[as.character(ST)]})
  
  mlst_col = unlist(mlst_col)
  
  par(new = TRUE)
  par(xpd = TRUE)
  plot(x = rep(110, length(path_tree$tip.label)), y = 1:length(path_tree$tip.label),
       col = mlst_col, pch = 19, cex = 0.5, xlim = c(0, 100), xlab = "", ylab = "",
       frame.plot = FALSE, xaxt='n', yaxt='n')
  text(110, length(path_tree$tip.label) + 3, "ST")
  legend(110, length(path_tree$tip.label) + 5,  names(mlst_cols[mlst_cols %in% unique(mlst_col)]),
         fill =   mlst_cols[mlst_cols %in% unique(mlst_col)], xpd = TRUE, horiz = FALSE, title = expression(bold("")),
         bty = "n", x.intersp = 0.25, text.width = 0.05, title.col = "#606060", border = FALSE)
  
  org_tree[[org]] = list("All isolates" = tree, "Putative intra-facility transmission\nisolates removed" = drop.tip(tree, pruned_tree_tips))
}

saveRDS(object = org_tree, file = "all_vs_pruned_trees.RDS")