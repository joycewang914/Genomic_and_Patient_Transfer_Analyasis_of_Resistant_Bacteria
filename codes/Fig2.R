# Plot Fig 2 showing the patient transfer network among regional nursing homes and acute care hospitals. 

library(igraph)

# As hospital transfer information is linked to sensitive clinical data, here we provide pre-processed patient transfer data for downstream analysis.

# All 6 nursing homes in study are included.
path_facils = as.character(1:6)

# Only include acute-care hospitals that have a total discharge of more than 10 patients into nursing homes
hospitals_gt10 = as.character(c(1, 2, 8, 9, 10, 16, 19, 20, 25, 26, 29))

# Make a matrix recording number of patients discharging from acute-care hospitals to nursing homes
all_pt_mat = matrix(0, 
                    nrow = length(path_facils) + length(hospitals_gt10), 
                    ncol = length(path_facils) + length(hospitals_gt10),
                    dimnames = list(c(paste0("NF ", path_facils), 
                                      paste0("Hospital ", hospitals_gt10)),
                                    c(paste0("NF ", path_facils), 
                                      paste0("Hospital ", hospitals_gt10))))
all_pt_mat[grep("Hospital", rownames(all_pt_mat)), grep("NF", colnames(all_pt_mat))] = c(3, 1, 0, 0, 0, 1, 2, 73, 21, 0, 9, 2, 4, 17, 0, 17, 0, 27, 0, 0, 0, 0, 2, 6, 2, 21, 39, 11, 10, 6, 4, 24, 3, 6, 4, 1, 2, 5, 5, 8, 8, 10, 65, 3, 0, 0, 23, 0, 16, 1, 1, 0, 0, 1, 0, 0, 4, 12, 0, 24, 0, 13, 0, 0, 0, 0)

# Save patient transfer matrix for other figures and analyses
# saveRDS(object = all_pt_mat, file = "patient_transfer_hospital_NF.RDS")

# Remove edges with fewer than 5 transfers
all_pt_mat[all_pt_mat <= 5] = 0

hospital_facility_graph = graph.adjacency(all_pt_mat, mode = "directed", diag = FALSE, weighted = TRUE)

E(hospital_facility_graph)$width = E(hospital_facility_graph)$weight * (10/max(E(hospital_facility_graph)$weight))
V(hospital_facility_graph)$color = ifelse(unlist(sapply(V(hospital_facility_graph)$name, FUN = function(x){if(length(grep("NF",x)) > 0) {1} else {0}})) == 1, "lightblue", "pink")
V(hospital_facility_graph)$shape = ifelse(unlist(sapply(V(hospital_facility_graph)$name, FUN = function(x){if(length(grep("NF",x)) > 0) {1} else {0}})) == 1, "square", "circle")        
V(hospital_facility_graph)$name = gsub(" ", "\n", V(hospital_facility_graph)$name)

# Assign position of each node using a previously-defined coordinate file generated with tkplot function in "igraph"

coord = matrix(as.numeric(c(349,496,365,121, 186,374, 536, 373, 165, 127, 548, 125, 637, 425, 84, 372, 367, 0, 93, 469, 238,244,94, 281, 477, 252, 221, 461, 500, 460, 365, 350, 348,592)), ncol = 2, byrow = T)

par(mar = c(2, 2, 2, 2))
plot.igraph(hospital_facility_graph, 
            vertex.label.font = 2, 
            vertex.label.color = "black", 
            vertex.frame.color = NA,
            edge.label = E(hospital_facility_graph)$weight, 
            edge.width = log10(E(hospital_facility_graph)$weight),
            edge.label.cex = 0.75, 
            edge.label.color = "red", 
            edge.label.font = 2, 
            edge.label.family="Helvetica", 
            vertex.label.family = "Helvetica",
            vertex.label.font = 2,
            vertex.label.cex = 0.75, 
            edge.arrow.size = 0.75, 
            layout = coord)
text(-1, 1, pos = 2, labels = expression(bold("A")), cex = 1.5, xpd = TRUE)

