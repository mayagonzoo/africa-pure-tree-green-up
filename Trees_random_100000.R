library(dplyr)
library(lubridate)

############################################################
## File paths
############################################################
data_dir <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"
out_dir  <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"

paired_tr = read.csv(file.path(data_dir, "FULL_paired_tree_pts.csv"))

# Shuffle the dataset once (without replacement)
shuffled = paired_tr[sample(nrow(paired_tr)), ]

samp1 = shuffled[1:20000, ]
samp2 = shuffled[20001:40000, ]
samp3 = shuffled[40001:60000, ]
samp4 = shuffled[60001:80000, ]
samp5 = shuffled[80001:100000, ]

# Export the five samples (total of 100,000 points)
write.csv(samp1, file = file.path(out_dir, "TREE_BATCH_1.csv"),
          row.names = FALSE)
write.csv(samp2, file = file.path(out_dir, "TREE_BATCH_2.csv"),
          row.names = FALSE)
write.csv(samp3, file = file.path(out_dir, "TREE_BATCH_3.csv"),
          row.names = FALSE)
write.csv(samp4, file = file.path(out_dir, "TREE_BATCH_4.csv"),
          row.names = FALSE)
write.csv(samp5, file = file.path(out_dir, "TREE_BATCH_5.csv"),
          row.names = FALSE)