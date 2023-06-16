# Load libraries
library(matrixStats)

# Read in matrix text file
mydata <- read.delim("HK2_pause_matrix.txt", header = TRUE, row.names = 1)

# Subset dataframe and add 0.001 to all values
reads <- mydata[8:208] + 0.001

# Calculate scaling factor for each 3'end and normalize data matrix
scalefactor <- 1 / reads$V101
mynorm <- reads[1:201] * scalefactor

# Get mean value for upstream and downstream edge and precision score
upstream <- rowMeans(mynorm[1:100])
downstream <- rowMeans(mynorm[102:201])
flankmean = (upstream + downstream) / 2 
mynorm$score <- 1 / flankmean

# Write out pause precision stats
precision <- cbind(mydata[1:7], mynorm[202])
write.table(precision, "HK2_precision_scores.txt", sep = "\t", row.names = TRUE)