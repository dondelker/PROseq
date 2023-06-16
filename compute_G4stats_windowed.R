# Computes GC skew and percent G
# Load libraries
library(stringr)
library(matrixStats)

# Read in sequence data file
dnaseqs <- read.delim("Q4_precision_1kb.txt", header = TRUE)

# Load G percentage, GC skew and DNA windowing functions
# Computes G percentage of strings of DNA
gc.percentage = function(s)
  str_count(s, "[G]") / str_count(s, "[ACGT]")

# Computes GC skew of strings of DNA.
gc.skew = function(s) {
  C = str_count(s, "[C]")
  G = str_count(s, "[G]")
  (G - C) / (G + C)
}

# Computes G4 stat function over 50 bp windows
gc.windowed.per.gene = function(s, f, window.size=50) {
  
  # we assume these are all the same length
  n = nchar(s[1])
  
  r = NULL
  
  for(i in 1:(n + 1 - window.size)) {
    if (i %% 100 == 0)
      cat(i, "")
    y = f(substr(s, i, i + window.size - 1))
    r = rbind(r, y)
  }
  
  as.matrix(r)
}

# Calculate %G and GC skew by windowing
percentG_matrix <- gc.windowed.per.gene(dnaseqs$sequence, f = gc.percentage)
G4skew_matrix <- gc.windowed.per.gene(dnaseqs$sequence, f = gc.skew)

# Replace NaN and NA values with 0
G4skew_matrix[is.nan(G4skew_matrix)]<-0
percentG_matrix[is.nan(percentG_matrix)]<-0
G4skew_matrix[is.na(G4skew_matrix)]<-0
percentG_matrix[is.na(percentG_matrix)]<-0

# Calculate means and stdev for all coordinate positions
G4_pctG <- as.data.frame(rowMeans(percentG_matrix))
G4_G4skew <- as.data.frame(rowMeans(G4skew_matrix))
G4_pctG$std <- rowSds(percentG_matrix)
G4_G4skew$std <- rowSds(G4skew_matrix)

# Print to file 
write.table(G4_pctG,"Q4_meanPCT_G.txt", sep = "\t", row.names = FALSE)
write.table(G4_G4skew,"Q4_mean_GCskew.txt", sep = "\t", row.names = FALSE)
