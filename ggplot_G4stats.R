#Make metaplots of %G and GC skew
#Load libraries
library(ggplot2)
library(tidyverse)

# Read in and combine G4 stats files 
quad4 <- read.delim("Quad4_mean_GCskew.txt", header = TRUE)
single4 <- read.delim("SingG4_mean_GCskew.txt", header = TRUE)
myG4 <- cbind(single4[1],quad4[1])
colnames(myG4) <- c("singleG4", "quadG4")

# Reshape the dataframe so it is 'ggplot'-able
df.gg <- myG4 %>% gather(G4_type, GC_skew, singleG4:quadG4)
df.gg$Index <- rep(1:951, 2)

# Add a column for standard error
newstd <- merge(single4, quad4, all = TRUE)
df.gg$se <- newstd$std / 55

# Make metaplot
ggplot(df.gg, aes(x=Index, y=GC_skew, Group=factor(G4_type))) +
  geom_ribbon(aes(ymin=GC_skew-se, ymax=GC_skew+se),
              alpha=0.2) +
  geom_line(aes(colour=factor(G4_type))) +
  scale_colour_manual(values = c("#F8766D", "grey10")) + ###7CAE00 ###F8766D
  scale_x_continuous(breaks = c(0, 238, 475, 713, 950),
                     labels = c("-500", "-250", "G4 center", "+250", "+500")) +
    geom_vline(xintercept=475, colour="red", linetype="longdash") +
  labs(colour="G4_Type") + xlab("Position(BP)") +
  ggtitle("Metagene plot of GC_skew") +
  theme_bw()        