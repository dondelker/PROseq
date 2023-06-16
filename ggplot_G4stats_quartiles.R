#Make metaplots of %G and GC skew
#Load libraries
library(ggplot2)
library(tidyverse)

# Read in and combine G4 stats files 
quart4 <- read.delim("Q4_meanPCT_G.txt", header = TRUE)
quart3 <- read.delim("Q3_meanPCT_G.txt", header = TRUE)
quart2 <- read.delim("Q2_meanPCT_G.txt", header = TRUE)
quart1 <- read.delim("Q1_meanPCT_G.txt", header = TRUE)

myG4 <- cbind(quart1[1],quart2[1],quart3[1],quart4[1])
colnames(myG4) <- c("quartile1", "quartile2", "quartile3", "quartile4")

# Reshape the dataframe so it is 'ggplot'-able
df.gg <- myG4 %>% gather(quartile, percentG, quartile1:quartile4)
df.gg$Index <- rep(1:951, 2)

# Add a column for standard error
newstd <- merge(quart1, y=c(quart2,quart3,quart4), all = TRUE)
df.gg$se <- newstd$std / 55

# Make metaplot
ggplot(df.gg, aes(x=Index, y=percentG, Group=factor(quartile))) +
  geom_ribbon(aes(ymin=percentG-se, ymax=percentG+se),
              alpha=0.2) +
  geom_line(aes(colour=factor(quartile))) +
  scale_colour_manual(values = c("#F8766D", "orange2", "blue3", "grey10")) +
  scale_x_continuous(breaks = c(0, 238, 475, 713, 950),
                     labels = c("-500", "-250", "Pause Max", "+250", "+500")) +
  geom_vline(xintercept=475, colour="red", linetype="longdash") +
  labs(colour="quartile") + xlab("Position(BP)") +
  ggtitle("Metagene plot of percent G") +
  theme_bw()         