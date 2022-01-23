
## Graphing test

#set for the Prism folder
setwd("~/Bioinformatics Work/Meth & RNA/Prism");
# Load the package
library(WGCNA);
library(flashClust)
library(dplyr)
library(tidyr)
library(pzfx)
library(ggplot2)


test <- read.table("Graph_test.txt", sep = "\t", header = TRUE)


p <- ggplot(test, aes(x=DNMT1, y=Value), colour = LncRNA) + 
    geom_point() + 
    theme_minimal()
p
