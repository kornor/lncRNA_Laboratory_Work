
#set for the Prism folder
setwd("~/Bioinformatics Work/Meth & RNA/Prism");
# Load the package
library(WGCNA);
library(flashClust)
library(dplyr)
library(plyr)
library(tidyr)
library(pzfx)
library(reshape2)
library(ggplot2)

#read in the data table
dat <- read.table("Cells_vs_Lncs.txt", sep = "\t", header = TRUE)
#melt, retaining the cell line info, one coding gene and all the lnc info
mdat <- melt(dat, id.vars = c( "Cell.Line","DNMT1_MEAN"), measure.vars = 16:35)


## there is SD data in that set, let's get rid of that using "select" from dplyr
df_sub <- select(dat, -contains("_SD"))
mdat <- melt(df_sub, id.vars = c( "Cell.Line","DNMT1_MEAN"), measure.vars = 9:18)


p <- ggplot(mdat, aes(x=DNMT1_MEAN, y=value, colour = variable)) + 
  geom_point() + 
  theme_minimal()
p 

p + scale_y_continuous(trans='log10')
# zero values will present as "infinite"

## different dnmt

mdat3a <- melt(df_sub, id.vars = c( "Cell.Line","DNMT3A_MEAN"), measure.vars = 9:18)


p <- ggplot(mdat3a, aes(x=DNMT3A_MEAN, y=value, colour = variable)) + 
  geom_point() + 
  theme_minimal()
p 

p + scale_y_continuous(trans='log10')
# zero values will present as "infinite"

#read in the data table
dat <- read.table("Cells_Lncs_v3.txt", sep = "\t", header = TRUE)

p<- ggplot(dat, aes(x=lncRNA, y=Mean, group = Cell.line, color = lncRNA)) + 
  facet_wrap(~Line.type) +
  geom_point()+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(0.05))
p



p <- ggplot(dat, aes(x=lncRNA, y= Mean, fill=Cell.line)) + 
  facet_wrap(~Line.type) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") + theme_minimal()


p <- ggplot(dat, aes(x=Cell.line, y= Mean, fill=lncRNA)) + 
  facet_wrap(~Line.type) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") + theme_minimal()



p <- dat %>% group_by(Line.type) %>% ggplot(aes(x=Cell.line, y= Mean, fill=lncRNA)) + 
  facet_wrap(~Line.type) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") + theme_minimal()

p <- dat %>% group_by(Line.type) %>% ggplot(aes(x=Cell.line, y= Mean, fill=lncRNA)) + 
  facet_wrap(~Line.type) +
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") + theme_minimal()


#read in the data table
datA <- read.table("qPCR_BasalA.txt", sep = "\t", header = TRUE)

datB <- read.table("qPCR_BasalB.txt", sep = "\t", header = TRUE)



p <- ggplot(datA, aes(x=Cell.line, y= Mean, fill=lncRNA)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") + theme_minimal()



## B

p <- ggplot(datB, aes(x=Cell.line, y= Mean, fill=lncRNA)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") + theme_minimal() +
  theme(text = element_text(size = 22), legend.position = "none") +
  labs(x = "Cell Line", y = "Mean Transcript \n Expression", tag = "B", title = "lncRNA Transcript Expression \nBasal B Cell lines " )


## A
p <- ggplot(datA, aes(x=Cell.line, y= Mean, fill=lncRNA)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") + theme_minimal() +
  theme(text = element_text(size = 22), legend.position = "none") +
  labs(x = "Cell Line", y = "Mean Transcript \n Expression", tag = "A", title = "lncRNA Transcript Expression \nBasal A Cell lines " )



