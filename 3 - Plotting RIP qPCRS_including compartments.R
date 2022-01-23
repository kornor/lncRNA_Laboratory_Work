

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
library()

######## plotting post RIP FC


data <- read.table("SVCT_Combined_postRIP.txt", sep = "\t", header = TRUE)

## group by Primer and Antibody 
by_grp <- data %>% group_by(Primer, Antibody)

# from here you can calculate the mean of the IgG Delta (the control delta)
y <- by_grp %>% summarise(disp = mean(Delta_IgG))

# you can make the Delta Delta too, right off the bat
x <-ddply(data,c("Primer","Antibody"),summarise,DeltaDelta =Delta_RBG - mean(Delta_IgG))

#and now the fold change
x$FC <- 2^-x$DeltaDelta

## wash delta deltas
x1 <-ddply(data,c("Primer","Antibody"),summarise,DeltaDelta =Delta_wash - mean(Delta_IgG))

#and now the fold change
x1$FC <- 2^-x1$DeltaDelta

y1 <-ddply(x1,c("Primer","Antibody"),summarise,meanFC =mean(FC), sdFC = sd(FC))
y1$log <- log2(y1$meanFC)


## now create a table with the Mean and SD FC
# this is a smaller table so we gotta split it
y <-ddply(x,c("Primer","Antibody"),summarise,meanFC =mean(FC), sdFC = sd(FC))
y$log <- log2(y$meanFC)


y <- y %>% mutate(Primer = factor(Primer, levels=c("LINC00022", "LINC00511", "LINC00665", "LINC01006", "LINC01315", "ec_CEBPA", "HOTAIR", "RNA18S5", "GAPDH")))
levels(y$Primer)

#### this is quite nice
p <- ggplot(y, aes(x=Primer, y= log, fill = Primer)) + 
  scale_x_discrete(limits = rev(levels(y$Primer))) +
  facet_wrap(~Antibody) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()


p + scale_fill_brewer(palette="Paired") + theme_minimal() +
  theme(text = element_text(size = 18)) +
  labs(x = "lncRNA", y = "Log2 Fold Change", tag = "A", title = "SVCT " )

write.table(y, "SVCT_FC_Calc.txt", sep = "\t")
############ SUM 159


data <- read.table("SUM_Combined_postRIP.txt", sep = "\t", header = TRUE)

## group by Primer and Antibody 
by_grp <- data %>% group_by(Primer, Antibody)

# from here you can calculate the mean of the IgG Detal (the control delta)
y <- by_grp %>% summarise(disp = mean(Delta_IgG))

# you can make the Delta Delta too, right off the bat
x <-ddply(data,c("Primer","Antibody"),summarise,DeltaDelta =Delta_RBG - mean(Delta_IgG))

#and now the fold change
x$FC <- 2^-x$DeltaDelta


## now create a table with the Mean and SD FC
# this is a smaller table so we gotta split it
y <-ddply(x,c("Primer","Antibody"),summarise,meanFC =mean(FC), sdFC = sd(FC))
y$log <- log2(y$meanFC)


y <- y %>% mutate(Primer = factor(Primer, levels=c("LINC00022", "LINC00511", "LINC00665", "LINC01006", "LINC01315", "ec_CEBPA", "HOTAIR", "RNA18S5", "GAPDH")))
levels(y$Primer)

#### this is quite nice
p <- ggplot(y, aes(x=Primer, y= log, fill = Primer)) + 
  scale_x_discrete(limits = rev(levels(y$Primer))) +
  facet_wrap(~Antibody) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()


p + scale_fill_brewer(palette="Paired") + theme_minimal() +
  theme(text = element_text(size = 18)) +
  labs(x = "lncRNA", y = "Log2 Fold Change", tag = "B", title = "SUM159" )


write.table(y, "SUM_FC_Calc.txt", sep = "\t")




####
#### Hmmm... the fold change looks different for the GAPDH etc because they were bigger before ---- 
# what about just plotting the Ct??
y1 <- y1 %>% mutate(Primer = factor(Primer, levels=c("LINC00022", "LINC00511", "LINC00665", "LINC01006", "LINC01315", "ec_CEBPA", "HOTAIR", "RNA18S5", "GAPDH")))

p <-  ggplot(y1, aes(x= Primer, y= meanFC, fill=Primer)) + 
  facet_wrap(~Antibody) +
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_errorbar(aes(ymin=meanFC-sdFC, ymax=meanFC+sdFC), width=.2,position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") + theme_minimal() +
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), legend.position = "none") +
  labs(x = "lncRNA", y = "Fold Change", tag = "A", title = "SVCT" )


p <-  ggplot(x1, aes(x= Primer, y= FC, fill=Primer)) + 
  facet_wrap(~Antibody) +
  geom_dotplot(binaxis='y', stackdir='center')

p

p <- ggplot(y1, aes(x=Primer, y= log, fill = Primer)) + 
  #scale_x_discrete(limits = rev(levels(y1$Primer))) +
  facet_wrap(~Antibody) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip()
p



####

foo <- data %>% group_by(Primer, Antibody)
ya <-ddply(data,c("Primer","Antibody"),summarise,Ct = mean(Mean.Ct.wash.sample), SD = sd(Mean.Ct.wash.sample))

ya <- ya %>% mutate(Primer = factor(Primer, levels=c("LINC00022", "LINC00511", "LINC00665", "LINC01006", "LINC01315", "ec_CEBPA", "HOTAIR", "RNA18S5", "GAPDH")))

p <-  ggplot(ya, aes(x= Primer, y= Ct, fill=Primer)) + 
  facet_wrap(~Antibody) +
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_errorbar(aes(ymin=Ct-SD, ymax=Ct+SD), width=.2,position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") +
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), legend.position = "none") +
  labs(x = "lncRNA", y = "Ct", tag = "A", title = "Wash Sample qPCR - SVCT" )


