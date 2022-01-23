

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


data <- read.table("SVCT_Control_postRIP.txt", sep = "\t", header = TRUE)

## group by Primer and Antibody 
by_grp <- data %>% group_by(Primer, Antibody)

# from here you can calculate the mean of the IgG Detal (the control delta)
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


y <- y %>% mutate(Primer = factor(Primer, levels=c("ec_CEBPA", "RNA18S5", "GAPDH")))
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

write.table(y, "SVCT_FC_Ctl_Calc.txt", sep = "\t")
############ SUM 159


data <- read.table("SUM_Control_postRIP.txt", sep = "\t", header = TRUE)

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


y <- y %>% mutate(Primer = factor(Primer, levels=c("ec_CEBPA", "RNA18S5", "GAPDH")))
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


write.table(y, "SUM_FC_Ctl_Calc.txt", sep = "\t")
