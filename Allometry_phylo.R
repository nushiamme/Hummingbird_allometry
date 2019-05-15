## Hummingbird FMR allometry
## Paper authors: Anusha Shankar, Catherine H Graham, Donald R Powers
## Code by: Anusha Shankar, github/nushiamme; 
# contact: anusha<dot>shankar<at>stonybrook<dot>edu
## MCMCglmm models, accounting for both the phylogenetic structure and 
# the repeated-measures per species

### Contents
## Setup, read files in, format data
## Figures

library(MCMCglmm)
library(nlme)
library(ape)
library(geiger) # for treedata() function
library(caper)
#library(coda) # only for autocorr function
library(phytools)
#library(tibble) # To add columns to datasets with control
library(RColorBrewer)
library(ggplot2)

#### Setup ####
setwd("C:\\Users\\nushi\\Dropbox\\DLW_paper\\Data")
fmr_data <- read.csv("DLW_TableS1.csv") #Compiled data from this paper and literature. Each row is an individual

## Read in McGuire et al. 2014 hummingbird phylogeny
tree_dlw<-read.tree("hum294.tre")
#tre_ou_edited <- read.tree("OU_hummer_tree_FMR_edit.txt")

## General plotting functions
## Generic theme
my_theme <- theme_classic(base_size = 30) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) 

## To add linear regression equation to plot
lm_eqn <- function(y, x){
  m <- lm(y ~ x);
  eq <- substitute(italic(y) == 
                     a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

## Aggregating dataset by Species and site, to get species means of mass and daily energy expenditure (DEE)
dlw_mean <- data.frame()
dlw_mean <- aggregate(fmr_data$kJ_day, by=list(fmr_data$Species, fmr_data$Big_site), FUN="mean", na.omit=T)
dlw_mass <- aggregate(fmr_data$Mass_g, by=list(fmr_data$Species, fmr_data$Big_site), FUN="mean", na.omit=T)
dlw_mean <- merge(dlw_mean, dlw_mass, by = c("Group.1", "Group.2"))
names(dlw_mean) <- c("Species", "Region", "kJ_day", "Mass_g")

## Trimming tree to DLW dataset
## Manually replacing because it's a manageable number
tree_dlw$tip.label[1]<-"FLME"
tree_dlw$tip.label[15]<-"PHYA"
tree_dlw$tip.label[83]<-"URBE"
tree_dlw$tip.label[92]<-"HEIM"
tree_dlw$tip.label[93]<-"HERU"
tree_dlw$tip.label[95]<-"HEJA"
tree_dlw$tip.label[128]<-"AGCO"
#tree_dlw$tip.label[154]<-"PAGI" ## Doing this one separately later
tree_dlw$tip.label[156]<-"EUFU"
tree_dlw$tip.label[163]<-"LACL"
tree_dlw$tip.label[185]<-"ARAL"
tree_dlw$tip.label[188]<-"CAAN"
tree_dlw$tip.label[219]<-"CYLA"
tree_dlw$tip.label[230]<-"CHUR"
tree_dlw$tip.label[234]<-"THCO"
tree_dlw$tip.label[235]<-"THFA"
tree_dlw$tip.label[269]<-"AMTZ"


## Tree without the Giant hummingbird
tree_no_Pgigas <- tree_dlw

tree_dlw$tip.label[154]<-"PAGI"

tips<-data.frame(levels(fmr_data$Species))
colnames(tips) <- "tips"
rownames(tips)<-tips$tips

#match tree to data, prune tree, species names should be in rownnames of "data" 
tre1<-treedata(tree_dlw, tips)$phy
#To check that the relationships between species in the trimmed tree look right
plot(tre1, cex=1.5, edge.width = 3) 

## Matching tree without P. gigas and trimming
tips2<-data.frame(levels(droplevels(fmr_data$Species[fmr_data$Species != "PAGI"])))
colnames(tips2) <- "tips"
rownames(tips2)<-tips2$tips

#match tree to data, prune tree, species names should be in rownnames of "data" 
tre1_noPgigas<-treedata(tree_no_Pgigas, tips2)$phy
#To check that the relationships between species in the trimmed tree look right
plot(tre1_noPgigas, cex=1.5, edge.width = 3) 



## GLS models with OU vs. Brownian motion
## https://www.r-phylo.org/wiki/HowTo/PGLS
## PGLS can only take one value per species, so aggregate by mean DEE first, and then run the pgls model.
fmr.agg <- aggregate(fmr_data$kJ_day,by = list(fmr_data$Species), FUN='mean', na.rm=T)
names(fmr.agg) <- c('Species', 'kJ_day')
mass.agg <- aggregate(fmr_data$Mass_g,by = list(fmr_data$Species), FUN='mean', na.rm=T)
names(mass.agg) <- c('Species', 'Mass_g')
dee.agg <- merge(fmr.agg, mass.agg,by="Species")

## Making FMR and Mass separate objects
fmr<-dee.agg$kJ_day
mass_g<-dee.agg$Mass_g
DF.fmr<-data.frame(fmr,mass_g,row.names=dee.agg$Species)
DF.fmr <-  DF.fmr[tre1$tip.label,]
DF.fmr

## Running brownian motion tree GLS model
bm.fmr<-corBrownian(phy=tre1)
bm.gls<-gls(log(fmr)~log(mass_g),correlation=bm.fmr,data=DF.fmr)
summary(bm.gls)
plot(bm.gls)

## Running GLS model with Ornstein-Uhlenbeck tree
ou.fmr<-corMartins(1,phy=tre1)
ou.gls<-gls(log(fmr)~log(mass_g),correlation=ou.fmr,data=DF.fmr)
summary(ou.gls)
plot(ou.gls$residuals)
plot(ou.gls)


#### Models ####
## Now, to run Bayesian models with repeated measures per species (i.e. multiple individuals per species), 
#we setup an inverse matrix and set up a prior
#Using a Bayesian rather than a maximum likelihood model because with an ML model we could include 
#repeated measures, OR we could include a phylogenetic structure. 
#But to get a hierarchy, with both a phylogeny and then repeated measures 
#within the phylogeny, we need turn the phylogeny into an inverse matrix
inv.phylo<-inverseA(tre1, nodes="TIPS", scale=TRUE)
inv.phylo_noPgigas <- inverseA(tre1_noPgigas, nodes="TIPS", scale=TRUE)

## Make OU tree
tre_ou <- rescale(tre1, model = "OU", alpha=48.13674) ## Alpha from running OU gls model above
plot(tre_ou, cex=1.5, edge.width = 3)
## Inverse matrix of the OU tree - doesn't work, edge lengths are zero => star phylogeny.
#inv.phylo_ou <-inverseA(tre_ou_edited,nodes="TIPS",scale=T)

#set up a prior for a phylogenetic mixed model
#Setting priors to be very uninformative
prior<-list(G=list(G1=list(V=0.02,nu=0.02)),R=list(V=0.02,nu=0.02)) 
#run the hierarchical phyogenetic model, the name of the species 
#(repeated across rows of observations) 

#### Models ####

## kJ_dayg ~ Mass_g + Big_site + Site ; random = Species 
## Making site a dummy categorical variable
fmr_data$Site.f <- as.factor(fmr_data$Site)
fmr_data$Big_site.f <- as.factor(fmr_data$Big_site)
levels(fmr_data$Big_site.f)[match("AZ",levels(fmr_data$Big_site.f))] <- 0
levels(fmr_data$Big_site.f)[match("CH",levels(fmr_data$Big_site.f))] <- 1
levels(fmr_data$Big_site.f)[match("CR",levels(fmr_data$Big_site.f))] <- 2
levels(fmr_data$Big_site.f)[match("EC",levels(fmr_data$Big_site.f))] <- 3

## Make a variable to compare temperate and tropical measurements 
fmr_data$Temptrop <- as.factor(fmr_data$Big_site)
levels(fmr_data$Temptrop)[match("AZ",levels(fmr_data$Temptrop))] <- 0
levels(fmr_data$Temptrop)[match("CH",levels(fmr_data$Temptrop))] <- 1
levels(fmr_data$Temptrop)[match("CR",levels(fmr_data$Temptrop))] <- 1
levels(fmr_data$Temptrop)[match("EC",levels(fmr_data$Temptrop))] <- 1

## Non-log model, just out of interest. Do not use!
DEE_full_raw <-MCMCglmm(kJ_day~Mass_g+Temptrop, 
             random=~Species, ginverse = list(Species=inv.phylo$Ainv), 
             prior=prior, data=fmr_data, verbose=FALSE, nitt = 5000000, thin = 1000)
summary(DEE_full_raw)
plot(DEE_full_raw) 

## Log-log full model
DEE_log <-MCMCglmm(log(kJ_day)~log(Mass_g)+Temptrop, 
                    random=~Species, ginverse = list(Species=inv.phylo$Ainv), 
                    prior=prior, data=fmr_data, verbose=FALSE, nitt = 5000000, thin = 1000)
summary(DEE_log)
plot(DEE_log) 

## Allowing two slopes and two intercepts; earlier it was + Temptrop
DEE_full_noTree <-MCMCglmm(log(kJ_day)~log(Mass_g)*Temptrop, 
                   random=~Species, 
                   prior=prior, data=fmr_data, verbose=FALSE, nitt = 5000000, thin = 1000)
summary(DEE_full_noTree)
plot(DEE_full_noTree) 

## DEE vs Mass, log-log, with tree
DEE_log_mass <-MCMCglmm(log(kJ_day)~log(Mass_g), 
                   random=~Species, ginverse = list(Species=inv.phylo$Ainv), 
                   prior=prior, data=fmr_data, verbose=FALSE, nitt = 5000000, thin = 1000)
summary(DEE_log_mass)
plot(DEE_log_mass)

## DEE-mass, log-log, no tree, and using only species means
DEE_log_mass_means <-MCMCglmm(log(kJ_day)~log(Mass_g), 
                        random=~Species, ginverse = list(Species=inv.phylo$Ainv), 
                        prior=prior, data=dlw_mean, verbose=FALSE, nitt = 5000000, thin = 1000)
summary(DEE_log_mass_means)
plot(DEE_log_mass_means)

## Now run the MCMCglmm model without the tree, because the OU tree yields a star phylogeny
DEE_log_mass_noTree <-MCMCglmm(log(kJ_day)~log(Mass_g), 
                               random=~Species, 
                               prior=prior, data=fmr_data, verbose=FALSE, nitt = 5000000, thin = 1000)
summary(DEE_log_mass_noTree)
plot(DEE_log_mass_noTree)
confint(DEE_log_mass_noTree)


## Derive R2 from MCMCglmm results
## Source: https://www.int-res.com/articles/suppl/m561p001_supp2.pdf
R2 <- function(mod){
  fixed_eff <- colMeans(mod$Sol)
  fixed_var_comp <- var(as.vector(fixed_eff %*% t(mod$X)))
  all_randoms <- colMeans(mod$VCV)
  residual <- all_randoms[["units"]]
  random_var_comp <- sum(all_randoms) - residual
  R2 <- (fixed_var_comp + random_var_comp)/(sum(all_randoms) + fixed_var_comp)
  round(R2,3)
}

R2(DEE_log_mass_noTree)
R2(DEE_log_mass_means)

DEE_log_mass_noTree_noPgigas <-MCMCglmm(log(kJ_day)~log(Mass_g), 
                               random=~Species, 
                               prior=prior, data=fmr_data[fmr_data$Species != "PAGI",], verbose=FALSE, nitt = 5000000, thin = 1000)
summary(DEE_log_mass_noTree_noPgigas)
plot(DEE_log_mass_noTree_noPgigas)


## Plot temp and tropical individuals
colourCount <- length(unique(dlw_mean$Species))
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))


## Good graph of species means and individual points, with regression line through them. 
ggplot(NULL, aes(log(Mass_g), log(kJ_day))) + 
  geom_point(data=dlw_mean, aes(col=Species), size=6, shape=19) + 
  geom_smooth(data=fmr_data, method=lm, alpha=0.3) + 
  geom_point(data=fmr_data, aes(col=Species), shape = 19, size=4, alpha=0.5) + 
  my_theme + xlab("Log(Mass (g))") +
  scale_colour_manual(values = getPalette(colourCount)) +
  ylab("Log(kJ per day)")  + 
  theme(legend.key.height=unit(2,"line"))

## Good graph of just individual points, with regression line through them. 
ggplot(NULL, aes(log(Mass_g), log(kJ_day))) + 
  geom_smooth(data=fmr_data, method=lm, alpha=0.7) + 
  geom_point(data=fmr_data, aes(col=Species), size=4, alpha=0.7) + my_theme + xlab("Log(Mass (g))") +
  scale_colour_manual(values = getPalette(colourCount)) +
  ylab("Log(kJ per day)")  + 
  theme(legend.key.height=unit(2,"line"))


## t-test testing temperate vs. tropical sites:
## Raw DEE significant, but includes PAGI (t(60) = 5.03, p-value = 4.65e-06)
an.fmr <- aov(log(kJ_day)~Temptrop+log(Mass_g), data=fmr_data)
summary(an.fmr)  
plot(an.fmr,1) # residuals vs. fitted. Shows that residuals don't have clear non-linear pattern
plot(an.fmr,2) # Q-Q plot. Shows that residuals are normally fitted for the most part. Few outliers to watch for in next few plots
plot(an.fmr,3) # Scale-location plots. Line is horizontal, no big slope. Means that residuals are spread evenly among predictors. Shows homoscedasticity
plot(an.fmr,5) # Reisduals vs. leverage. All points are well enough away from the dashed red line; shows that no single point is overly influential.
aov_residuals <- residuals(object = an.fmr)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals)
