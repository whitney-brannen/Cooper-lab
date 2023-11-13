###-------------------------------------
###
###-------------------------------------
setwd("/Users/whitneybrannen/CooperLab/EnvironmentalAssociations/")

library(LandGenCourse)
library(vegan)   
library(lfmm)     
library(qvalue)   

# read genotype data
gen <- read.csv("ctarsalis_genotypes.csv")

# remove NAs
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp))

# read environment data
env <- read.csv("gps_annual_monthly_climate_samples_edited.csv")

# make sure they line up
env$Sample <- as.character(env$Sample)
identical(rownames(gen.imp), env[,1]) 

# subset only environmental predictors
pred <- env[,6:18]

# run pca using rda from vegan for environmental
pred.pca <- rda(pred, scale=T)
summary(pred.pca)$cont

# scree plot
screeplot(pred.pca, main = "Screeplot: Eigenvalues of C. tarsalis Predictor Variables")

# correlations between the PC axis and predictors:
round(scores(pred.pca, choices=1:10, display="species", scaling=0), digits=3)

# store synthetic PC axis predictor as pred.PC1 for use in LFMM
pred.PC1 <- scores(pred.pca, choices=1, display="sites", scaling=0)

###--------------------------------------------
### Latent Factor Mixed Models (LFMM) Univariate GEA
###--------------------------------------------

# determine K
screeplot(pred.pca, main = "Screeplot of C. tarsalis Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")

# run pca for genetic 
gen.pca <- rda(gen.imp, scale=T)
screeplot(gen.pca, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")

K <- 3

# run lfmm
ctarsalis.lfmm <- lfmm_ridge(Y=gen.imp, X=pred.PC1, K=K) ## change K as you see fit

ctarsalis.pv <- lfmm_test(Y=gen.imp, X=pred.PC1, lfmm=ctarsalis.lfmm, calibrate="gif")

names(ctarsalis.pv) 

ctarsalis.pv$gif

# check pvalues
hist(ctarsalis.pv$pvalue[,1], main="Unadjusted p-values")        
hist(ctarsalis.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

# change the GIF and readjust the p-values
zscore <- ctarsalis.pv$score[,1]   
(gif <- ctarsalis.pv$gif[1])    

# Manual adjustment of the p-values
# ???? if needed
# 

# q values
ctarsalis.qv <- qvalue(ctarsalis.pv$calibrated.pvalue)$qvalues

# how many SNPs have an FDR < 10%?
length(which(ctarsalis.qv < 0.1)) 

# identify the snps!
(ctarsalis.FDR.1 <- colnames(gen.imp)[which(ctarsalis.qv < 0.1)]) 
# these are displaying the row name, have to go back to original df to get contig


###----------------------------------------------
### Redundancy Analysis (RDA): multivariate GEA
###----------------------------------------------

# run rda
ctarsalis.rda <- rda(gen.imp ~ ., data=pred, scale=T)
ctarsalis.rda

RsquareAdj(ctarsalis.rda)

# summary of variance by pcs
summary(ctarsalis.rda)$concont

screeplot(ctarsalis.rda)

# check for variance inflation factors
vif.cca(ctarsalis.rda)

plot(ctarsalis.rda, scaling=3) 

# identify rda candidates
load.rda <- summary(ctarsalis.rda)$species[,1:3]

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

# identify which snps call in the tails of above figures
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) ## f.nd loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}

# apply to first 3 constrained axes
cand1 <- outliers(load.rda[,1], 3) ## 3.
cand2 <- outliers(load.rda[,2], 3) ## 6.
cand3 <- outliers(load.rda[,3], 3) ## 3.

ctarsalis.rda.cand <- c(names(cand1), names(cand2), names(cand3)) 

length(ctarsalis.rda.cand[duplicated(ctarsalis.rda.cand)])

ctarsalis.rda.cand <- ctarsalis.rda.cand[!duplicated(ctarsalis.rda.cand)] 

# Set up the color scheme for plotting:
bgcol  <- ifelse(colnames(gen.imp) %in% ctarsalis.rda.cand, 'gray32', '#00000000')
snpcol <- ifelse(colnames(gen.imp) %in% ctarsalis.rda.cand, 'red', '#00000000')

## a.es 1 & 2 - zooming in to just the SNPs here...
plot(ctarsalis.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), main="ctarsalis RDA, axes 1 and 2")
points(ctarsalis.rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(ctarsalis.rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(ctarsalis.rda, scaling=3, display="bp", col="#0868ac", cex=1)

## a.es 2 & 3
plot(ctarsalis.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3), main="ctarsalis RDA, axes 2 and 3")
points(ctarsalis.rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3, choices=c(2,3))
points(ctarsalis.rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3, choices=c(2,3))
text(ctarsalis.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(2,3))


# which environmental factors are most strongly correlated with first 3 RDA
intersetcor(ctarsalis.rda)[,1:3]


###---------------------------------
### Compare LFMM and RDA candidates
###---------------------------------

# check overlap
intersect(ctarsalis.FDR.1, ctarsalis.rda.cand)

# unique to lfmm
setdiff(ctarsalis.FDR.1, ctarsalis.rda.cand)   
