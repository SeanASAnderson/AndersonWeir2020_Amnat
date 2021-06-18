##########################################################
### This script contains code for running the analyses ###
### in the Empirical Examples section of Anderson and  ###
### Weir (2020), The American Naturalist               ###
##########################################################

# see model performance script for instructions on installing the r-package 'diverge'
require(diverge)
require(ape)

# NOTES ON DATA
# 1. The ‘lawson_weir14’ dataset is a curated versions of the dataset originally published by 
# Lawson and Weir (2014). This curated version is available on the dryad repository for this paper,
# link: http://doi.org/10.5061/dryad/.0p2ngflxc
# The complete dataset is available in the supplementary material (S2 file) of the original paper 
# at https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.12346
# Citation: Weir, J.T., and A. Lawson. 2014. Latitudinal gradients in climatic‐niche
# evolution accelerate trait evolution at high latitudes. Ecology Letters 14:1427-1436.
# 2. The phylogeny of mammals is taken from Upham et al. (2019)
# Citation: Upham, S.H., Esselstyn, J.A., and W. Jetz. 2019. Inferring the mammal tree: species-
# level sets of phylogenies for questions in ecology, evolution, and conservation. PLoS
# Biology 17(12)::e3000494.
# 3. Body size in mammals is taken from Smith et al. (2003)
# Citation: Smith, A.F. et al. 2003. Body mass of late quarternary mammals. Ecology 84:3403-3403.
# The mammal tree and body size datasets are available in the dryad repository for this paper,
# link: http://doi.org/10.5061/dryad/.0p2ngflxc
# as well as at the relevant web locations for the original papers.

# =========================================================
### 1. Climate niche evolution in new world birds ##
# load data
clim = read.csv("lawson_weir14_climate.csv")
head(clim)

# extract key vectors and divide species pairs into tropical and temperate categories
# note that these categories are coded as 1 or 0 in the 'cats' vector
ages = clim[,"ages"]
grad = clim[,"abs_centroid_lat"]
cats = rep(NA, length(grad))
for(i in 1:length(grad)) {
	if(grad[i]>23.5) cats[i]=1
	if(grad[i]<23.5) cats[i]=0
}

# model testing for PC1
# note 1: when testing a continuous variable, the domain over which the variable can vary must be define
# here, testing the variable of latitude, we define a domain of 0-60 degrees (domain=c(0,60))
# note 2: the diverge package can test for an OU_linear model in which either alpha or sig2 varies continuously
# OU_linear_sig denotes the latter and was used in the original Lawson and Weir (2014) analysis
pc1 = clim[,"clim_div_pc1"]
respc1 = model_select(div = pc1, ages = ages, GRAD=grad, cats=cats, domain=c(0, max(grad)), models=c("BM_null","BM_linear", "OU_null","OU_linear_sig", "DA_null","DA_linear","DA_cat"))
respc1

# model testing for PC2
pc2 = clim[,"clim_div_pc2"]
respc2 = model_select(div = pc2, ages = ages, GRAD=grad, cats=cats, domain=c(0, max(grad)), models=c("BM_null","BM_linear", "OU_null","OU_linear_sig", "DA_null","DA_linear","DA_cat"))
respc2

# model testing for PC3
pc3 = clim[,"clim_div_pc3"]
respc3 = model_select(div = pc3, ages = ages, GRAD=grad, cats=cats, domain=c(0, max(grad)), models=c("BM_null","BM_linear", "OU_null","OU_linear_sig", "DA_null","DA_linear","DA_cat"))
respc3

# =========================================================
### 2. Body size evolution in new world birds ##
# load data
bodsize = read.csv("lawson_weir14_bodysize.csv")
head(bodsize)

# extract key vectors and divide species pairs into tropical and temperate categories
# note that these categories are coded as 1 or 0 in the 'cats' vector
ages = bodsize[,"ages"]
grad = bodsize[,"abs_centroid_lat"]
cats = rep(NA, length(grad))
for(i in 1:length(grad)) {
	if(grad[i]>23.5) cats[i]=1
	if(grad[i]<23.5) cats[i]=0
}

# model testing for body size
# note that when testing a continuous variable, the domain over which the variable can vary must be define
# Here, testing the variable of latitude, we define a domain of 0-60 degrees (domain=c(0,60))
bs = bodsize[,"ed_mass"]
resbs = model_select(div = bs, ages = ages, GRAD=grad, cats=cats, domain=c(0, 60), models=c("BM_null","BM_linear", "OU_null","OU_linear_sig", "DA_null","DA_linear","DA_cat"))
resbs

# =========================================================
# 3. Body size evolution in mammals
# NOTE: in this example, rather than providing a ready-to-go species pair dataset, 
# we demonstrate how to generate one from a published phylogeny and a published 
# set of body size data in mammals. 

## 3A. Isolating sister species pairs from a mammalian supertree ##
# load the dated phylogeny from Upham et al. (2019)
tre = read.nexus("mammal_dna_tree.nex")
#plot(tre, type="fan", cex=0.15, no.margin=T)

# extract the identity and age of each sister species pair in the tree
# note: since branch lengths of this tree are in units of ma, there's no need to designate 
# a molecular clock rate as an argument in 'extract_pairs'
sp_pairs = extract_sisters(tree=tre, sis_age=T)

# view the formatting of sp_pairs
head(sp_pairs)

# rename the tip labels to be in the simpler form "Genus_species", which matches the format for taxa ID in the body size dataset
sp_pairs$sp1 = sapply(strsplit(as.character(sp_pairs$sp1), "_"), function(x) paste(x[1],x[2],sep="_"))
sp_pairs$sp2= sapply(strsplit(as.character(sp_pairs$sp2), "_"), function(x) paste(x[1],x[2],sep="_"))
head(sp_pairs)

## 3B. Finding data for sister species pairs from a global dataset of body size in mammals ##
# load body size data
dat = read.csv("mammal_size.csv")

# grab body size data for each species pair from the global dataset
# some pairs are present more than once in the dataset, so here we average their body size
sp_pairs = cbind(sp_pairs, sp1_mass = rep(NA, nrow(sp_pairs)), sp2_mass = rep(NA, nrow(sp_pairs)))
for(i in 1:nrow(sp_pairs)) {
	if(sp_pairs$sp1[i] %in% dat$gen_sp) {
      if(length(which(dat$gen_sp == sp_pairs$sp1[i])>1)) {
		sp_pairs[i,"sp1_mass"] = mean(dat[dat$gen_sp==sp_pairs$sp1[i],"log10mass"])
		} else {
		sp_pairs[i,"sp1_mass"] = dat[dat$gen_sp==sp_pairs$sp1[i],"log10mass"]
	    }
	}
	if(sp_pairs$sp2[i] %in% dat$gen_sp) {
      if(length(which(dat$gen_sp == sp_pairs$sp2[i])>1)) {
		sp_pairs[i,"sp2_mass"] = mean(dat[dat$gen_sp==sp_pairs$sp2[i],"log10mass"])
		} else {
	    sp_pairs[i,"sp2_mass"] = dat[dat$gen_sp==sp_pairs$sp2[i],"log10mass"]
	    }
    }
}

# trim the sister pair dataset to include just pairs for which size data are available for both spp.
# note that this includes 'NA's for species that weren't in the larger dataset 
# as well as missing data in that dataset which were denoted "-999.00" by the original authors
sp_pairs_trimmed = sp_pairs[-which(is.na(sp_pairs$sp1_mass)|is.na(sp_pairs$sp2_mass)),]
sp_pairs_trimmed = sp_pairs_trimmed[-which(sp_pairs_trimmed$sp1_mass == -999.00| sp_pairs_trimmed$sp2_mass == -999.00),]
head(sp_pairs_trimmed)

## 3C. Calculate trait differentiation and run model selection ##
# calculate differentiation in mass and add to the sp_pairs_trimmed dataframe
sp_pairs_trimmed = cbind(sp_pairs_trimmed, mass_diff = abs(sp_pairs_trimmed$sp1_mass-sp_pairs_trimmed$sp2_mass))
head(sp_pairs_trimmed)

# remove two very strange outliers (the two pairs containing Dorcopsulus_vanheurni & Eliomys_quercinus)
sp_pairs_trimmed = sp_pairs_trimmed[-order(sp_pairs_trimmed$mass_diff, decreasing=T)[1:2],]

# observe the distribution of body mass differentiation
hist(sp_pairs_trimmed$mass_diff, freq=F, breaks=20)

# run model selection on body mass in mammals
res_mam = model_select(div = sp_pairs_trimmed$mass_diff, ages = sp_pairs_trimmed$pair_age, models=c("BM_null","OU_null","DA_null"))
res_mam

## SCRIPT END ##





