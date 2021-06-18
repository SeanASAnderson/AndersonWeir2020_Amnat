#########################################################
###   This file contains code for running the model   ###
###   performance tests in Anderson and Weir (2020)   ###
#########################################################

# Code here requires the package 'diverge' version 1.0.4 or higher
require(diverge)

# The following tests were run in Anderson & Weir (2020):
# 1. Test of accuracy/precision in parameter estimation as a function of sample size
# A) DA_null (alpha, sig2, psi)
# B) DA_linear (psi_slope, psi_int)
# C) DA_cat (psi1, psi2)
#
# 2. Test of type I error rate
# A) DA_null selected over correct BM_null as function of sample size
# B) DA_null selected over correct OU_null as function of sample size
# C) DA_linear OR DA_cat selected over correct DA_null as function of sample size
#
# 3. Test of type II error rate
# A) OU_null or BM_null selected over correct DA_NULL
# i) as function of sample size
# ii) as function of magnitude of psi
# iii) as function of alpha/sig2 ratio
# B) DA_null selected over correct DA_cat
# i) as function of sample size
# ii) as function of evenness of distribution of pairs among the two categories
# iii) as function of the magnitude of the difference in psi between the two categories
# C) DA_null selected over correct DA_linear as function of gradient severity
#
# 4. Tests for the effect of unacknowledged measurement error
# A) DA_NULL selected over correct BM_null as function of hidden measurement error
# B) DA_null selected over correct OU_null as function of hidden measurement error
# C) BM_null or OU_null selected over correct DA_null as function of hidden measurement error

# WARNING: Re-estimating models and parameters for a large number of datasets is computationally intensive
# I STRONGLY RECOMMEND USING PARALLEL COMPUTING RESOURCES WHERE AVAILABLE when re-estimating large numbers of replicate datsets
# The actual analysis in Anderson & Weir (2020) was split into chunks and run on a high performance computer cluster
# The computer-intensive commands below are written for a parallel application on a single multi-core machine
# These commands can be slightly modified by the user for serial applications or for multi-node parallel runs
# I've set Nsets, the no. of datasets to simulate to 10 in these analyses so code can be tested in reasonable time
# This value was Nsets=1000 for all analyses in the paper

# ==============================================
# ==============================================

# Preliminary setup

# Define age categories representative of sister-species datasets to be used in all analyses
# Spread these ages evenly among 5 different gradient values for tests of DA-linear
# Here we assume a continuous gradient of 0-60, which corresponds to a latitudinal gradient of 0-60 degrees absolute latitude
# Define category for tests of DA-cat vector such that each age has the same number of cateogry=0 and category=1
# note that these are category codes for a discrete or categorical variable, such as allopatric versus sympatric
age_cats=c(0.5, 1, 1.5, 2, 4, 8)
grad_cats=c(0, 15, 30, 45, 60)
cat_cats = c(0,1)

# Create ages and gradient values for datasets of 30, 60, 90, 120, 150, and 300 species pairs
# 'ages' will be a list containing vectors of sister pair ages of length equal to the number of sisters
# 'grad' will be a list containing vectors of gradient positions of length equal to the number of sisters
# 'cats' will be a list containing the categorical value (encoded as either a 1 or a 0) for each pair in the dataset
ages=list()
grad=list()
cats=list()
for(i in 1:6) {
	if(i==6){
	    ages[[i]]=rep(age_cats, 50)
	    grad[[i]]=c(rep(grad_cats[1], 60), rep(grad_cats[2], 60), rep(grad_cats[3], 60), rep(grad_cats[4], 60), rep(grad_cats[5], 60))
	    cats[[i]]=c(rep(cat_cats[1],150),rep(cat_cats[2],150))
	} else {
		ages[[i]]=rep(age_cats,i*5)
		grad[[i]]=c(rep(grad_cats[1], i*6), rep(grad_cats[2],i*6), rep(grad_cats[3],i*6), rep(grad_cats[4],i*6), rep(grad_cats[5],i*6))
	    cats[[i]]=c(rep(cat_cats[1],i*15), rep(cat_cats[2], i*15))
	}
}

# Define 'domain', which is the minimum and maximum value of the gradient (i.e. the continuous variable over which parameters are hypothesized to vary)
# Here we set domain to be 0-60, which would be appropriate, for instance, when using latitude as a continuous variable
domain=c(0,60)

# 'Nsets' is the number of datasets to simulate in each test
Nsets=10

## UNLESS OTHERWISE SET, THESE VALUES ARE USED IN ALL TESTS BELOW

# =============================================
# =============================================

# TEST 1: Accuracty/precision in parameter estimation

# The code here assesses the effect of dataset size on the accuracy of parameter re-estimation
# First, replicate datasets of trait divergence are simulated under a set of parameters
# Next, we estimate the maximum likelihood parameters for each dataset
# The distribution of parameters estimated should be centred around the value under which data were simulated

# =============================================
# =============================================

## TEST 1. A)

# set model parameters
alpha=0.8
sig2=0.2
psi=0.6
params= c(alpha, sig2, psi)

# simulate replicate datasets of EACH SIZE for DA_null
datasets_DAnull = lapply(X=ages, FUN=simulate_div, model = "DA_null", pars=params, Nsets=Nsets)

# Re_estimate parameters for each replicate of each dataset size
estimates_DAnull = list()
for(i in 1:length(datasets_DAnull)) {
	estimates_DAnull[[i]] = re_estimator(datasets_DAnull[[i]], model="DA_null", ages=ages[[i]], parallel=T)
}

# Find median estimate of alpha, sig2, and psi for each dataset size
meds = do.call(rbind,lapply(estimates_DAnull, function(x) c(median(x[,"alpha"]), median(x[,"sig2"]), median(x[,"psi"]))))
colnames(meds) = c("alpha", "sig2", "psi")
meds

# Find interquartile ranges for estimates of the three parameters from each dataset
iqs = do.call(rbind,lapply(estimates_DAnull, function(x) c(quantile(x[,"alpha"], probs=c(0.25, 0.75)), quantile(x[,"sig2"], probs=c(0.25,0.75)), quantile(x[,"psi"], probs=c(0.25,0.75)))))
colnames(iqs) = c("alpha.25","alpha.75","sig2.25","sig2.75","psi.25","psi.75")
iqs

## TEST 1. B)

# set model parameters
alpha=0.8
sig2=0.2
psi_slope= 0.03
psi_int = 3.5
params= c(alpha, sig2, psi_slope, psi_int)

# simulate replicate datasets of EACH size for DA_linear
datasets_DAlin = list()
for(i in 1:length(ages)){
    datasets_DAlin[[i]] = simulate_div(model="DA_linear", ages=ages[[i]], GRAD=grad[[i]], pars=params, Nsets=Nsets)
}

# re-estimate model parameters
estimates_DAlin = list()
for(i in 1:length(datasets_DAlin)) {
	estimates_DAlin[[i]] = re_estimator(datasets_DAlin[[i]], model="DA_linear", ages=ages[[i]], GRAD=grad[[i]], domain=domain, parallel=T)
}

# Find median estimate of alpha, sig2, and psi for each dataset size
meds = do.call(rbind,lapply(estimates_DAlin, function(x) c(median(x[,"psi_slope"]), median(x[,"psi_int"]))))
colnames(meds) = c("psi_slope","psi_int")
meds

# Find interquartile ranges for estimates of the three parameters from each dataset
iqs = do.call(rbind,lapply(estimates_DAlin, function(x) c(quantile(x[,"psi_slope"], probs=c(0.25, 0.75)), quantile(x[,"psi_int"], probs=c(0.25,0.75)))))
colnames(iqs) = c("psi_sl.25","psi_sl.75","psi_int.25","psi_int.75")
iqs

## TEST 1. C)

# Define parameters, parameter vector, and N
alpha=0.8
sig2=0.2
psi1=0.2
psi2=0.8
params= c(alpha, sig2, psi1, psi2)

# simulate replicate datasets for each sized dataset under DA_cat
datasets_cat = list()
for(i in 1:length(ages)) {
	datasets_cat[[i]] = simulate_div(model="DA_cat", pars=params, ages=ages[[i]], cats=cats[[i]], Nsets=Nsets)
}

# re-estimate model parameters for each dataset of each size
estimates_cat = list()
for(i in 1:length(datasets_cat)) {
	estimates_cat[[i]] = re_estimator(datasets_cat[[i]], model="DA_cat", ages=ages[[i]], cats=cats[[i]], parallel=T)
}

# Find median estimate of alpha, sig2, and psi for each dataset size
meds = do.call(rbind,lapply(estimates_cat, function(x) c(median(x[,"psi1"]), median(x[,"psi2"]))))
colnames(meds) = c("psi1","psi2")
meds

# Find interquartile ranges for estimates of the three parameters from each dataset
iqs = do.call(rbind,lapply(estimates_cat, function(x) c(quantile(x[,"psi1"], probs=c(0.25, 0.75)), quantile(x[,"psi2"], probs=c(0.25,0.75)))))
colnames(iqs) = c("psi1.25","psi1.75","psi2.25","psi2.75")
iqs


# =============================================
# =============================================

# TEST 2: Rates of Type I error

# TEST 2. A)

# define BM params
sig2=0.3
params=sig2

# simulate the datasets (NOTE: 'datasets' is a list in which each element is a matrix where nrow == no.replicate datasets and ncol == no. species pairs)
datasets_BM = lapply(X=ages, FUN=simulate_div, model = "BM_null", pars=params, GRAD=NULL, breakpoint=NULL, Nsets=Nsets)

# find the type I error rate (proportion of times a two-process model is favored over the single process BM) for datasets of different size
type1_BM_DAnull=rep(0, length(datasets_BM))
for(i in 1:length(datasets_BM)) {
	type1_BM_DAnull[i] = model_error_rate(datasets=datasets_BM[[i]], ages=ages[[i]], sim_model="BM_null", alternatives="DA_null", type=1, parallel=T)
}
type1_BM_DAnull

# TEST 2. B)

# define OU params
sig2=0.2
alpha=0.8
params=c(alpha, sig2)

# simulate the datasets
datasets_OU = lapply(X=ages, FUN=simulate_div, model = "OU_null", pars=params, GRAD=NULL, breakpoint=NULL, Nsets=Nsets)

# find the type I error rate (proportion of times a two-process model is favored over the single process OU) for datasets of different size
type1_OU_DAnull=rep(0, length(datasets_OU))
for(i in 1:length(datasets_OU)) {
	type1_OU_DAnull[i] = model_error_rate(datasets=datasets_OU[[i]], ages=ages[[i]], sim_model="OU_null", alternatives="DA_null", type=1, parallel=T)
}
type1_OU_DAnull

# TEST 2. C) 

# define DA_null params
sig2=0.2
alpha=0.8
psi=0.8
params=c(alpha, sig2, psi)

# simulate the datasets
datasets_DA = lapply(X=ages, FUN=simulate_div, model = "DA_null", pars=params, GRAD=NULL, breakpoint=NULL, Nsets=Nsets)

# find the type I error rate (proportion of times a more complex model DA model is favored over DA_null) for datasets of different size
type1_DAnull_others = rep(0, length(datasets_DA))
for(i in 1:length(datasets_DA)) {
	type1_DAnull_others[i] = model_error_rate(datasets=datasets_DA[[i]], ages=ages[[i]], sim_model="DA_null", alternatives=c("DA_linear", "DA_cat"), type=1, GRAD=grad[[i]], cats=cats[[i]], domain=domain, parallel=T)
}
type1_DAnull_others

# =============================================
# =============================================

## TEST 3: RATES OF TYPE II ERROR

## TEST 3 A)

# i) type II error rate with DAnull as function of sample size
# null models are BM_null and OU_null
# define DA_null params
sig2=0.2
alpha=0.8
psi=0.8
params=c(alpha, sig2, psi)

# simulate the datasets
datasets_DA = lapply(X=ages, FUN=simulate_div, model = "DA_null", pars=params, GRAD=NULL, breakpoint=NULL, Nsets=Nsets)

# find the type II error rate for datasets of different size
type2_DAnull = rep(0, length(datasets_DA))
for(i in 1:length(datasets_DA)) {
	type2_DAnull[i] = model_error_rate(datasets=datasets_DA[[i]], ages=ages[[i]], sim_model="DA_null", alternatives=c("BM_null","OU_null"), type=2, parallel=T)
}
type2_DAnull


# ii) type II error rate with DAnull as function of the magnitude of psi
# define a constant set of ages (10 age categories, each with 15 pairs, for a total of 150 pairs)
age_cat=c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8)
ages_psi=do.call(c, lapply(age_cat, function(x) rep(x,15)))

# define DA_null params
alpha = 0.8
beta = 0.2

# set psi
# psi = c(0.25, 0.5, 0.65, 0.75, 1, 3) # these are psis used in first round of testing; almost zero error rate for all so made them smaller to show at least some error
psi = c(0.1, 0.25, 0.5, 0.65, 0.75, 1)

# define param list
params = lapply(psi, function(x) c(alpha, beta, x))

# simulate datasets with the different parameter combos
data_DA_psi = lapply(params, FUN=simulate_div, model="DA_null", ages=ages_psi, Nsets=Nsets)

# find the proportion of times a simpler model is selected for each psi value
type2_DA_psi = unlist(lapply(data_DA_psi, FUN=model_error_rate, ages=ages_psi, sim_model="DA_null", alternatives=c("BM_null", "OU_null"), type=2, parallel=T), use.names=F)
type2_DA_psi


# iii) type II error rate with DAnull as function of alpha/sig2 ratio
# define a constant set of ages (10 age categories, each with 15 pairs, for a total of 150 pairs)
age_cat=c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8)
ages_ratio=do.call(c, lapply(age_cat, function(x) rep(x,15)))

# define params
alpha = 0.8
psi = 0.8

# set alpha/beta ratio levels (note, alpha is almost always quite a bit (often an order of magnitude) larger than beta)
ratio = c(0.1, 0.5, 1, 2, 5, 10)

# calculate betas
sig2 = alpha/ratio

# define param list
params = lapply(sig2, function(x) c(alpha, x, psi))

# simulate datasets with the different parameter combos
data_DA_ratio = lapply(params, FUN=simulate_div, model="DA_null", ages=ages_ratio, Nsets=Nsets)

# find the proportion of times a simpler model is selected for each a/b ratio
type2_DA_ratio = unlist(lapply(data_DA_ratio, FUN=model_error_rate, ages=ages_ratio, sim_model="DA_null", alternatives=c("BM_null", "OU_null"), type=2, parallel=T), use.names=F)
type2_DA_ratio

## TEST 3 B)

# i) type II error rate of DA_cat as a function of sample size
# set params
sig2=0.2
alpha=0.8
psi1=0.2
psi2=0.8
params=c(alpha, sig2, psi1, psi2)

# simulate the datasets
datasets_DAcat = list()
for(i in 1:length(ages)) {
  datasets_DAcat[[i]] = simulate_div(model="DA_cat", ages = ages[[i]], cats=cats[[i]], pars=params, Nsets=Nsets)
}

# find the type II error rate for datasets of different size
type2_DAcat_null=rep(0, length(datasets_DAcat))
for(i in 1:length(datasets_DAcat)) {
	type2_DAcat_null[i] = model_error_rate(datasets=datasets_DAcat[[i]], ages=ages[[i]], sim_model="DA_cat", alternatives="DA_null", type=2, cats=cats[[i]], parallel=T)
}
type2_DAcat_null


# ii) type II error rate with DA_cat as function of proportion of pairs belonging to category 1
# Note that tis is simply a measure of the evenness of the distribution of pairs among the two categories
# set params 
alpha = 0.8
sig2 = 0.2
psi1=0.2
psi2=0.8

# define ages (the same 6 age categories as above, each with 25 pairs, for a total of 150 pairs)
age_cat=c(0.5, 1, 1.5, 2, 4, 8)
ages_prop=rep(age_cat,25)

# define proportions to test (these are proportion of pairs that belong to category 1)
props = c(.04, .12, .2, .32, .4, .48)

# create a list of category vectors corresponding to those proportions
cat_cat=c(0, 1)
cats_prop=list()
for(i in 1:length(props)) {
	    x=150*props[i]
		cats_prop[[i]]=c(rep(cat_cat[1],x), rep(cat_cat[2], (150-x)))
}

# simulate datasets with the different category distributions
datasets_prop = list()
for(i in 1:length(cats_prop)) {
	datasets_prop[[i]] = simulate_div(model="DA_cat", pars=params, ages=ages_prop, cats=cats_prop[[i]], Nsets=Nsets)
}

# find the proportion of times DA_null is selected for each dataset
type2_DAcat_null_prop=rep(0, length(datasets_prop))
for(i in 1:length(datasets_prop)) {
	type2_DAcat_null_prop[i] = model_error_rate(datasets=datasets_prop[[i]], ages=ages_prop, sim_model="DA_cat", alternatives="DA_null", type=2, cats=cats_prop[[i]], parallel=T)
}
type2_DAcat_null_prop


# iii) type II error rate of DA_cat as function of the magnitude of difference between psi values in the two categories
# define ages (the same 6 age categories as above, each with 25 pairs, for a total of 150 pairs)
age_cat=c(0.5, 1, 1.5, 2, 4, 8)
ages_mag=rep(age_cat,25)

# define a category vector that evenly distributes pairs among the two categories
cat_cat=c(0, 1)
cats_mag = c(rep(cat_cat[1],75), rep(cat_cat[2],75))

# define params with range of psi2 values
alpha=0.8
sig2=0.2
psi1=0.2
psi2=c(0.25, 0.3, 0.4, 0.6, 1.0, 1.5)

# make param list with different disparities between psi1 and psi2
params= list()
for(i in 1:length(psi2)) params[[i]]=c(alpha, sig2, psi1, psi2[i])

# simulate replicate datasets for each sized dataset under a two-peak OU model with a wait time to an epoch shift
datasets_mag = list()
for(i in 1:length(params)) {
	datasets_mag[[i]] = simulate_div(model="DA_cat", pars=params[[i]], ages=ages_mag, cats=cats_mag, Nsets=Nsets)
}

# find the proportion of times DA_null is selected for each dataset
type2_DAcat_null_mag = sapply(X=datasets_mag, FUN=model_error_rate, ages=ages_mag, cats=cats_mag, sim_model="DA_cat", alternatives="DA_null", type=2, parallel=T)
type2_DAcat_null_mag

## TEST 3 C)

# type II error rate with DA_linear as a function of gradient severity
# define ages (the same 6 age categories as above, each with 25 pairs, for a total of 150 pairs)
age_cat=c(0.5, 1, 1.5, 2, 4, 8)
ages_sev=rep(age_cat,25)

# define a gradient vector that evenly distributes the pairs among 5 gradient positions
grad_sev=c(rep(grad_cats[1], 30), rep(grad_cats[2],30), rep(grad_cats[3],30), rep(grad_cats[4],30), rep(grad_cats[5],30))

# DA_linear params
alpha = 0.8
sig2 = 0.2

# set maximum value of gradient in the domain
dom=60

# set minimum value for psi to take over domain
minval = 0.01

# set vector slope values in order of severity for which to test error rate
psi_slope = c(-0.001, -0.005, -0.01, -0.05, -0.1, -0.5)

# create set of intercepts that, when combined with each slope, keep the psi value above the minimum value
psi_int = minval+(-dom*psi_slope)

# define param set (must be in the order alpha, sig2, psi_slope, psi_int) 
# this will be a matrix in which each column is a parameter and each row is a different set of parameters to be tried
params = cbind(rep(alpha,length(psi_slope)), rep(sig2, length(psi_slope)), psi_slope, psi_int)
# the function works faster when this is converted to a list
params = as.list(as.data.frame(t(params)))

# simulate N datasets for each param set
data_DAlin_severity = lapply(params, FUN=simulate_div, model="DA_linear", ages=ages_sev, GRAD=grad_sev, Nsets=Nsets)

# find the proportion of times the no-gradient model is selected for each level of gradient severity
type2_DAlin_severity = sapply(data_DAlin_severity, FUN=model_error_rate, ages=ages_sev, GRAD=grad_sev, domain=c(0,60), sim_model="DA_linear", alternatives="DA_null", type=2, parallel=T)
type2_DAlin_severity

# =============================================
# =============================================

## TEST 4: EFFECT OF UNACKNOWLEDGED MEASUREMENT ERROR ON MODEL SELECTION

# define ages (the same 6 age categories as above, each with 25 pairs, for a total of 150 pairs)
age_cat=c(0.5, 1, 1.5, 2, 4, 8)
ages_test=rep(age_cat,25)

# define different values of measurement error 
# these are proportions of the expected variance to which measurement error is equivalent
me = c(0.1, 0.2, 0.5, 0.65, 0.8, 1)

# define model(s) parameters
alpha = 0.8
sig2 = 0.3
psi = 0.8

# TEST 4A) Type I error when simulating under BM_null
# 1 BM-DA type 1 error rate
type1_BM_me = sapply(me, function(x) {
    sets = simulate_div(model="BM_null", pars=sig2, ages=ages_test, me_prop=x, Nsets=Nsets)
    model_error_rate(datasets = sets, ages=ages_test, sim_model="BM_null", alternatives=c("DA_null"), type = 1, parallel = TRUE)
  })
names(type1_BM_me) = paste("prop",me,sep="_")
type1_BM_me

# TEST 4B) Type I error when simulating under OU_null
# 3 OU-DA type 1 error rate
type1_OU_me = sapply(me, function(x) {
    sets = simulate_div(model="OU_null", pars=c(alpha, sig2), ages=ages_test, me_prop=x, Nsets=Nsets)
    model_error_rate(datasets = sets, ages=ages_test, sim_model="OU_null", alternatives=c("DA_null"), type = 1, parallel = TRUE)
  })
names(type1_OU_me) = paste("prop",me,sep="_")
type1_OU_me

# Test 4C) Type II error when simulating under DA_null (when either BM or OU can be the null alternative)
type2_DA_me = sapply(me, function(x) {
    sets = simulate_div(model="DA_null", pars=c(alpha, sig2, psi), ages=ages_test, me_prop=x, Nsets=Nsets)
    model_error_rate(datasets = sets, ages=ages_test, sim_model="DA_null", alternatives=c("BM_null","OU_null"), type = 2, parallel = TRUE)
  })
names(type2_DA_me) = paste("prop",me,sep="_")
type2_DA_me


### SCRIPT END ###
