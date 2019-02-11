## After loading custom functions,
## User can skip down to Run code section

## Load package
library(data.table)

## Load data
pPositivePlot <- fread("pPositivePlot.csv")


## Define custom functions

## this function calculates the probability of 
## detecting eDNA in at least one sample
## over a range of sample sizes given 
## a detection probability pPos
calcDetectOne <- function(minSample = 1,
                          maxSample = 1,
                          pPos = 0.5){
    samples <- minSample:maxSample
    pDetectOne <- 1 - (1 - pPos)^samples
    return(pDetectOne)
}
calcDetectOne(maxSample = 4)

## This function takes a vector, x, and finds the nearst value 
## to "desiredProb". It is largely an internal function used later.
nearestValue <- function(desiredProb = 1/3, x = seq(0,1, by = 0.25)){
    which(abs(x-desiredProb)==min(abs(x-desiredProb)))
}
nearestValue()

## this calculates the the number of samples needed 
## to get one detection given a desired detection probability
## it requires a minimum and maximum number of samples,
## a desired detection probability
## and a vector of possible detection probabilities
## This is largely an internal function.
samplesForOne <- function(minSample = 1, maxSample = 1000,
                          desiredProb = 0.8,
                          pOne = seq(0, 1, length.out = 1000)){
    return((minSample:maxSample)[nearestValue(desiredProb = desiredProb,
                                              x = pOne)])
}
samplesForOne()

## This is a wrapper funciton that calcluates the sample size 
## needed to detect eDNA in at least one sample given 
## a desired detection probability, a min and max number of possible samples
## and a probability of positive detection 
sampleSizeCalc <- function(minSample = 1,
                           maxSample = 250,
                           desiredProb = 0.8,
                           pPos = 0.5){
    pOne = calcDetectOne(minSample = minSample, maxSample = maxSample,
                         pPos = pPos)
    nSamples <- samplesForOne(minSample = minSample, maxSample = maxSample,
                              desiredProb = desiredProb,
                              pOne = pOne)
    return(nSamples)
}
sampleSizeCalc()

## This is a wrapper funciton for the USFWS UMR monitoring data. 
## it requires input from the fitted model (currently not publicly available).
## It takes the same inputs as above, but uses a specific sampling month and site. 
siteMonthWrapper <- function(month = "November",
                             site = "Boston Bay backwater",
                             data = pPositivePlot,
                             minSample = 1,
                             maxSample = 1000,
                             desiredProb = 0.8
                             ){
    
    pPos <- pPositivePlot[
        WATERBODY == site &
        MONTH == month, median]
    
    pPos_L95 <- pPositivePlot[
        WATERBODY == site &
        MONTH == month, l95]

    pPos_U95 <- pPositivePlot[
        WATERBODY == site &
        MONTH == month, u95]

    ## Find number of samples needed to have 80% detection prob
    nSamplesNeeded <- sampleSizeCalc(minSample = 1,
                                     maxSample = 1000,
                                     desiredProb = 0.8,
                                     pPos = pPos)
    
    l95 <- calcDetectOne(nSamplesNeeded, nSamplesNeeded, pPos = pPos_L95)
    u95 <- calcDetectOne(nSamplesNeeded, nSamplesNeeded, pPos = pPos_U95)

    return(list(nSamplesNeeded = nSamplesNeeded,
                l95 = l95,
                u95 = u95))
}


####################################################3
## Run code 
siteMonthWrapper(month = "November",
                 site = "Boston Bay backwater",
                 data = pPositivePlot,
                 minSample = 1,
                 maxSample = 1000,
                 desiredProb = 0.8
                 )
