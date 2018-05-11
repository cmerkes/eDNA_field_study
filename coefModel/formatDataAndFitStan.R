## Load required packages and set options 
library(data.table)
library(rstan)
options(mc.cores = parallel::detectCores())
options(width = 100)


#########################################################################################
#########################################################################################
## Data wrangling section
## Load in data summary and merge
d1 <- fread("../30000_Data_Summary1_shared.csv")
d2 <- fread("../30001_Data_Summary1_shared.csv")
d3 <- fread("../30005_Data_Summary1_shared.csv")
dAll <- rbind(d1, d2, d3)

## Load in metadata and merge
m1 <- fread("../Iowa_30000.csv")
m2 <- fread("../Iowa_30001.csv")
m3 <- fread("../Iowa_30005.csv")
mAll <- rbind(m1, m2, m3)

## Merge data summary and metadata
nrow(mAll) == nrow(dAll)
setkey(mAll, "RUID")
setkey(dAll, "UniqueID")
dataUse <- mAll[ dAll]

## Remove blanks 
dataUse <- dataUse[ BLANK == "No", ] 

## Create ID Columns 
dataUse[ , WATERBODY := factor(WATERBODY)]
dataUse[ , WATERBODYid := as.numeric(WATERBODY)]
dataUse[ , MONTH := factor(gsub("(\\d{1,2})/(\\d{1,2})/(\\d{1,2})", "\\1", DATE_COLL))]
dataUse[ , MONTHid := as.numeric(MONTH)]
dataUse[ , sampleEvent := factor(paste(WATERBODY, MONTH, sep = "-"))]
dataUse[ , sampleEventID := as.numeric(sampleEvent)]

## Create sample level detection columns
dataUse[ , A1 := ifelse(AC1.FAM.Hits > 0, 1, 0)]
dataUse[ , A3 := ifelse(AC3.HEX.Hits > 0, 1, 0)]

## Important assumption
dataUse[ , A := ifelse(A1 > 0 | A3 > 0, 1, 0)]

## Look at summary of data and merge in sample level Z
dataSummary <- dataUse[ , .(Posistive = sum(A)),
                       by = .(  WATERBODY, WATERBODYid,
                              MONTH, MONTHid,
                              sampleEventID)][ order(WATERBODY),]
dataSummary[ , Z := ifelse(Posistive > 0, 1, 0)]
print(dataSummary)

setkey(dataSummary, "sampleEventID")
setkey(dataUse, "sampleEventID")

dataUse <- dataUse[dataSummary[ , .(sampleEventID, Z)]]


sampleEventKey <- copy(dataUse[ , .(ID = mean(sampleEventID),
                               hits = sum(AC1.FAM.Hits + AC3.HEX.Hits),
                               nPositive = sum(A),
                               nSites = .N),
                               by = .(sampleEvent, MONTH, WATERBODY)])
                               

sampleEventKey[ , MONTH := factor(MONTH, levels = c('11', '5', '4'), 
                                  labels = c('November', "May", 'April'))]

sampleEventKey[ , WATERBODY := factor(WATERBODY, levels = levels(WATERBODY)[rev(c(2, 4, 1, 3))])]


sampleEventKey[ , pPos := nPositive / nSites ]
sampleEventKey[ , pHits := hits / (nSites * 8)]

colnames(dataUse)


#######################################################
## Coefficients for predictions 
## Note that X is not used because we do not have enough replication at that level
## X = dataUse[ , model.matrix( ~ WATERBODY - 1 )]

## Specify equations for model 
W = dataUse[ , model.matrix ( ~ sampleEvent - 1 + DEPTH + TEMP_F )]
V = dataUse[ , model.matrix ( ~ sampleEvent - 1 + DEPTH + TEMP_F)]

## site-event matrix also includes the weighting for numberof samples per site 
siteEventID <- dataUse[ , model.matrix ( ~ sampleEvent - 1)]
nPerSampleEventIDdt <- dataUse[ , 1.0 / .N, by = sampleEvent]
siteEventID <- siteEventID * nPerSampleEventIDdt[ , V1]

head(dataUse)
dataUse[ A >0,]
dataUse[ , Aboth := ifelse(A1 >0 & A3 >0, 1, 0)]

write.csv(file = "dataUse.csv", x = dataUse)
dataUse[ , summary(WATERBODY)]
dataUse[ , summary(MONTH)]

summary(W[, "TEMP_F"])
summary(V[, "DEPTH"])

## Create list of data for Stan
dataUseStanData <- list(
    siteEventID = siteEventID,
    nSiteEvents = dim(siteEventID)[2],
    nPerSampleEvent = nPerSampleEventIDdt[ , V1],
    W = W,
    nAlpha = dim(W)[2],
    V = V,
    nDelta = dim(V)[2],
    nObs = nrow(dataUse),
    nPsi = dataUse[ , length(unique(WATERBODY))],
    nTheta = dataUse[ , length(unique(sampleEventID))],
    nPposistive = dataUse[ , length(unique(sampleEventID))],
    AC1 = dataUse[ , AC1.FAM.Hits],
    AC3 = dataUse[ , AC3.HEX.Hits],
    A = dataUse[ , A],
    Z = dataUse[ , Z],
    ZseID = dataUse[, .(mean(Z >  0)), by = sampleEventID][ , V1],
    psiID = dataUse[ ,  WATERBODYid],
    thetaID = dataUse[ , sampleEventID],
    pID = dataUse[ ,sampleEventID],
    K = dataUse[ , min(IPC.Cy5.Hits)]
)
dataUse[ A >0, .(AC1 = mean(AC1.FAM.Hits), AC3 = mean(AC3.HEX.Hits)), by = .(WATERBODY, MONTH)]



### Fit model, uncomment as needed
## stanOut <- stan("positiveSampleCoef.stan", chains = 4,
##                 iter = 10000, data = dataUseStanData)
## save(file = "stanOut.RData", stanOut)


load("stanOut.RData")
summaryOut <- summary(stanOut)$summary
hist(summaryOut[ , "Rhat"])

## print(stanOut, c("alpha", "pPsi", "deltaAC1", "deltaAC3", "lp__"))
print(stanOut, c("pPsi", "pTheta", "pDetectAC1", "pDetectAC3", "lp__"))

traceplot(stanOut, c("alpha"), inc_warmup = FALSE)
traceplot(stanOut, c("alpha[13]", "alpha[14]"), inc_warmup = FALSE)

traceplot(stanOut, c("deltaAC1"))
traceplot(stanOut, c("deltaAC3"))
traceplot(stanOut, "pPsi")

traceplot(stanOut, c("deltaAC1[13]", "deltaAC1[14]"), inc_warmup = FALSE)
quartz()
traceplot(stanOut, c("deltaAC3[13]", "deltaAC3[14]"), inc_warmup = FALSE)

## Lookat raw data
## For Z
dataUse[ , .(Z = mean(Z)), by = .(MONTH, WATERBODY)]
dataUse[ ,  .N, by = .(Z >0)]
dataUse[ , .(Z = mean(Z)), by = .(MONTH, WATERBODY)]

dataUse[ WATERBODY == "Iowa River" & MONTH == '5' & A > 0,]
dataUse[ , .(A = sum(A)), by = .(MONTH, WATERBODY)]

colnames(dataUse)

either <- dataUse[ , .(AC1 = sum(AC1.FAM.Hits),
                       AC3 = sum(AC3.HEX.Hits)), by = .(WATERBODY, MONTH)]
either
write.csv(file = "AC1_AC3_table.csv", x = either, row.names  = FALSE)
both <- dataUse[ AC1.FAM.Hits >0 & AC3.HEX.Hits, .N, by = .(WATERBODY, MONTH)]
both
write.csv(file = "posistive_table.csv", x = both, row.names  = FALSE)

## For A
dataUse[Z > 0 , .(A = mean(A)), by = .(MONTH, WATERBODY)]

## For Y
OrCompare <- dataUse[A1 > 0 | A3 >0 ,
        .(pAC1 = mean(AC1.FAM.Hits/8),
          pAC3 = mean(AC3.HEX.Hits/8)),
        by = .(MONTH, WATERBODY)]
OrCompare

mean(OrCompare[ , dbinom( 0, 8, pAC1)])
mean(OrCompare[ , dbinom( 0, 8, pAC3)])


AndCompare <- dataUse[A1 > 0 & A3 >0 ,
        .(pAC1 = mean(AC1.FAM.Hits/8),
          pAC3 = mean(AC3.HEX.Hits/8)),
        by = .(MONTH, WATERBODY)]
AndCompare

A1Only <- dataUse[A1 > 0,
        .(pAC1 = mean(AC1.FAM.Hits/8)),
        by = .(MONTH, WATERBODY)]


dataUse[ A1 != A3, .N, by = .(A1 >0, MONTH, WATERBODY)]
dataUse[ , .N, by = .(A1 >0, MONTH, WATERBODY)]


A1Only
A1Only[ , dbinom( 0, 8, pAC1)]

A3Only <- dataUse[A3 > 0 ,
        .(pAC3 = mean(AC3.HEX.Hits/8)),
        by = .(MONTH, WATERBODY)]

A3Only[ , dbinom( 0, 8, pAC3)]


ggplot(dataUse[Z > 0 & A1 > 0,], aes(x = AC1.FAM.Hits)) + geom_histogram() + facet_grid(MONTH ~ WATERBODY)
## Extract and plot nice


extractEst <- function(dataIn = stanOut,
                       parameter = "pPsi",
                       rowNames = NULL
                       ){
    dataOut <-  data.table(t(apply(extract(dataIn, pars = c(parameter))[[1]],
                                    2, quantile, c(0.975, 0.9, 0.5, 0.1, 0.025))))
    
    setnames( dataOut, c("u95", "u80", "median", "l80", "l95"))
    if(is.null(rowNames)){
        dataOut[ , ID := 1:.N]
    } else {
        dataOut[ , ID := rowNames]
    }
    return(dataOut)
}

## Extract out Psis to plot
pPsiPlot <- extractEst(dataIn = stanOut,
                       parameter = "pPsi")

pPsiPlot
siteKey <- dataUse[ , .(ID = mean(WATERBODYid)), by = WATERBODY]
setkey(siteKey, "ID")
setkey(pPsiPlot, "ID")

pPsiPlot <- siteKey[ pPsiPlot]

pPsiPlot[ , WATERBODY := factor(WATERBODY, levels = levels(WATERBODY)[rev(c(2, 4, 1, 3))])]

pPsiPlotFig <- ggplot(data = pPsiPlot, aes(x = WATERBODY, y = median)) +
    geom_linerange(aes(ymin = l95, ymax = u95)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4) + 
    geom_point(size = 3) +
    coord_flip() +
    ylim(c(0, 1))+
    ylab(expression("Estimated "*psi)) +
    xlab("Site") +
    theme_minimal()

pPsiPlotFig 
ggsave(paste0("pPsiPlot.pdf"), pPsiPlotFig, width = 6, height = 4)
ggsave(paste0("pPsiPlot.jpg"), pPsiPlotFig, width = 6, height = 4)

## Extract out thetas to plot
pThetaPlot <- extractEst(dataIn = stanOut,
                       parameter = "pTheta")
pThetaPlot[ , ID := 1:.N]
dataUse
eventKey <- dataUse[ , .(ID = mean(sampleEventID)), by = .(WATERBODY, MONTH) ]
eventKey

setkey(sampleEventKey, "ID")
setkey(pThetaPlot, "ID")

pThetaPlot <- sampleEventKey[ pThetaPlot]

pThetaPlotFig <- ggplot(data = pThetaPlot, aes(x = WATERBODY, y = median,
                                               color = MONTH,
                                               alpha = pHits)) +
    scale_alpha(expression("Proportion K"), breaks = c( 0, 0.075,0.15))+
    scale_size(expression("Proportion J"), breaks = c(0, 0.15, 0.3)) +
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 0.5)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 0.5)) + 
    geom_point(aes(size = pPos), position = position_dodge(width = 0.5)) + 
    coord_flip() +
    ylab(expression("Estimated "*theta)) +
    xlab("Site") +
    theme_bw() +
    scale_color_manual("Month", values = c("blue", "red", "black"),
				       breaks = c("April", "May", "November")) +
	facet_grid( WATERBODY ~., scales = 'free') +
	theme(
  		strip.background = element_blank(),
		strip.text.x = element_blank(),
		strip.text.y = element_blank()
		)

pThetaPlotFig
ggsave(paste0( "pThetaPlot.pdf"), pThetaPlotFig, width = 6, height = 5)
ggsave(paste0( "pThetaPlot.jpg"), pThetaPlotFig, width = 6, height = 5)


## Extract out detectionProbs to plot
pAC1Plot <- extractEst(dataIn = stanOut,
                       parameter = "pDetectAC1")
pAC3Plot <- extractEst(dataIn = stanOut,
                       parameter = "pDetectAC3")

pAC1Plot[ , ID := 1:.N]
pAC3Plot[ , ID := 1:.N]

setkey(sampleEventKey, "ID")
setkey(pAC1Plot, "ID")
setkey(pAC3Plot, "ID")


pAC1Plot <- sampleEventKey[ pAC1Plot]
pAC1Plot[ , marker := "AC1"]

pAC3Plot <- sampleEventKey[ pAC3Plot]
pAC3Plot[ , marker := "AC3"]

pACPlot <- rbind(pAC3Plot, pAC1Plot)
pACPlot <- pACPlot
pACPlot

pACPlot[ , marker := factor(marker, levels = c("AC3", "AC1"))]

pACPlotFig <- ggplot(data = pACPlot, aes(x = WATERBODY, y = median,
                                         color = MONTH,
                                         shape = marker,
                                         alpha = pHits)) +
    scale_shape("Marker", breaks = c("AC1", "AC3")) +
    scale_alpha(expression("Proportion K"), breaks = c( 0, 0.075,0.15))+
    scale_size(expression("Proportion J"), breaks = c(0, 0.15, 0.3)) +
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 1)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 1)) + 
    geom_point(aes(size = pPos), position = position_dodge(width = 1)) + 
    coord_flip() +
    ylab(expression("Estimated "*italic(p))) +
    xlab("Site") +
    theme_bw() +
    scale_color_manual("Month", values = c("blue", "red", "black"),
                       breaks = c("April", "May", "November"))+
    facet_grid( WATERBODY ~., scales = 'free') +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank()
    )



pACPlotFig
ggsave(paste0( "pACPlot.pdf"), pACPlotFig, width = 6, height = 7)
ggsave(paste0( "pACPlot.jpg"), pACPlotFig, width = 6, height = 7)


dataUse[Z >0 & A > 0,
        .(pAC1 = mean(AC1.FAM.Hits/8),
          pAC3 = mean(AC3.HEX.Hits/8)),
        by = .(MONTH, WATERBODY)]

head(dataUse)
dataUse[ ,
        .(pAC1 = sum(AC1.FAM.Hits),
          pAC3 = sum(AC3.HEX.Hits)),
        by = .(MONTH, WATERBODY)]

dataUse[ MONTH == 11 ,
        .(N = .N),
        by = .(MONTH, WATERBODY, A1>0, A3 > 0)]


dataUse[ WATERBODY == "Dam 17 Spillway" & A > 0  & MONTH == 11,]
dataUse[ WATERBODY == "Boston Bay"  & MONTH == 11,]

        
## Extract out pPosistive
pPosistivePlot <- extractEst(dataIn = stanOut,
                          parameter = "pPosistive")
pPosistivePlot[ , ID := 1:.N]

setkey( sampleEventKey, "ID")
setkey(pPosistivePlot, "ID")

pPosistivePlot <- sampleEventKey[ pPosistivePlot]
head(pPosistivePlot)
tail(pPosistivePlot)

pPosistivePlotFig <- ggplot(data = pPosistivePlot, aes(x = WATERBODY, y = median,
                                                 color = MONTH,
                                                 alpha = pHits)) +
    scale_alpha(expression("Proportion K"))+
    scale_size(expression("Proportion J")) +
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 0.5)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 0.5)) + 
    geom_point(aes(
        size =  pPos), 
               position = position_dodge(width = 0.5)) + 
    coord_flip() +
    ## ylim(c(0, 1))+
    ylab(expression("Probability of positive detection in one sample")) +
    xlab("Site") +
    scale_color_manual("Month", values = c("blue", "red", "black"),
				       breaks = c("April", "May", "November"))+
    facet_grid( WATERBODY ~., scales = 'free') +
    theme_bw() +
	theme(
  		strip.background = element_blank(),
		strip.text.x = element_blank(),
		strip.text.y = element_blank()
		)

## x11()
pPosistivePlotFig
ggsave(paste0( "pPosistivePlot.pdf"), pPosistivePlotFig, width = 8, height = 6)
ggsave(paste0( "pPosistivePlot.jpg"), pPosistivePlotFig, width = 8, height = 6)

##########################################################################
## Extract and plot regression coefficients 
alphaNames <- gsub("sampleEvent", "", colnames(W))
alphaPlot  <- extractEst(dataIn = stanOut, parameter = "alpha", rowNames = alphaNames)
alphaPlot
ggplot(alphaPlot, aes(x = ID, y = median)) +
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 0.5)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 0.5)) + 
    geom_point() + 
    coord_flip() 

deltaAC1Names <- gsub("sampleEvent", "", colnames(V))
deltaAC1Plot  <- extractEst(dataIn = stanOut, parameter = "deltaAC1", rowNames = deltaAC1Names)
print(deltaAC1Plot)

ggplot(deltaAC1Plot, aes(x = ID, y = median)) +
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 0.5)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 0.5)) + 
    geom_point() + 
    coord_flip() 


deltaAC3Names <- gsub("sampleEvent", "", colnames(V))
deltaAC3Plot  <- extractEst(dataIn = stanOut, parameter = "deltaAC3", rowNames = deltaAC3Names)
print(deltaAC3Plot)

ggplot(deltaAC3Plot, aes(x = ID, y = median)) +
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 0.5)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 0.5)) + 
    geom_point() + 
    coord_flip() 


alphaPlot[ , level := "Sample"]
deltaAC1Plot[ , level := "AC1"]
deltaAC3Plot[ , level := "AC3"]

allPlot <- rbind(alphaPlot, deltaAC1Plot, deltaAC3Plot)
coefPlot <- allPlot[ grep("DEPTH|TEMP", ID),]

coefPlot[ , ID := factor(ID, labels = c("Depth", "Temperature"))]

coefGG <- ggplot( coefPlot, aes(x = ID, y = median)) +
    geom_hline(yintercept = 0, color = 'red', size = 1.5) + 
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 0.5)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 0.5)) + 
    geom_point() +
    facet_grid(  level ~ . ) + 
    coord_flip() +
    ylab("Estimated distribution on raw logit scale") +
    xlab("Model covariate coefficients") + theme_minimal()
print(coefGG)

ggsave("coefGG.pdf", coefGG, width = 6, height = 4)
ggsave("coefGG.jpg", coefGG, width = 6, height = 4)

depthGG <- ggplot(dataUse, aes(x = DEPTH, y = AC1.FAM.Hits + AC3.HEX.Hits)) + geom_point() +
    stat_smooth()+ facet_grid( ~ WATERBODY)
depthGG
ggsave("depthGG.pdf", depthGG, width = 6, height = 4)

dataUse[ , summary(TEMP_F)]
tempMargin <- seq(39, 64, by = 1)

alphaMarginPlot <- copy(alphaPlot)
alphaMarginPlot[ , tempSlope := alphaMarginPlot[ ID == "TEMP_F", median]]
alphaMarginPlot <- copy(alphaMarginPlot[ ID != 'TEMP_F' & ID != 'DEPTH',
                                        .(median, ID, tempSlope)])


alphaTempValues = data.table(expand.grid(ID = alphaMarginPlot[ , ID], temp = tempMargin))

setkey(alphaTempValues, "ID")
setkey(alphaMarginPlot, "ID")
alphaGGplot <- alphaTempValues[ alphaMarginPlot]
invLogit <- function(x){ return(exp(x) / (1 + exp(x)))}

alphaGGplot[ , predicted := invLogit( median + temp * tempSlope)]

dataSummary[ , ID := paste( WATERBODY, MONTH, sep = "-")] 

setkey(dataSummary, "ID")
setkey(alphaGGplot, 'ID')
alphaGGplot <- dataSummary[ alphaGGplot]

tempMarginPlot <- ggplot(alphaGGplot, aes(x = temp, y = predicted, color = MONTH)) +
    geom_jitter(data = dataUse, aes(x = TEMP_F, y = A, color = MONTH),
                width = 0, height = 0.1)  + 
    geom_line() +
    facet_grid( ~ WATERBODY) + ylab(expression("estiamted "*theta)) +
    xlab(expression("Temperature ("*degree*"F)")) + theme_minimal() +
    scale_color_manual("Month", values = c("red", "blue", "black"))

print(tempMarginPlot)
ggsave("TemperatureMarginPlot.pdf", tempMarginPlot, width = 6, height = 4)


ggplot(data = dataUse, aes(x = TEMP_F, y = A, color = MONTH))  +
    geom_jitter(width = 0, height = 0.1) + 
    facet_grid( WATERBODY ~ MONTH, scales = "free_x") +
    ylab(expression("estiamted "*theta)) +
    xlab(expression("Temperature ("*degree*"F)")) + theme_bw() +
    scale_color_manual("Month", values = c("red", "blue", "black")) + 
    stat_smooth(method = 'glm', method.args = list(family = "binomial"))
