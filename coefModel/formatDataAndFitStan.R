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

## Important assumption about what sample positive detection
dataUse[ , A := ifelse(A1 > 0 | A3 > 0, 1, 0)]

## Look at summary of data and merge in sample level Z
dataSummary <- dataUse[ , .(Posistive = sum(A)),
                       by = .(  WATERBODY, WATERBODYid,
                              MONTH, MONTHid,
                              sampleEventID)][ order(WATERBODY),]
## calculate Z 
dataSummary[ , Z := ifelse(Posistive > 0, 1, 0)]

## Merge dataUse and dataSummary 
setkey(dataSummary, "sampleEventID")
setkey(dataUse, "sampleEventID")
dataUse <- dataUse[dataSummary[ , .(sampleEventID, Z)]]

## Create summary table, used for Stan and Plots 
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

### Look at summary of water depths and temps
dataUse[ , mean(DEPTH) * 0.3048, by = WATERBODY]
dataUse[ , .(mean = round(mean(TEMP_F),1), min = min(TEMP_F), max = max(TEMP_F)), by = MONTH]

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

dataUse[ , Aboth := ifelse(A1 >0 & A3 >0, 1, 0)]

write.csv(file = "dataUse.csv", x = dataUse)

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
    K = dataUse[ , min(IPC.Cy5.Hits),],
    nSimSamples = 200
)
dataUse[ A >0, .(AC1 = mean(AC1.FAM.Hits), AC3 = mean(AC3.HEX.Hits)), by = .(WATERBODY, MONTH)]

#########################################################################################
#########################################################################################
## Fit model,
##
## Uncomment as needed
## Model takes ~ 45 minutes to run
## stanOut <- stan("positiveSampleCoef.stan", chains = 4,
##                 iter = 10000, data = dataUseStanData) ## use 10,000 for actual run 2588s, 100 takes ~210s 
## save(file = "stanOut.RData", stanOut)

## load saved data
load("stanOut.RData")

#########################################################################################
#########################################################################################
## Examine model summary output
summaryOut <- summary(stanOut)$summary
summaryOut

hist(summaryOut[ , "Rhat"])
## ID which parameters did not converge
summaryOut[ summaryOut[, "Rhat"] > 1.01, ] 

print(stanOut, c("pPsi", "pTheta", "pDetectAC1", "pDetectAC3", "lp__"))

## ## Examine traceplots
## traceplot(stanOut, c("alpha"), inc_warmup = FALSE)
## traceplot(stanOut, c("alpha[13]", "alpha[14]"), inc_warmup = FALSE)

## traceplot(stanOut, c("deltaAC1"))
## traceplot(stanOut, c("deltaAC3"))
## traceplot(stanOut, "pPsi")

## traceplot(stanOut, c("deltaAC1[13]", "deltaAC1[14]"), inc_warmup = FALSE)
## traceplot(stanOut, c("deltaAC3[13]", "deltaAC3[14]"), inc_warmup = FALSE)

## Lookat raw data
## For Z
dataUse[ , .(Z = mean(Z)), by = .(MONTH, WATERBODY)]
dataUse[ ,  .N, by = .(Z >0)]
dataUse[ , .(Z = mean(Z)), by = .(MONTH, WATERBODY)]

dataUse[ , .(either = round(mean(A), 2), both = round(mean(Aboth), 2)),
        by = .(MONTH, WATERBODY)]
dataUse[ , .(either = sum(A), both = sum(Aboth), A1 = sum(A1), A3 = sum(A3)),
        by = .(MONTH, WATERBODY)]

dataUse[ , .(either = round(mean(A), 2), both = round(mean(Aboth), 2)),
        by = .(WATERBODY)]

## Create tables for GIS layers
either <- dataUse[ , .(AC1 = sum(AC1.FAM.Hits),
                       AC3 = sum(AC3.HEX.Hits)), by = .(WATERBODY, MONTH)]
either
write.csv(file = "AC1_AC3_table.csv", x = either, row.names  = FALSE)
both <- dataUse[ AC1.FAM.Hits >0 & AC3.HEX.Hits, .N, by = .(WATERBODY, MONTH)]
both
write.csv(file = "posistive_table.csv", x = both, row.names  = FALSE)


## Extract and plot coefficients

## Extract out Psis to plot
pPsiPlot <- data.frame(summary(stanOut,
                   par = "pPsi",
                   probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary)
colnames(pPsiPlot)[4:8] <- c("l95", "l80", "median", "u80", "u95")

pPsiPlot$parameter <- rownames(pPsiPlot)
pPsiPlot <- data.table(pPsiPlot)
pPsiPlot[ , ID := as.numeric(gsub("pPsi\\[(\\d{1,2})\\]", "\\1", parameter))]

siteKey <- dataUse[ , .(ID = mean(WATERBODYid)), by = WATERBODY]
setkey(siteKey, "ID")
setkey(pPsiPlot, "ID")

pPsiPlot <- siteKey[ pPsiPlot]

pPsiPlot[ , WATERBODY := factor(WATERBODY, levels = levels(WATERBODY)[rev(c(2, 4, 1, 3))])]

levels(pPsiPlot$WATERBODY) <- c("Dam 18 spillway", "Boston Bay backwater", "Iowa River tributary", "Dam 17 spillway")

pPsiPlotFig <- ggplot(data = pPsiPlot, aes(x = WATERBODY, y = median)) +
    geom_linerange(aes(ymin = l95, ymax = u95)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4) + 
    geom_point(size = 3) +
    coord_flip() +
    ylim(c(0, 1))+
    ylab(expression("Estimated "*psi)) +
    xlab("Study site") +
    theme_minimal()

pPsiPlotFig 
ggsave(paste0("pPsiPlot.pdf"), pPsiPlotFig, width = 7, height = 4)
ggsave(paste0("pPsiPlot.jpg"), pPsiPlotFig, width = 7, height = 4)

## Extract out thetas to plot
pThetaPlot <- data.frame(summary(stanOut,
                   par = "pTheta",
                   probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary)
colnames(pThetaPlot)[4:8] <- c("l95", "l80", "median", "u80", "u95")

pThetaPlot$parameter <- rownames(pThetaPlot)
pThetaPlot <- data.table(pThetaPlot)
pThetaPlot[ , ID := as.numeric(gsub("pTheta\\[(\\d{1,2})\\]", "\\1", parameter))]

eventKey <- dataUse[ , .(ID = mean(sampleEventID)), by = .(WATERBODY, MONTH) ]

setkey(sampleEventKey, "ID")
setkey(pThetaPlot, "ID")

pThetaPlot <- sampleEventKey[ pThetaPlot]

levels(pThetaPlot$WATERBODY) <- c("Dam 18 spillway", "Boston Bay backwater", "Iowa River tributary", "Dam 17 spillway")

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
    xlab("Study site") +
    ylim(c(0, 1)) + 
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
ggsave(paste0( "pThetaPlot.pdf"), pThetaPlotFig, width = 7, height = 5)
ggsave(paste0( "pThetaPlot.jpg"), pThetaPlotFig, width = 7, height = 5)

## Extract out detectionProbs to plot
pAC1Plot <- data.frame(summary(stanOut,
                   par = "pDetectAC1",
                   probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary)
colnames(pAC1Plot)[4:8] <- c("l95", "l80", "median", "u80", "u95")

pAC1Plot$parameter <- rownames(pAC1Plot)
pAC1Plot <- data.table(pAC1Plot)
pAC1Plot[ , ID := as.numeric(gsub("pDetectAC1\\[(\\d{1,2})\\]", "\\1", parameter))]


pAC3Plot <- data.frame(summary(stanOut,
                   par = "pDetectAC3",
                   probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary)
colnames(pAC3Plot)[4:8] <- c("l95", "l80", "median", "u80", "u95")

pAC3Plot$parameter <- rownames(pAC3Plot)
pAC3Plot <- data.table(pAC3Plot)
pAC3Plot[ , ID := as.numeric(gsub("pDetectAC3\\[(\\d{1,2})\\]", "\\1", parameter))]


setkey(sampleEventKey, "ID")
setkey(pAC1Plot, "ID")
setkey(pAC3Plot, "ID")


pAC1Plot <- sampleEventKey[ pAC1Plot]
pAC1Plot[ , marker := "ACTM1"]

pAC3Plot <- sampleEventKey[ pAC3Plot]
pAC3Plot[ , marker := "ACTM3"]

pACPlot <- rbind(pAC3Plot, pAC1Plot)
pACPlot <- pACPlot
pACPlot
colnames(pACPlot)
pACPlot[ , .(sampleEvent, marker, l95, median, u95)]
pACPlot[ , marker := factor(marker, levels = c("ACTM3", "ACTM1"))]

levels(pACPlot$WATERBODY) <- c("Dam 18 spillway", "Boston Bay backwater", "Iowa River tributary", "Dam 17 spillway")


pACPlotFig <- ggplot(data = pACPlot, aes(x = WATERBODY, y = median,
                                         color = MONTH,
                                         shape = marker,
                                         alpha = pHits)) +
    scale_shape("Marker", breaks = c("ACTM1", "ACTM3")) +
    scale_alpha(expression("Proportion K"), breaks = c( 0, 0.075,0.15))+
    scale_size(expression("Proportion J"), breaks = c(0, 0.15, 0.3)) +
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 1)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 1)) + 
    geom_point(aes(size = pPos), position = position_dodge(width = 1)) + 
    coord_flip() +
    ylab(expression("Estimated "*italic(p))) +
    ylim(c(0, 1)) + 
    xlab("Study site") +
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
ggsave(paste0( "pACPlot.pdf"), pACPlotFig, width = 7, height = 7)
ggsave(paste0( "pACPlot.jpg"), pACPlotFig, width = 7, height = 7)


dataUse[Z >0 & A > 0,
        .(pAC1 = mean(AC1.FAM.Hits/8),
          pAC3 = mean(AC3.HEX.Hits/8)),
        by = .(MONTH, WATERBODY)]

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
pPositivePlot <-  data.frame(summary(stanOut,
                                     par = "pPositive",
                                     probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary)

colnames(pPositivePlot)[4:8] <- c("l95", "l80", "median", "u80", "u95")

head(pPositivePlot)
pPositivePlot$rawStanID <- rownames(pPositivePlot)

pPositivePlot <- data.table(pPositivePlot)
pPositivePlot
pPositivePlot[ , ID := as.numeric(gsub("(pPositive)\\[(\\d{1,2})\\]", "\\2", rawStanID))]

setkey( sampleEventKey, "ID")
setkey(pPositivePlot, "ID")

pPositivePlot<- sampleEventKey[ pPositivePlot]
head(pPositivePlot)

tail(pPositivePlot)
levels(pPositivePlot$WATERBODY) <- c("Dam 18 spillway", "Boston Bay backwater", "Iowa River tributary", "Dam 17 spillway")

pPositivePlotFig <- ggplot(data = pPositivePlot, aes(x = WATERBODY, y = median,
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
    ylim(c(0, 1))+
    ylab(expression("Probability of positive detection in \u2265 1 sample")) +
    xlab("Study site") +
    scale_color_manual("Month", values = c("blue", "red", "black"),
				       breaks = c("April", "May", "November"))+
    facet_grid( WATERBODY ~., scales = 'free') +
    theme_bw() +
	theme(
  		strip.background = element_blank(),
		strip.text.x = element_blank(),
		strip.text.y = element_blank()
		)

pPositivePlotFig
ggsave(paste0( "pPositivePlot.pdf"), pPositivePlotFig, width = 8, height = 6)
ggsave(paste0( "pPositivePlot.jpg"), pPositivePlotFig, width = 8, height = 6)

##########################################################################
## Extract and plot regression coefficients 
alphaNames <- gsub("sampleEvent", "", colnames(W))

alphaPlot <- data.frame(summary(stanOut,
                                par = "alpha",
                                probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary)
colnames(alphaPlot)[4:8] <- c("l95", "l80", "median", "u80", "u95")

alphaPlot$parameter <- rownames(alphaPlot)
alphaPlot <- data.table(alphaPlot)
alphaPlot[ , ID := alphaNames]


ggplot(alphaPlot, aes(x = ID, y = median)) +
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 0.5)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 0.5)) + 
    geom_point() + 
    coord_flip() 

deltaAC1Names <- gsub("sampleEvent", "", colnames(W))

deltaAC1Plot <- data.frame(summary(stanOut,
                                par = "deltaAC1",
                                probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary)
colnames(deltaAC1Plot)[4:8] <- c("l95", "l80", "median", "u80", "u95")

deltaAC1Plot$parameter <- rownames(deltaAC1Plot)
deltaAC1Plot <- data.table(deltaAC1Plot)
deltaAC1Plot[ , ID := deltaAC1Names]

ggplot(deltaAC1Plot, aes(x = ID, y = median)) +
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 0.5)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 0.5)) + 
    geom_point() + 
    coord_flip() 



deltaAC3Names <- gsub("sampleEvent", "", colnames(W))

deltaAC3Plot <- data.frame(summary(stanOut,
                                par = "deltaAC3",
                                probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary)
colnames(deltaAC3Plot)[4:8] <- c("l95", "l80", "median", "u80", "u95")

deltaAC3Plot$parameter <- rownames(deltaAC3Plot)
deltaAC3Plot <- data.table(deltaAC3Plot)
deltaAC3Plot[ , ID := deltaAC3Names]


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
coefPlot[  , .(l95 = exp(l95), median = exp(median), u95 = exp(u95)), by = .(ID, level)]

coefPlot[ , level := factor(level)]
levels(coefPlot$level) <- c("ACTM1", "ACTM3", "Collection sample")

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
    xlab("Model variable coefficients") +
    theme_minimal()

print(coefGG)

ggsave("coefGG.pdf", coefGG, width = 6, height = 4)
ggsave("coefGG.jpg", coefGG, width = 6, height = 4)

## Probability of detection

rownames(summaryOut)

pPos <- data.frame(summary(stanOut, par = "pDetectOne", probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary)

colnames(pPos)[4:8] <- c("l95", "l80", "median", "u80", "u95")
head(pPos)
pPos$stanName <- rownames(pPos)
pPos <- data.table(pPos)
pPos

pPos[ , site := as.numeric(gsub("(pDetectOne\\[)(\\d{1,2}),(\\d{1,3})\\]", "\\2", stanName))]
pPos[ , nSamples := as.numeric(gsub("(pDetectOne\\[)(\\d{1,2}),(\\d{1,3})\\]", "\\3", stanName))]
pPos

sampleEventKey
setkey( sampleEventKey, "ID")
setkey( pPos, "site")

pPosPlot <- sampleEventKey[ pPos]
pPosPlot
levels(pPosPlot$WATERBODY) <- c("Dam 18 spillway", "Boston Bay backwater",
                                      "Iowa River tributary", "Dam 17 spillway")

pPosPlot[ , MONTH := factor(MONTH)]
pPosPlot[ , MONTH := factor(MONTH, levels = c("April", "May", "November"))]
pPosPlot
library(scales)
pPosPlot
detectOne <- ggplot(pPosPlot[ nPositive > 0, ], aes(x = nSamples, y = median, color = MONTH, fill = MONTH)) +
    facet_grid(WATERBODY ~ MONTH)  + theme_minimal() +
    scale_color_manual("Month", values = c("black", "red", "blue")) + 
    scale_fill_manual("Month", values = c("black", "red", "blue")) + 
    geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.25, color = NA) +
    geom_ribbon(aes(ymin = l80, ymax = u80), alpha = 0.25, color = NA) +
    geom_line() +
    xlab("Number of samples per site collected") +
    ylab( "Probability of detecting \u22651 positive sample per site")
print(detectOne)

ggsave("detectOne.jpg", detectOne, width = 6, height = 6)
ggsave("detectOne.pdf", detectOne, width = 6, height = 6)

