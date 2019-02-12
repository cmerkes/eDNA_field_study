library(tidyverse)
library(data.table)
library(rstan)
options(mc.cores = parallel::detectCores())

## Load in data summary and merge
d1 <- fread("30000_Data_Summary1_shared.csv")
d2 <- fread("30001_Data_Summary1_shared.csv")
d3 <- fread("30005_Data_Summary1_shared.csv")
dAll <- rbind(d1, d2, d3)

## Load in metadata and merge
m1 <- fread("Iowa_30000.csv")
m2 <- fread("Iowa_30001.csv")
m3 <- fread("Iowa_30005.csv")
mAll <- rbind(m1, m2, m3)

## Merge data summary and metadata
nrow(mAll) == nrow(dAll)
setkey(mAll, "RUID")
setkey(dAll, "UniqueID")

dataUse <- mAll[ dAll]
dataUse <- dataUse[ BLANK == "No", ] 


# dataUse <- dataUse[ WATERBODY == "Iowa River", ]

## Create ID Columns 
dataUse[ , WATERBODY := factor(WATERBODY)]
dataUse[ , WATERBODYid := as.numeric(WATERBODY)]
dataUse[ , MONTH := factor(gsub("(\\d{1,2})/(\\d{1,2})/(\\d{1,2})", "\\1",
                                DATE_COLL))]
dataUse[ , MONTHid := as.numeric(MONTH)]

dataUse[ , sampleEvent := factor(paste(WATERBODY, MONTH, sep = "-"))]
dataUse[ , sampleEventID := as.numeric(sampleEvent)]

## Create sample level detection columns
dataUse[ , A1 := ifelse(AC1.FAM.Hits > 0, 1, 0)]
dataUse[ , A3 := ifelse(AC3.HEX.Hits > 0, 1, 0)]

dataUse[ , A := ifelse(A1 > 0 | A3 > 0, 1, 0)]

dataUse[ , .(mean(TEMP_F),
             mean(DEPTH))]
        
## Center temp and depth
dataUse[ , TEMP_2 := scale(TEMP_F)]
dataUse[ , DEPTH_2 := scale(DEPTH)]


## Look at summary of data and merge in sample level Z
dataSummary <- dataUse[ , .(Prelim = sum(A), .N,
                            Abar = mean(A)),
                       by = .(  WATERBODY, WATERBODYid,
                              MONTH, MONTHid,
                              sampleEventID)][ order(sampleEventID),]
dataSummary[ , Z := ifelse(Prelim > 0, 1, 0)]
print(dataSummary)

setkey(dataSummary, "sampleEventID")
setkey(dataUse, "sampleEventID")

dataUse <- dataUse[dataSummary[ , .(sampleEventID, Z)]]


sampleEventKey <- copy(dataUse[ , .(ID = mean(sampleEventID),
                               hits = sum(AC1.FAM.Hits + AC3.HEX.Hits),
                               nPositive = sum(A)),
                               by = .(sampleEvent, MONTH, WATERBODY)])
sampleEventKey[ , MONTH := factor(MONTH, levels = c("4", "5", "11"))]

## Create site visit index
dataUse[ , index := 1:nrow(dataUse)]

start_stop <- dataUse[ , .(start_indx = min(index),
                           end_indx = max(index)),
                      sampleEventID ]

start_index <- as.array(start_stop[ , start_indx])
end_index <- as.array(start_stop[ , end_indx])

n_samples <- as.array(dataUse[ , .N, by = sampleEventID][ , N])
n_samples
head(dataUse)

## calculate number of samples per site
dataSummarySite <- dataUse[ , .(.N,
                                zObs = ifelse(sum(Z) > 0, 1, 0),
                                TEMP_F = mean(TEMP_F),
                                DEPTH = mean(DEPTH),
                                TEMP_F_2 = mean(TEMP_2), 
                                DEPTH_2 = mean(DEPTH_2)
                                ),
                           by = .(sampleEvent, sampleEventID,
                                  WATERBODY, WATERBODYid, MONTH)]

dataSummarySite
print(dataSummary, digits = 2)
site_detection <- dataSummary[ , Z]

## Create predictor matricies
Xpsi <- model.matrix( ~ WATERBODY - 1, dataSummarySite)
dim(Xpsi)
head(Xpsi)


## V is p-level coef, W is theta-level coef
head(dataUse)

Vp <- Wtheta <- model.matrix( ~  sampleEvent - 1 + TEMP_2 + DEPTH_2, dataUse)
dim(Vp)
head(Vp)

## create summary matrix for theta
head(dataSummarySite)
site_event_predict <- model.matrix( ~ sampleEvent - 1 + TEMP_F_2 + DEPTH_2, dataSummarySite)
head(site_event_predict)
dim(site_event_predict)

## Get number of sites, samples per site
n_sites <- dataUse[ , length(unique(sampleEvent))]
n_obs <- nrow(dataUse)
n_obs



## Save data in list for Stan
dataUseStanData <- list(
    n_sites = n_sites,
    n_obs = n_obs,
    n_psi_coef = ncol(Xpsi),
    X_psi = Xpsi,
    site_detection = site_detection,
    n_samples = n_samples,
    n_theta_coef = ncol(Wtheta),
    W_theta = Wtheta,
    n_p_coef = ncol(Vp),
    V_p = Vp,
    AC1 = dataUse[ , A1],
    AC3 = dataUse[ , A3],
    any_detected = dataUse[ , A],
    k = rep(8, nrow(dataUse)),
    start_index  = start_index,
    end_index = end_index,
    n_row_sep = nrow(site_event_predict),
    n_col_sep = ncol(site_event_predict),
    site_event_predic = site_event_predict,
    n_predict_samples = 200,
    delta_mean = 0,
    delta_sd = 1,
    alpha_mean = 0,
    alpha_sd = 1,
    beta_mean = 0,
    beta_sd = 1
)

## build model and sample from it 
build_model <- 
    stan_model("coefModel/coefModel.stan")

fit <- sampling(build_model,
                chains = 4,
                iter = 10000, data = dataUseStanData)

fitSummary <- summary(fit, probs = c(0.025, 0.1, 0.50, 0.9, 0.975))$summary

## Examine initial output
print(fit, pars = 'psi_site')
print(fit, pars = 'beta_psi')

print(fit, pars = 'theta_site')
print(fit, pars = 'p_AC1_site')
print(fit, pars = 'p_AC3_site')

print(fit, pars = c("beta_psi", "alpha_theta"))
print(fit, pars = c("delta_p_AC1", "delta_p_AC3"))


######################
## Examine traceplots 
traceplot(fit, pars = "psi_site", inc_warmup = TRUE)
traceplot(fit, pars = "psi_site", inc_warmup = FALSE)
traceplot(fit, pars = "theta_site", inc_warmup = TRUE)
traceplot(fit, pars = "theta_site", inc_warmup = FALSE)
traceplot(fit, pars = "p_AC1_site", inc_warmup = TRUE)
traceplot(fit, pars = "p_AC1_site", inc_warmup = FALSE)
traceplot(fit, pars = "p_AC3_site", inc_warmup = TRUE)
traceplot(fit, pars = "p_AC3_site", inc_warmup = FALSE)


traceplot(fit, pars = "lp__")

traceplot(fit, pars = "beta_psi", inc_warmup = TRUE)
traceplot(fit, pars = c("alpha_theta"), inc_warmup = TRUE)

traceplot(fit, pars = "delta_p_AC1", inc_warmup = TRUE)
traceplot(fit, pars = "delta_p_AC3", inc_warmup = TRUE)

pairs(fit, pars = c("alpha_theta", "delta_p_AC3"), inc_warmup = TRUE)


## Lookat raw data

## For Z
dataUse[ , .(Z = mean(Z)), by = .(MONTH, WATERBODY)]
dataUse[ ,  .N, by = .(Z >0)]

dataUse[ , .(Z = mean(Z)), by = .(MONTH, WATERBODY)]


## For A
dataUse[Z > 0 , .(A = mean(A)), by = .(MONTH, WATERBODY)]

## For Y
dataUse[Z > 0 & A > 0 ,
        .(pAC1 = mean(AC1.FAM.Hits/8),
          pAC3 = mean(AC3.HEX.Hits/8)),
        by = .(MONTH, WATERBODY)]

dataUse[Z > 0,
        .(pAC1 = mean(AC1.FAM.Hits/8),
          pAC3 = mean(AC3.HEX.Hits/8)),
        by = .(MONTH, WATERBODY)]

## Extract and plot nice
fit_summary <-
    fitSummary %>%
    as.data.frame() %>% 
    rownames_to_column('parameter') %>%
    as_tibble() %>%
    rename( 'l95' = `2.5%`, "l80" = `10%`,
           "median" = `50%`,
           "u80" = `90%`,  "u95" = `97.5%`)

fit_summary


## Extract out Psis to plot
pPsiPlot <-
    fit_summary %>%
    filter(grepl("psi_site", parameter)) %>%
    select(-sd) %>%
    mutate(site_id = as.numeric(gsub("psi_site\\[(\\d+)\\]", "\\1", parameter))) 

siteKey <-
    dataUse %>%
    group_by(WATERBODY) %>%
    summarize(site_id = mean(WATERBODYid))


pPsiPlot <-
    pPsiPlot %>%
    full_join(siteKey) %>%
    select(-parameter) %>%
    mutate(waterbody =
               factor(WATERBODY,
                      levels = rev(c("Dam 18 Spillway", "Boston Bay",
                                 "Iowa River", "Dam 17 Spillway")),
                      labels = rev(c("Dam 18 spillway", "Boston Bay backwater",
                                 "Iowa River tributary", "Dam 17 spillway"))))



pPsiPlot %>%
    select(waterbody, l95, median, u95) 

pPsiPlotFig <- ggplot(data = pPsiPlot, aes(x = waterbody, y = median)) +
    geom_linerange(aes(ymin = l95, ymax = u95)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4) + 
    geom_point(size = 3, shape = 15) +
    coord_flip() +
    ylim(c(0, 1))+
    ylab(expression("Estimated "*psi)) +
    xlab("Site") +
    theme_minimal()

pPsiPlotFig 
iomFolder = "./Feb_2015_figures/"
ggsave(paste0(iomFolder, "pPsiPlot.pdf"), pPsiPlotFig, width = 6, height = 4)


write.csv(file = paste0(iomFolder, "psi_prob.csv"), pPsiPlot, row.names = FALSE)

## Extract out thetas to plot
pThetaPlot <-
    fit_summary %>%
    filter(grepl("theta_site", parameter)) %>%
    select(-sd, -se_mean) %>%
    mutate(sampleEventID = as.numeric(gsub("theta_site\\[(\\d+)\\]", "\\1", parameter))) 

siteEventKey <-
    dataUse %>%
    group_by(WATERBODY, MONTH, sampleEventID) %>%
    summarize(nPositive = sum(A),
              Z = mean(Z),
              nSamples = n())


pThetaPlot <-
    pThetaPlot %>%
    full_join(siteEventKey) %>%
    select(-parameter) %>%
    filter(Z > 0) %>%
    mutate(Month= factor(MONTH,
                         levels = rev(c( 4, 5, 11)),
                         labels = rev(c( "April", "May", "November"))),
           waterbody =
               factor(WATERBODY,
                      levels = c("Dam 18 Spillway", "Boston Bay",
                                 "Iowa River", "Dam 17 Spillway"),
                      labels = c("Dam 18 spillway", "Boston Bay backwater",
                                 "Iowa River tributary", "Dam 17 spillway"))) %>%
    mutate(pPos = nPositive/nSamples)


pThetaPlot %>%
    select(waterbody, Month, l95, median, u95)


pThetaPlot %>%
    select(waterbody, Month, l95, median, u95) %>%
    arrange(Month)

write.csv(file = paste0(iomFolder, "theta_prob.csv"), pThetaPlot, row.names = FALSE)

pThetaPlotFig <-
    ggplot(data = pThetaPlot, aes(x = waterbody, y = median,
                                               color = Month)) +
    scale_size(expression("Proportion J"), trans = "sqrt") +
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 0.5)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 0.5)) + 
    geom_point(aes(size = pPos), position = position_dodge(width = 0.5), shape = 15) + 
    coord_flip() +
    ylab(expression("Estimate "*theta)) +
    xlab("Study site") +
    theme_bw() +
    scale_color_manual("Month", values = rev(c("black", "red", "blue")),
                       guide = guide_legend(reverse=TRUE)) +
    ylim(c(0, 1)) +
    facet_grid(waterbody ~ ., scales = 'free_y') +
    theme(strip.background = element_blank(), strip.text = element_blank())

pThetaPlotFig

ggsave(paste0(iomFolder, "pThetaPlot.pdf"), pThetaPlotFig, width = 6, height = 6)
ggsave(paste0(iomFolder, "pThetaPlot.jpg"), pThetaPlotFig, width = 6, height = 6)



## Extract out detectionProbs to plot
pAC1Plot <-
    fit_summary %>%
    filter(grepl("p_AC1_site", parameter)) %>%
    select(-sd, -se_mean) %>%
    mutate(sampleEventID = as.numeric(gsub("p_AC1_site\\[(\\d+)\\]", "\\1", parameter)))  %>%
    full_join(siteEventKey) %>%
    select(-parameter) %>%
    filter(Z > 0) %>%
    mutate(Month= factor(MONTH,
                         levels = rev(c( 4, 5, 11)),
                         labels = rev(c( "April", "May", "November"))),
           marker = "ACTM1") %>%
    mutate(pPos = nPositive/nSamples)

pAC3Plot <-
    fit_summary %>%
    filter(grepl("p_AC3_site", parameter)) %>%
    select(-sd, -se_mean) %>%
    mutate(sampleEventID = as.numeric(gsub("p_AC3_site\\[(\\d+)\\]", "\\1", parameter)))  %>%
    full_join(siteEventKey) %>%
    select(-parameter) %>%
    filter(Z > 0) %>%
    mutate(Month= factor(MONTH,
                         levels = rev(c( 4, 5, 11)),
                         labels = rev(c( "April", "May", "November"))),
           marker = "ACTM3") %>%
    mutate(pPos = nPositive/nSamples)

pACPlot <-
    rbind(pAC3Plot, pAC1Plot) %>%
    mutate(waterbody =
               factor(WATERBODY,
                      levels = c("Dam 18 Spillway", "Boston Bay",
                                 "Iowa River", "Dam 17 Spillway"),
                      labels = c("Dam 18 spillway", "Boston Bay backwater",
                                 "Iowa River tributary", "Dam 17 spillway")),
           marker = factor(marker,
                           levels = c("ACTM3", "ACTM1") ))

pACPlot %>%
    select(marker, waterbody, Month, l95, median, u95)


pACPlot %>%
    select(marker, waterbody, Month, l95, median, u95) %>%
    arrange(waterbody, Month)

write.csv(file = paste0(iomFolder, "p_ACTM_prob.csv"), pACPlot, row.names = FALSE)

pACPlotFig <- ggplot(data = pACPlot, aes(x = waterbody, y = median,
                                         color = Month,
                                         shape = marker)) +
    scale_size(expression("Proportion J"), trans = "sqrt") +
    scale_shape_manual("Marker", guide = guide_legend(reverse=TRUE),
                values = c(15, 18)) + 
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 0.9)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 0.9,
                   position = position_dodge( width = 0.9)) + 
    geom_point(aes(size = pPos), position = position_dodge(width = 0.9)) + 
    coord_flip() +
    ylim(c(0, 1))+
    ylab(expression("Estimate "*italic(p))) +
    xlab("Study site") +
    theme_bw() +
    scale_color_manual("Month", values = c("blue", "red", "black"),
                       guide = guide_legend(reverse=TRUE)) +
    facet_grid(waterbody ~ ., scales = 'free_y') +
    theme(strip.background = element_blank(), strip.text = element_blank())

pACPlotFig
ggsave(paste0(iomFolder, "pACPlot.pdf"), pACPlotFig, width = 6, height = 6)
ggsave(paste0(iomFolder, "pACPlot.jpg"), pACPlotFig, width = 6, height = 6)


## plot regression coef
regCoefPlot <-
    fit_summary %>%
    filter(grepl("(alpha_theta|delta_p_(AC1|AC3))\\[(13|14|15|16)\\]", parameter)) %>%
    mutate(lvl = gsub("(alpha_theta|delta_p_(AC1|AC3))\\[(13|14|15|16)\\]", "\\1", parameter),
           par = ifelse(grepl("13", parameter), "Temperature",  "Depth")) %>%
    mutate(lvl_plt = factor(lvl,
                            levels = c("delta_p_AC1", "delta_p_AC3", "alpha_theta"),
                            labels = c("ACTM1", "ACTM3","Collection sample")))

regCoefPlot %>%
    select(parameter, lvl, par, lvl_plt)

regCoefPlotFig <-
    ggplot(regCoefPlot, aes(x = par, y = median)) +
    geom_hline(yintercept = 0, color = "red", size = 2) + 
    geom_point(size = 3) +
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4) +
    geom_linerange(aes(ymin = l95, ymax = u95)) + 
    facet_grid( lvl_plt ~ .) +
    coord_flip() +
    theme_bw() +
    theme(strip.background = element_blank()) +
    xlab("Model variable coefficient") +
    ylab("Coefficient distribution on logit scale")

regCoefPlotFig
ggsave(paste0(iomFolder, "regCoefPlotFig.pdf"), regCoefPlotFig, width = 6, height = 6)
ggsave(paste0(iomFolder, "regCoefPlotFig.jpg"), regCoefPlotFig, width = 6, height = 6)


## Extract out pPrelim
pPrelimPlot <-
    fit_summary %>%
    filter(grepl("p_posistive", parameter)) %>%
    select(-sd, -se_mean) %>%
    mutate(sampleEventID = as.numeric(gsub("p_posistive\\[(\\d+)\\]", "\\1", parameter))) %>%
    full_join(siteEventKey) %>%
    select(-parameter) %>%
    filter(Z > 0) %>%
    mutate(Month= factor(MONTH,
                         levels = rev(c( 4, 5, 11)),
                         labels = rev(c( "April", "May", "November")))) %>%
    mutate(pPos = nPositive/nSamples) %>%
    mutate(waterbody =
               factor(WATERBODY,
                      levels = c("Dam 18 Spillway", "Boston Bay",
                                 "Iowa River", "Dam 17 Spillway"),
                      labels = c("Dam 18 spillway", "Boston Bay backwater",
                                 "Iowa River tributary", "Dam 17 spillway")))

           
pPrelimPlot %>%
    select(waterbody, Month, l95, median, u95)

pPrelimPlot %>%
    select(waterbody, Month, l95, median, u95) %>%
    filter(waterbody != "Boston Bay backwater") %>%
    pull(median) %>%
    range() %>%
    round(2)


write.csv(file = paste0(iomFolder, "pPrelimPlot.csv"), pPrelimPlot, row.names = FALSE)

pPrelimPlotFig <-
    ggplot(data = pPrelimPlot, aes(x = waterbody,
                                   y = median,
                                   color = Month)) +
    scale_size(expression("Proportion of samples\n(J) detected eDNA")) +
    geom_linerange(aes(ymin = l95, ymax = u95),
                   position = position_dodge( width = 0.5)) + 
    geom_linerange(aes(ymin = l80, ymax = u80), size = 1.4,
                   position = position_dodge( width = 0.5)) + 
    geom_point(aes(
        size =  pPos), 
        position = position_dodge(width = 0.5)) + 
    coord_flip() +
    ylab(expression("Probability of detecting eDNA if taking one sample")) +
    xlab("Site") +
    theme_bw() +
    scale_color_manual("Sample month",
                       values = c("blue", "red", "black"),
                       guide = guide_legend(reverse=TRUE)) +
    theme(strip.background = element_blank(),
          strip.text = element_blank()) +
    facet_grid(  waterbody ~ ., scale = "free_y") +
    ylim(c(0,1))  +
    ylab("Probability of positive detection in \u2265 1 sample") +
    xlab("Study site")
    

pPrelimPlotFig
ggsave(paste0(iomFolder, "pPrelimPlot.pdf"), pPrelimPlotFig, width = 8, height = 6)
ggsave(paste0(iomFolder, "pPrelimPlot.jpg"), pPrelimPlotFig, width = 8, height = 6)

## Extract and plot prob of detecting at least one 
prob_detect_one_Plot <-
    fit_summary %>%
    filter(grepl("prob_detect_one", parameter)) %>%
    select(-sd, -se_mean) %>%
    mutate(sampleEventID = as.numeric(gsub("prob_detect_one\\[(\\d+),(\\d+)\\]", "\\1", parameter)),
           sampleNumber = as.numeric(gsub("prob_detect_one\\[(\\d+),(\\d+)\\]", "\\2", parameter))) %>%
    full_join(siteEventKey) %>%
    select(-parameter) %>%
    filter(Z > 0) %>%
    mutate(Month= factor(MONTH,
                         levels = c( 4, 5, 11),
                         labels = c( "April", "May", "November")),
           nPositive2 = nPositive / nSamples,
           )  %>%
    mutate(waterbody =
               factor(WATERBODY,
                      levels = c("Dam 18 Spillway", "Boston Bay",
                                 "Iowa River", "Dam 17 Spillway"),
                      labels = c("Dam 18 spillway", "Boston Bay backwater",
                                 "Iowa River tributary", "Dam 17 spillway")))
           

prob_detect_one_Fig <-
    ggplot(data = prob_detect_one_Plot,
           aes(x = sampleNumber, y = median,
               color = Month, fill = Month)) +
    geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.25, color = NA) +
    geom_ribbon(aes(ymin = l80, ymax = u80), alpha = 0.25, color = NA) + 
    geom_line() +
    facet_grid(waterbody ~ Month) +
    theme_minimal() +
    scale_color_manual("Sample month", values = rev(c("blue", "red", "black")))  +
    scale_fill_manual("Sample month", values = rev(c("blue", "red", "black"))) +
    theme(legend.position="none") +
    xlab("Number of samples per site (J)") +
    ylab("Probability of positive detection in \u2265 1 sample")


prob_detect_one_Fig
ggsave(paste0(iomFolder, "prob_detect_one_Plot.pdf"), prob_detect_one_Fig, width = 8, height = 6)
ggsave(paste0(iomFolder, "prob_detect_one_Plot.jpg"), prob_detect_one_Fig, width = 8, height = 6)


