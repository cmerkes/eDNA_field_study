data {
  int<lower = 0> nObs; // number of observations 
  int<lower = 0> nPsi; // number of sites
  int<lower = 0> nTheta; // number of samples
  int<lower = 0> nPposistive; // number of positive detection probs
  int<lower = 0> AC1[nObs]; // number of AC1 detections per sample
  int<lower = 0> AC3[nObs]; // number of AC3 detections per sample
  int<lower = 0> A[nObs]; // sample level detection for AC1 OR AC3 OR Both
  int<lower = 0> Z[nObs]; // site level detection
  vector[nTheta] ZseID; // were carp detected during sampleing event (may need to change if nTheta != n sample events)
  int<lower = 0> psiID[nObs]; // dummy variable for psi, how many psi to estimate?
  int<lower = 0> thetaID[nObs]; // dummy variable for theta, how many theat to estimate?
  int<lower = 0> pID[nObs]; // dummy variable for p, how many p to estimate?
  int<lower = 0> K; // Number of PCR replicates
  int<lower = 0> nAlpha; // number of alpha parameters
  int<lower = 0> nDelta; // number of alpha parameters
  matrix[ nObs, nAlpha] W; // Predictor variables W
  matrix[ nObs, nDelta] V; // Predictor variables V
  int<lower = 0> nSiteEvents; // Number of sites per event
  matrix[ nObs, nSiteEvents] siteEventID;
  vector[nSiteEvents] nPerSampleEvent;
  int<lower = 0> nSimSamples; // Number of samples to simulate over for P(1) curve
}
parameters {
  vector[nPsi] muPpsi;
  vector[nObs] muPtheta;
  vector[nObs] muPdetectAC1;
  vector[nObs] muPdetectAC3;

  vector[nAlpha]      alpha;
  vector[nDelta]   deltaAC1;
  vector[nDelta]   deltaAC3;
 
}
model {
    
  // local variables to avoid recomputing log(psi) and log(1 - psi)
  vector[nObs] log_inv_muPpsi;
  vector[nObs] log1m_inv_muPpsi;

  vector[nObs] log_inv_muPtheta;
  vector[nObs] log1m_inv_muPtheta;

 
  for(sPsi in 1:nPsi){
    log_inv_muPpsi[sPsi]   = log_inv_logit(muPpsi[sPsi]);
    log1m_inv_muPpsi[sPsi] = log1m_inv_logit(muPpsi[sPsi]);
  }
  
  for(sT in 1:nObs){  
    log_inv_muPtheta[sT]   = log_inv_logit(muPtheta[sT]);
    log1m_inv_muPtheta[sT] = log1m_inv_logit(muPtheta[sT]);
  }

  // Priors 
  alpha    ~ normal(0, 2); 
  deltaAC1 ~ normal(0, 2);
  deltaAC3 ~ normal(0, 2);


  // regression coefficients 
  muPtheta     ~ normal( W * alpha, 2);

  muPdetectAC1 ~ normal( V * deltaAC1, 2);
  muPdetectAC3 ~ normal( V * deltaAC3, 2);
  

  // likelihood
  for (d in 1:nObs) {
    if (Z[d] > 0){ // Has DNA been found at the site?
      if (A[d] > 0) { // Has DNA been detected within a sample?
	target += // Yes, DNA is in both site and sample 
	  log_inv_muPpsi[ psiID[d]] + // yes at site
	  log_inv_muPtheta[d] + // yes in sample
	  binomial_logit_lpmf( AC1[d] | K, 
			       muPdetectAC1[d] ) + // yes detected by AC1
	  binomial_logit_lpmf( AC3[d] | K,
			       muPdetectAC3[d] ); // yes detected by AC3
      } else {  // Yes DNA is a the site, but not within sample
	target += log_sum_exp(
			      log_inv_muPpsi[psiID[d]] + // Yes at site
			      log_inv_muPtheta[d] + // Yes in sample
			      binomial_logit_lpmf( AC1[d] | K,
						   muPdetectAC1[d] ) + // Not in AC1/
			      binomial_logit_lpmf( AC3[d] | K,
						   muPdetectAC3[d] ), // Not in AC3
			      log1m_inv_muPtheta[d] ); // Or, not in sample
      }
    } else {
      target += log_sum_exp(
			    log_sum_exp( // No DNA at the site nor within sample
					log_inv_muPpsi[psiID[d]] + // Yes At site
					log_inv_muPtheta[d] + // Yes in sample 
					binomial_logit_lpmf( AC1[d] | K,
							     muPdetectAC1[d] ) + // Missed by AC1
					binomial_logit_lpmf( AC3[d] | K,
							     muPdetectAC3[d] ), // Missed by AC3
					log1m_inv_muPtheta[d]), // Or Missed by Sample
			    log1m_inv_muPpsi[psiID[d]]); // or really not in site 
    }
  }
}

generated quantities {
  vector<lower = 0, upper = 1>[nPsi] pPsi;
  vector<lower = 0, upper = 1>[nSiteEvents] pTheta;
  vector<lower = 0, upper = 1>[nSiteEvents] pDetectAC1;
  vector<lower = 0, upper = 1>[nSiteEvents] pDetectAC3;

  row_vector[nSiteEvents] muThetaSiteEvent;
  row_vector[nSiteEvents] muPAC1SiteEvent;
  row_vector[nSiteEvents] muPAC3SiteEvent;
  
  vector<lower = 0, upper = 1>[nSiteEvents] pPositive;

  matrix[ nSiteEvents, nSimSamples] pDetectOne;
  
  for(sPsi in 1:nPsi){
    pPsi[sPsi]   = inv_logit(muPpsi[sPsi]);
  }

  muThetaSiteEvent = muPtheta' * siteEventID;
  
  muPAC1SiteEvent  = muPdetectAC1' * siteEventID;
  muPAC3SiteEvent  = muPdetectAC3' * siteEventID;

  for( i in 1:nSiteEvents){
    pTheta[i]     = inv_logit(muThetaSiteEvent[i]);
    pDetectAC1[i] = inv_logit(muPAC1SiteEvent[i]);
    pDetectAC3[i] = inv_logit(muPAC3SiteEvent[i]);
  }


  for(ii in 1:nSiteEvents){
    //Calculate probability of posistive eDNA occurance for a water sample (i.e., does a water sample have eDNA?)
    pPositive[ii] = pTheta[ii]* (1.0 - (1.0 - pDetectAC1[ii])^(8.0) ) * (1.0 - ( 1.0 - pDetectAC3[ii])^(8.0));
    /* matrix[ nSiteEvents, nSimSamples] pDetectOne; */
    for(simSample in 1:nSimSamples){
      pDetectOne[ ii, simSample] = 1.0 - (1.0 - pPositive[ii]) ^ simSample;
    }
  }
}
