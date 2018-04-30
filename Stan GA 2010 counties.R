#----------------------------------------------------------------------------------
# MCMC sampling from posterior distribution of TFR, Georgia counties 2010
#----------------------------------------------------------------------------------

rm(list=ls())
graphics.off()
if (.Platform$OS.type == 'windows') windows(record=TRUE)

set.seed(6447100)

library(dplyr)
library(rstan)

save.fit = TRUE

#-------------- read data prepared in external programs
TFR.df = read.csv('GA Counties 2010 Census data.csv', stringsAsFactors = FALSE)
NCHS   = read.csv('NCHS fertility data.csv', stringsAsFactors = FALSE)
Oasis  = read.csv('GA Oasis Data.csv', stringsAsFactors = FALSE)

load(file='svd.constants.RData')
m = svd.constants$m
X = svd.constants$X

#------------- MORTALITY model ------------------

q5_hat = Oasis$q5  # estimates for 2006-2014 from GA public data website

# calculate a and b coeffs for each q5 (2 x 159)
ab = sapply(q5_hat, function(this.q) {
  LearnBayes::beta.select( list(x= this.q/2, p=.05), list(x=this.q*2, p=.95))
     })

q5_a = ab[1,]
q5_b = ab[2,]

##--- Wilmoth et al. coefficients from Pop Studies
wilmoth = 
  read.csv(text = '
           age,am,bm,cm,vm,af,bf,cf,vf
           0,  -0.5101, 0.8164,-0.0245,     0,-0.6619, 0.7684,-0.0277,     0
           1,      -99,    -99,    -99,   -99,    -99,    -99,    -99,   -99
           5,  -3.0435, 1.5270, 0.0817,0.1720,-2.5608, 1.7937, 0.1082,0.2788
           10, -3.9554, 1.2390, 0.0638,0.1683,-3.2435, 1.6653, 0.1088,0.3423
           15, -3.9374, 1.0425, 0.0750,0.2161,-3.1099, 1.5797, 0.1147,0.4007
           20, -3.4165, 1.1651, 0.0945,0.3022,-2.9789, 1.5053, 0.1011,0.4133
           25, -3.4237, 1.1444, 0.0905,0.3624,-3.0185, 1.3729, 0.0815,0.3884
           30, -3.4438, 1.0682, 0.0814,0.3848,-3.0201, 1.2879, 0.0778,0.3391
           35, -3.4198, 0.9620, 0.0714,0.3779,-3.1487, 1.1071, 0.0637,0.2829
           40, -3.3829, 0.8337, 0.0609,0.3530,-3.2690, 0.9339, 0.0533,0.2246
           45, -3.4456, 0.6039, 0.0362,0.3060,-3.5202, 0.6642, 0.0289,0.1774
           50, -3.4217, 0.4001, 0.0138,0.2564,-3.4076, 0.5556, 0.0208,0.1429
           55, -3.4144, 0.1760,-0.0128,0.2017,-3.2587, 0.4461, 0.0101,0.1190
           60, -3.1402, 0.0921,-0.0216,0.1616,-2.8907, 0.3988, 0.0042,0.0807
           65, -2.8565, 0.0217,-0.0283,0.1216,-2.6608, 0.2591,-0.0135,0.0571
           70, -2.4114, 0.0388,-0.0235,0.0864,-2.2949, 0.1759,-0.0229,0.0295
           75, -2.0411, 0.0093,-0.0252,0.0537,-2.0414, 0.0481,-0.0354,0.0114
           80, -1.6456, 0.0085,-0.0221,0.0316,-1.7308,-0.0064,-0.0347,0.0033
           85, -1.3203,-0.0183,-0.0219,0.0061,-1.4473,-0.0531,-0.0327,0.0040
           90, -1.0368,-0.0314,-0.0184,     0,-1.1582,-0.0617,-0.0259,     0
           95, -0.7310,-0.0170,-0.0133,     0,-0.8655,-0.0598,-0.0198,     0
           100,-0.5024,-0.0081,-0.0086,     0,-0.6294,-0.0513,-0.0134,     0
           105,-0.3275,-0.0001,-0.0048,     0,-0.4282,-0.0341,-0.0075,     0
           110,-0.2212,-0.0028,-0.0027,     0,-0.2966,-0.0229,-0.0041,     0
           ')

af = wilmoth$af[1:11]  # keep age 0,1,...45 
bf = wilmoth$bf[1:11]  # keep age 0,1,...45 
cf = wilmoth$cf[1:11]  # keep age 0,1,...45 
vf = wilmoth$vf[1:11]  # keep age 0,1,...45 

########################## STAN MODEL ##############################################

ncounty = nrow(TFR.df)  

## construct the matrix of constants for cumulative hazard calculations
n = c(1,5, rep(5,9))    # widths of life table age intervals for x=0,1,5,10...45
cs_constants = matrix(0, 11, 12)
for (j in 1:11) cs_constants[1:j,j+1] = head(n,j)

## construct the constants for the trapezoidal approx of L0...L45 from a row of l0,l1,l5,l10,...,l50
trapez_constants = matrix(0, 12, 10, 
                          dimnames=list(paste0('l', c(0,1,seq(5,50,5))), paste0('L', seq(0,45,5))))
trapez_constants[c('l0','l1','l5'), 'L0'] = c( 1/2, 5/2, 4/2)
for (j in 2:10) trapez_constants[j+1:2, j] = 5/2

stanDataList = list(
  n = ncounty,
  C = TFR.df$C,
  W = as.matrix( select(TFR.df, contains('W'))),
  q5_a = q5_a,
  q5_b = q5_b,
  wilmothABCV = rbind(af,bf,cf,vf),
  cs_constants = cs_constants,
  trapez_constants = trapez_constants,
  m = matrix(m, nrow=ncounty, ncol=7, byrow=TRUE),
  X = X
)

stanModelText = '
data {
  int n;              // number of counties
  int C[n];           // number of children 0-4 (both sexes)
  matrix[n,7] W;      // number of women 15-49 
  vector[n]   q5_a;   // prior will be q(5) ~ beta( q5_a, q5_b), by county 
  vector[n]   q5_b;         

  matrix[11,12]   cs_constants;     // constants for cumulative sums (mx * cs_contants) is cumul haz
  matrix[12,10]   trapez_constants; // constants for trapezoidal approx of Lx from lx

  matrix[4,11] wilmothABCV;         // Wilmoth mortality coefs for 11 ages x=0,1,5,10,...,45 
                                    // (a coefs on row 1, b on 2, ... )  

  matrix[n,7]   m;    // each row = empirical mean of gamma in the (HFD+IDB) sample of schedules (repeated n times)
  matrix[2,7]   X;    // rows contain largest 2 principal components of (gamma - mu_gamma) = beta %*% X, 
                      // scaled so that beta ~ iid N(0, 1)
  }                      

parameters {
    real<lower=0>              mu_TFR;     // average county TFR
    real<lower=0>              sigma_TFR;  // sd of county TFRs

    vector<lower=0,upper=1>[n]  q5;    // h = log(q5/1000) in the Wilmoth et al. system  PER 1000
    vector[n]                   k;     // k = Wilmoth et al. shape parameter

    matrix<lower=-2,upper=2>[n,2]   beta;   // fertility shape parameters
    vector<lower=0>[n]              TFR;    // total fertility
}
transformed parameters {
  matrix[n,7] phi;     // proportion of total fertility by age group 
  matrix[n,7] gamma;
  vector[n]   colsums_phi;

  gamma =  m + beta * X ;
  colsums_phi = exp(gamma) * rep_vector(1, 7);

  phi = diag_pre_multiply( 1 ./ colsums_phi, exp(gamma) );

}

model {

  vector[n]             h;       # log(q5)
  matrix[n,8]           Fx;      # age-group fertility rates F10...F45  (F10=0)
  matrix[n,11]          tmp;       
  matrix[n,11]          mx;      # mortality rates for age groups starting at 0,1,5,10,...45
  matrix[n,12]          lx;      # life table {lx} values for 0,1,5...,50
  matrix[n,10]          Lx;      # life table {5Lx} values for 0,5...,45
  matrix[n,7]           Sx;      # life table {5Lx} values for 0,5...,45
  
  matrix[n,7]           Kx;          # expected surviving children 0-4 per woman in age group x to x+4
  vector[n]             expected_C;  # expected surviving total number of children               

  h = log(q5);
  
  tmp     = exp( append_col( append_col( append_col( rep_vector(1.0, n), h), square(h) ), k) * wilmothABCV );
  mx      = tmp;
  mx[,2]  = -0.25 * tmp[,1] - 0.25 * log(1-q5);
  
  lx = exp( -mx * cs_constants) ;
  
  Lx = lx * trapez_constants;
  
  // numerator in Sx is Lx[, 10 ... 40], denominator is Lx[, 15 ... 45]
  Sx = block(Lx, 1, 3, n, 7 ) ./ block(Lx, 1, 4, n, 7 );  
  
  Fx = 0.20 * append_col( rep_vector(0.0, n), diag_pre_multiply( TFR, phi));
  
  for (j in 1:7) {
    Kx[,j] = (Sx[,j] .* Fx[,j] + Fx[,j+1]) .* Lx[,1] / 2. ; 
  }
  
  for (i in 1:n) expected_C[i] = dot_product( W[i,] , Kx[i,]) ;

// DISTRIBUTIONAL ASSUMPTIONS 

  # LIKELIHOOD
  C ~ poisson(expected_C);
  
  # PRIORS
  TFR ~ normal( mu_TFR, sigma_TFR);    # vectorized iid

  beta[,1] ~ normal(0, 1);
  beta[,2] ~ normal(0, 1);

  q5       ~ beta(q5_a, q5_b);   # vectorized. 90% prior prob. that q5 is between 1/2 and 2x estimated q5

  k        ~ normal(0,1);

}
'

#------------ MCMC --------------------

## for beta inits, sample randomly from the first 2 cols of the svd V matrix calculated above

stanInits = function(nchains=1) {
  L = vector('list',nchains)

  for (i in seq(L)) {
    L[[i]] =   list(
      mu_TFR    = runif(1,0,4),
      sigma_TFR = runif(1,0,2),
      q5 = rbeta(n=ncounty, shape1=q5_a, shape2=q5_b),
      k = rnorm(n=ncounty, mean=0,sd=1),
      beta = matrix(runif(ncounty*2, min=-.10, max=.10) , ncounty,2),
      TFR = pmax( .10, TFR.df$iTFR + rnorm(ncounty, 0, sd=.10))
    )
  }
  return(L)
} # stanInits

nchains = 4

  fit = stan(model_name = 'GA 2010 counties',
             model_code = stanModelText,
             data       = stanDataList,
             pars       = c('mu_TFR','sigma_TFR','TFR','beta','q5','k'),
             init       = stanInits(nchains),
             seed       = 6447100,
             iter       = 600,
             warmup     = 100,
             thin       = 1,
             chains     = nchains,
             control    = list(stepsize=0.05, max_treedepth=15)
  )

  
  now = format(Sys.time(), "%a %d%b%y %H%M")
  
  ## save the stan model results
if (save.fit) {
  save(fit,TFR.df,stanDataList, 
       file=paste0('Stan GA 2010 counties fit ', now, '.RData'))
}  
  

