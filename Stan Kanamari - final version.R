#######################################################################
# Stan version of C/W , for 5-year age groups
# Carl Schmertmann
# started 30 Oct 16
# edited 28 Apr 18 (Kanamari data loaded directly, rather than assembled from
#                   multiple files)
#######################################################################

graphics.off()
if (.Platform$OS.type == 'windows') windows(record=TRUE)
rm(list=ls() )

library(rstan)
library(dplyr)

save.fit  = FALSE
save.plot = FALSE

  #------------ STAN DATA  --------------------
  stanDataList = 
  list(
    # age pyramid data
       C    = 191, 
       W    = c(40, 34, 29, 19, 14, 9, 8),
    # local mortality estimates   
       q5_a = 3.99, 
       q5_b = 114.26, 
    # Wilmoth et al. model constants
       af = c(-0.6619,     -99, -2.5608, -3.2435, -3.1099, 
              -2.9789, -3.0185, -3.0201, -3.1487, -3.2690,
              -3.5202), 
       bf = c(0.7684,      -99,  1.7937,  1.6653,  1.5797, 
              1.5053,   1.3729,  1.2879,  1.1071,  0.9339, 
              0.6642), 
       cf = c(-0.0277,     -99,  0.1082,  0.1088,  0.1147, 
               0.1011,  0.0815,  0.0778,  0.0637,  0.0533, 
              0.0289), 
       vf = c(      0,     -99,  0.2788,  0.3423,  0.4007, 
               0.4133, 0.3884,   0.3391,  0.2829,  0.2246, 
               0.1774),
    # fertility model constants
       m =        c(   0, 1.3921, 1.5924, 1.2291, 0.4523, -0.8868, -3.4372),
       X = matrix(c(   0, 0.2744, 0.5408, 0.7319, 0.8774,  1.0405,  1.5172, 
                       0, 0.3178, 0.5108, 0.5107, 0.3483,  0.0529, -0.7236),
                    nrow=7,ncol=2)
    )  
  

  #------------ STAN MODEL --------------------
  stan.model.text = '
  data {
    int<lower=0>        C;  # observed number of children 0-4
    vector<lower=0>[7]  W;  # observed numbef of women 15-19...45-49
  
    real q5_a;              # prior will be q(5) ~ beta( q5_a, q5_b) 
    real q5_b;         
  
    real af[11];        # 11 age groups starting at x=0,1,5,10,...,45
    real bf[11];        # 11 age groups starting at x=0,1,5,10,...,45
    real cf[11];        # 11 age groups starting at x=0,1,5,10,...,45
    real vf[11];        # 11 age groups starting at x=0,1,5,10,...,45
  
    vector[7]     m;       # mean vector for gamma
    matrix[7,2]   X;       # covariance matrix for gamma
  
  }                      
  
  parameters {
    real<lower=0,upper=1> q5;    # h = log(q5) in the Wilmoth et al. system
    real                   k;    # k = Wilmoth et al. shape parameter
  
    vector[2]     beta; 
    real<lower=0> TFR; 
  
  }
  
  transformed parameters {
    vector[7]             gamma;
    real<upper=0>         h;        # log(q5)
    simplex[7]            phi;      # proportion of total fertility by age group 
    vector[8]             Fx;       # age-group fertility rates F10...F45  (F10=0)
    real<lower=0>         mx[11];   # mortality rates for age groups starting at 0,1,5,10,...45
    real<lower=0,upper=1> lx[12];   # life table {lx} values for 0,1,5...,50
    real<lower=0,upper=5> Lx[10];   # life table {5Lx} values for 0,5...,45
  
    vector[7] Kx;                   # expected surviving children 0-4 per woman 
                                    # in age group x to x+4
  
    real Kstar;                     # expected surviving total number of children               
  
  #--- child mortality index for Wilmoth model
    h = log(q5);
  
  #--- fertility rates

    gamma = m + X * beta;
    for (i in 1:7) phi[i] = exp(gamma[i]) / sum( exp(gamma));

    Fx[1] = 0;                                                # F10
    for (i in 2:8) Fx[i]   = TFR * phi[i-1] / 5;              # F15...F45 
  
  #--- mortality rates, life table survival probs, and big L values
  
    for (i in 1:11) {  mx[i]  = exp( af[i] + bf[i]*h + cf[i]*square(h) + vf[i]*k ); }
    
    mx[2]  = -0.25 * (mx[1] + log(1-q5) );               #  recalculate 1_mu_4 = -1/4 log(l[5]/l[1])

    lx[1]  = 1;                                          # x=0
    lx[2]  = lx[1] * exp(-mx[1]);                        # x=1
    lx[3]  = lx[2] * exp(-4*mx[2]);                      # x=5
    for (i in 3:12) lx[i]  = lx[i-1] * exp(-5*mx[i-1]);  # x=5,10,...50
  
    Lx[1]                  = 1* (lx[1]+lx[2])/2 + 4*(lx[2]+lx[3])/2 ;   # 5L0
    for (i in 2:10) Lx[i]  = 5* (lx[i+1]+lx[1+2])/2 ;                   # 5Lx
  
  #--- main result: expected surviving 0-4 yr olds per woman in each age group
  #      indexing is a bit complicated: 
  #       Lx[1:10] is for age groups 0,5,...,45
  #       Fx[1:8]  is for age groups 10,15,...,45
  #       Kx[1:7]  is for age groups 15,...,45
  
    for (i in 1:7) Kx[i]  = (Lx[i+2]/Lx[i+3] * Fx[i] + Fx[i+1]) * Lx[1]/2 ;
  
    Kstar  = dot_product(W, Kx);
  }
  
  
  model {
  
  # LIKELIHOOD
    C ~ poisson(Kstar);
  
  # PRIORS
    beta  ~ normal(0,1); 
  
    q5 ~ beta(q5_a, q5_b);   # 90% prior prob. that q5 is between 1/2 and 2x estimated q5
    
    k  ~ normal(0,1);
  
  }'
  
  
  
  #------------ MCMC --------------------
  stanInits = function(nchains=1) {
    L = vector('list',nchains)
    for (i in seq(L)) {
      L[[i]] =   list(
                    q5 = rbeta(n=1, shape1=3.99, shape2=114.26),
                    k = rnorm(n=1, mean=0,sd=1),
                    beta = rnorm(2),
                    TFR = runif(n=1,min=1, max=5)
                  )
    }
    return(L)
  } # stanInits
  
  
  nchains = 2
  
    fit = stan(model_code = stan.model.text,
             data       = stanDataList,
             pars       = c('h','k','beta','TFR','q5','phi','Fx','Lx','Kx','Kstar'),
             init       = stanInits(nchains),
             iter       = 10100,
             warmup     = 100,
             thin       = 4,
             chains     = nchains)
  
  
if (save.fit) {
  save(fit, file='Kanamari fit.RData')
}
    
ti = 'Kanamari do Rio Juruá'

TFRdist = summary(fit,'TFR', probs=c(0,10,25,50,75,90,100)/100)[[1]]      

TFRdist

iTFR = 7*stanDataList$C/sum(stanDataList$W)

G = stan_dens(fit,'TFR', fill='lightsteelblue',col='darkblue',alpha=.60) + 
           geom_vline(xintercept=iTFR, col='black', lwd=1, lty=2) +
           theme_bw() +
           scale_x_continuous(breaks=5:13) +
           scale_y_continuous(expand=c(0,0)) +
           geom_text(label="7*C/W = 8.74", 
                     x= iTFR+1, y=0.4, size=4,col='darkblue' ) 

G + 
  geom_point(x=TFRdist[1,'50%'], y=.01, shape=16, size=3, color='darkblue') +
  geom_segment(x=TFRdist[1,'10%'],xend=TFRdist[1,'25%'], y=.01, yend=.01, lwd=1, color='darkblue') +
  geom_segment(x=TFRdist[1,'75%'],xend=TFRdist[1,'90%'], y=.01, yend=.01, lwd=1, color='darkblue') +
  geom_text(aes(x=TFRdist[1,'10%'],y=.025), label='10%', size=3, col='darkblue')  +
  geom_text(aes(x=TFRdist[1,'25%'],y=.025), label='25%', size=3, col='darkblue')  +
  geom_text(aes(x=TFRdist[1,'50%'],y=.025), label='50%', size=3, col='darkblue')  +
  geom_text(aes(x=TFRdist[1,'75%'],y=.025), label='75%', size=3, col='darkblue')  +
  geom_text(aes(x=TFRdist[1,'90%'],y=.025), label='90%', size=3, col='darkblue') 
      
if (save.plot) {      
  ggsave(file='Stan Kanamari terra posterior density TFR.eps', device=cairo_ps,
             width=7.5, height=4.5)
}      
      