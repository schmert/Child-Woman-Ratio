##############################################################
# approximate posterior of TFR by simulation, as
# described in Appendix
# 
# Caution: this program has a long run time. 
##############################################################
graphics.off()
if (.Platform$OS.type == 'windows')  windows(record=TRUE)

rm(list=ls())

library(dplyr)
library(ggplot2)
library(HIV.LifeTables)

#############################################################
delta  = .005
TFRvals = seq(0.5+delta/2,  12.5-delta/2, delta)  # grid of possible TFR values
############################################################

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

############################################################
# grab matrices of empirical LX and PHIX schedules 
# from HMD and HFD
############################################################

## construct PHIX matrix (2054 x 7) of HFD age patterns

HFD = read.table(file='asfrRR.txt', skip=2, as.is=TRUE, header=TRUE)
HFD$Age = 12:55
HFD = HFD %>%
       filter(Age %in% 15:45) %>%
       mutate(age.group= 5 * (Age %/% 5)) %>%
       group_by(Code,Year,age.group) %>%
       summarize(Fx=mean(ASFR)) %>%
       mutate(phix = Fx/sum(Fx), sched = paste0(Code,Year)) %>%
       ungroup

HFD_PHIX = matrix(HFD$phix, ncol=7, byrow=TRUE,
             dimnames=list(unique(HFD$sched), paste0('phi',seq(15,45,5))))

## construct L matrix ( x ) of HMD Lx values

#  all of the 5x5 period life tables (both sexes) from HMD, AUS.bltper_5x5.txt, etc.
file.names = paste0('bltper_5x5/', dir('./bltper_5x5'))  

Lvalues = data.frame()

q5_limit = 1

for (this.file in file.names) {
  tmp = read.table(this.file, skip=2, stringsAsFactors = FALSE, header=TRUE) %>%
          group_by(Year) %>%
          summarize(Code=substr(this.file,14,16),
                    L0  = Lx[1]+Lx[2],
                    L5  = Lx[3],
                    L10 = Lx[4],
                    L15 = Lx[5],
                    L20 = Lx[6],
                    L25 = Lx[7],
                    L30 = Lx[8],
                    L35 = Lx[9],
                    L40 = Lx[10],
                    L45 = Lx[11],
                    q5  = 1-lx[3]/100000) %>%
          mutate(sched = paste0(Code,Year)) %>%
          ungroup %>%
          filter(q5 <= q5_limit)
  
  Lvalues = rbind(Lvalues, tmp)

}  

HMD_LX = as.matrix(Lvalues[,3:13], dimnames=list(NULL, names(Lvalues[3:13])))
rownames(HMD_LX) = Lvalues$sched


############################################################
# fit Dirichlet parameters for International Database schedules
############################################################
library(dirmult)

IDB = read.csv(file='http://schmert.net/calibrated-spline/data/IDB%205fx.csv',
               stringsAsFactors = FALSE) 

X = as.matrix(IDB[, paste0('f', seq(15,45,5))])

X = round(10000 * prop.table(X,1))

fit = dirmult(X, epsilon=1e-9)

IDBgamma = fit$gamma   # 7x1 parameter vector for Dirichlet simulation


#####################################################################
# functions to draw phi and L vectors, depending on prior
#####################################################################
load("svd.constants.RData")

m = svd.constants$m
XMAT = t(svd.constants$X)  # note transpose: this X is 7x2

## returns an nx7 matrix of random phi values (rowSums = 1)
draw_phi = function(nsamp=1, prior='text') {
  if (prior=='text') {
    beta = matrix(rnorm(2*nsamp), nrow=2)
    gamma = m + XMAT %*% beta
    phi   = prop.table( exp(t(gamma)), 1)
  } else
  if (prior=='empirical') {
    sel = sample(nrow(HFD_PHIX), nsamp, replace=TRUE)
    phi = HFD_PHIX[sel,]
  } else 
  if (prior=='dirichlet') {  
     phi = rdirichlet(n=nsamp, IDBgamma)    
  } 
  
  return( matrix(phi, nrow=nsamp, 
            dimnames=list(NULL,paste0('phi',seq(15,45,5)))))
} # draw_phi


## returns an nx10 matrix of Lx values
draw_L = function(nsamp=1, prior='text', q5min, q5max,
                  prev.min=0, prev.max=5) {
  
    ab = LearnBayes::beta.select(list(p=.05,x=q5min/2),
                                 list(p=.95,x=q5max*2))
    
    if (prior=='text') {
      Lx = matrix(NA, nsamp, 10, 
                  dimnames=list(NULL, paste0('L',seq(0,45,5))))
      this.h = log( rbeta(nsamp,ab[1],ab[2]))
      this.k  = rnorm(nsamp,0,1)
      for (i in seq(this.h)) {
         logmx = af + bf*this.h[i] + cf*this.h[i]^2 + vf*this.k[i]
         mx = exp(logmx)
         mx[2] = -1/4 * (mx[1] + log(1-exp(this.h[i])))
         
         x = c(0, 1, seq(5,45,5))
         n = c(diff(x),5)
         lx = c(1, exp(-cumsum(n*mx)))
         ax = c( .07+1.5*mx[1], 1.5, rep(2.5, length(n)-2)) 
         tmp = (n-ax)*tail(lx,-1) + ax*head(lx,-1)
         Lx[i,] = c( sum(tmp[1:2]), tmp[3:11])
      } # for i 
    }  else
    if (prior=='empirical') {
        relp = dbeta(HMD_LX[,'q5'], ab[1], ab[2])
        sel = sample(nrow(HMD_LX), nsamp, prob=relp, replace=TRUE)
        Lx  = HMD_LX[sel, paste0('L', seq(0,45,5))] / 1e5
      } else 
    if (prior=='HIV') {  
         Lx = matrix(NA, nsamp, 10, 
                dimnames=list(NULL, paste0('L',seq(0,45,5))))
         this.q5 = rbeta(nsamp,ab[1],ab[2])
         this.prev = runif(nsamp, prev.min, prev.max)
         for (i in seq(this.q5)) {
           mx = mortmod.5q0(child.mort = this.q5[i],
                            prev = this.prev[i],
                            region=0,
                            sex=1)
           x = c(0, 1, seq(5,100,5))
           n = c(diff(x),Inf)
           lx = c(1, exp(-cumsum(n*mx)))
           ax = c( .07+1.5*mx[1], 1.5, rep(2.5, 20)) 
           tmp = (n-ax)*tail(lx,-1) + ax*head(lx,-1)
           Lx[i,] = c( sum(tmp[1:2]), tmp[3:11])
         } # for i
      } # if prior==HIV 
  
  return(Lx)
} # draw_L



#-----------------------------------------------------------------------
# select a random row of PHI and another random row of LX, and calculate
# the Ka(PHI,L) vector

random.kappa = function(fert_prior='text',
                        mort_prior='text',
                        minq5,maxq5,
                        prev.min=0, prev.max=5) {
  
  selected_PHIX = as.vector( draw_phi(nsamp=1, prior=fert_prior) )
  selected_LX   = as.vector( draw_L  (nsamp=1, prior=mort_prior, 
                                      minq5, maxq5,
                                      prev.min,prev.max))
  
  PP = c(0, selected_PHIX)
  names(PP) = paste0('phi', seq(10,45,5))

  LL = selected_LX
  names(LL) = paste0('L', seq(0,45,5))
  
  surv.num  = LL[paste0('L', seq(10,40,5))]
  surv.den  = LL[paste0('L', 5+seq(10,40,5))]
  surv      = surv.num/surv.den
  phi.early = PP[paste0('phi', seq(10,40,5))]
  phi.later = PP[paste0('phi', 5+seq(10,40,5))]
  kappa     = LL['L0']/5 * 1/2 * (surv * phi.early + phi.later)
  names(kappa) = paste0('kappa',seq(15,45,5))
  return(kappa)
}


random.kappa(fert_prior='text', mort_prior='text', minq5=.02, maxq5 = .02,
             prev.min=0, prev.max=5)

#################################################################
# read Kanamari and GA data
#################################################################

tmp = read.csv('test data with all counties.csv', stringsAsFactors = FALSE) %>%
  mutate(LTFR = iTFR / 2,
         HFTR = iTFR * 2) 

G = read.csv(text="
Fprior,Mprior
text,text
empirical,empirical
dirichlet,HIV
", stringsAsFactors=FALSE)

## get the GA county Stan fits  and the Kanamari fit for comparison
library(rstan)
load('Stan GA 2010 counties fit Fri 11Aug17 1410.RData')
GA_fit = fit

load('Kanamari fit.RData')
Kanamari_fit = fit

# alphabet order of GA counties 1=Appling Cnty...159=Worth Cnty,
# then 160=Kanamari do Rio Juruá

##########################################################################
# estimate posteriors for each county + Kanamari under three alternative
# prior combinations
#
# (1) as in text: log-odds fertility model + Wilmoth et al. mortality
# (2) empirical: select empirical HFD Phix and HMD Lx schedules at random (q5-wtd mortal)
# (3) exotic: phi ~ Dirichlet,  Lx ~ HIV schedule with prev of 10% (q5-wtd)
#
# The big loop below draws plots and creates a data frame of quantiles
##########################################################################

df = data.frame()

for (this.county in 1:nrow(tmp) )  {  

  P = list()
  
  for (g in 1:nrow(G)) {

    fert_prior = G$Fprior[g]
    mort_prior = G$Mprior[g]
    
    minq5 = tmp$minq5[this.county]
    maxq5 = tmp$maxq5[this.county]
    
    W = tmp[this.county, paste0('W', seq(15,45,5))]
    C = tmp$C[this.county]
    
    nsim = 1000
    
    logP = matrix(NA, length(TFRvals), nsim, 
                  dimnames=list(TFRvals,NULL))
    
    for (sample_num in 1:nsim) {
      Kstar = sum( W * random.kappa(fert_prior, mort_prior, 
                                    minq5, maxq5,
                                    prev.min=10, prev.max=10))
    
      logP[,sample_num] = -Kstar * TFRvals  + C * log(TFRvals)   
    }
    
# rescale
    logP = sweep(logP, 2, apply(logP,2,max), '-') 
    
    relP = prop.table( exp(logP), margin=2)  # standardized so the each colsum=1
    
#    relP[!is.finite(relP)] = 0
    
    wt = rowSums(relP)
    wt = wt / (delta * sum(wt))
    
    P[[g]] = list(Fprior=fert_prior,
                  Mprior=mort_prior,
                  density = list(x=TFRvals, y=wt))
##--- save quantiles in data frame
    cumDist =  cumsum(delta * P[[g]]$density$y )   
    Q = approx(x =cumDist , y=TFRvals, xout=c(.10,.50,.90))$y
    
    new.info = expand.grid( priors      = paste(fert_prior,mort_prior,sep='+'),
                            county      = this.county,
                            county.name = tmp$county.name[this.county],
                            quantile    = c(10,50,90))
    new.info$value = Q
                            
    df = rbind(df, new.info)
    
  } # for g
  
##--- plots  
  plot(  P[[1]]$density, col='black', lwd=3, type='l',main=tmp$county.name[this.county],
         xlim=c(tmp$LTFR[this.county], tmp$HFTR[this.county]))
  lines( P[[2]]$density, col='red', lwd=3)
  lines( P[[3]]$density, col='dodgerblue', lwd=3)
  
  if (this.county==160) {
       z = extract(Kanamari_fit, 'TFR' )[[1]]
     } else {
       z = extract(GA_fit,paste0('TFR[',this.county,']'))[[1]]
     }
  
  lines(density(z), col='forestgreen', lwd=4)
  
  
  
  
} # for this.county


write.csv(df, file='big test of alternative priors.csv')

##############################################################
#  FACET PLOT OF QUANTILES
##############################################################

# objective: rearrange for plotting:
# one row per (county, quantile, alternative prior)
# with county, quantile, Qtext, Qalt, cat of alternative

base = df %>% 
        group_by(county,quantile) %>% 
        summarize(Q = value[1]) 

head(base)

## for legacy reasons, create a copy of the priors factor with difft levels
df$cat = df$priors
levels(df$cat) = c('text','EMPIRICAL','IDBHIV')

alt1 = df %>%
         filter(cat=='EMPIRICAL') %>%
         select(county, quantile, value,cat) %>%
         rename(Qalt = value)

new1 = left_join(base, alt1, by=c('county','quantile'))
         
alt2 = df %>%
  filter(cat=='IDBHIV') %>%
  select(county, quantile, value,cat) %>%
  rename(Qalt = value)

new2 = left_join(base, alt2, by=c('county','quantile'))

plot.df = rbind( new1, new2) %>%
            arrange(county,quantile,cat) %>%
          mutate(quantile = paste0(quantile,'%ile'),
                 cat = paste('Alternative prior =',cat))

## summary plot of quantiles under alternative priors

pdf(file='big-test-of-alternative-priors.pdf', width=11, height=8)

ggplot(data=filter(plot.df, county < 160), aes(x=Q, y=Qalt)) + 
  geom_point(col='red',shape=16, size=2, fill='red',alpha=.50) + 
  geom_abline(intercept=0, slope=1) + 
  theme_bw() +
  theme(axis.title=element_text(face='bold',size=14),
        strip.text=element_text(face='bold',size=12)) +
  facet_wrap(cat ~ quantile) +
  labs(x='Posterior Quantile: Proposed prior', 
       y='Posterior Quantile: Alternative prior')

dev.off()
