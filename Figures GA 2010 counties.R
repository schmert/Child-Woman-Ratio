### Figures from the Stan model output for GA county TFRs 

rm(list=ls())
graphics.off()
windows(record=TRUE)

library(rstan)
library(ggplot2)
library(ggmap)
library(dplyr)
library(maptools)
library(rgdal)
library(rgeos)
library(grid)


stan_file = 'Stan GA 2010 counties fit Fri 11Aug17 1410.RData'
load(stan_file)

ncounty = stanDataList$n

fit.summary = summary(fit,'TFR', probs=c(10,50,90)/100)$summary

TFR.df$Bayes   = fit.summary[,'50%']
TFR.df$sdBayes = fit.summary[,'sd']
TFR.df$Q10     = fit.summary[,'10%']
TFR.df$Q90     = fit.summary[,'90%']

TFR.df$county.fips = as.character(TFR.df$county.fips)

#######################

stan_trace( fit,  c('mu_TFR','sigma_TFR'))
stan_trace( fit,  paste0('TFR[', sort(sample(ncounty,10)), ']')) + ylim(1,4) + ggtitle('10 Random Areas')
stan_dens( fit,  paste0('TFR[', sort(sample(ncounty,10)), ']'), separate_chains = TRUE) + xlim(1,4) + ggtitle('10 Random Areas')

R = summary(fit)$summary[,'Rhat']
stan_trace(fit, names(tail(sort(R),20))) + ylim(0,4) + ggtitle('Highest 20 Rhat values')

RTFR = R[paste0('TFR[', 1:ncounty, ']')]
stan_trace(fit, names(tail(sort(RTFR),20))) + ylim(0,4) + ggtitle('Highest 20 Rhat TFR values')
stan_dens(fit, names(tail(sort(RTFR),20)), separate_chains = TRUE) + xlim(0,4) + ggtitle('Highest 20 Rhat TFR values')


ggplot(data=TFR.df, aes(x=iTFR, y=Bayes, size=F1549)) + 
  geom_point(shape=18) + 
  geom_abline(intercept=0, slope=1, col='red') + 
  theme_bw() +
  ggtitle('iTFR versus Bayesian posterior means: GA counties')

ggsave(file='GA counties iTFR vs Bayes.eps', device=cairo_ps)
ggsave(file='GA counties iTFR vs Bayes.pdf', device=cairo_pdf)


##############################################################################
# map the results
##############################################################################


## read shape file

needs.download = !file.exists('gz_2010_us_050_00_500k.shp')

if (needs.download) {
  download.file(url = "http://www2.census.gov/geo/tiger/GENZ2010/gz_2010_us_050_00_500k.zip",
                destfile = "./gz_2010_us_050_00_500k.zip")
  
  unzip("gz_2010_us_050_00_500k.zip")
}

## convert to a spatial polygons data frame 

USA  = readOGR(dsn='.', layer='gz_2010_us_050_00_500k')  # this keeps projection info
GA  = USA[USA$STATE==13,]

GA@data = GA@data %>%
  mutate( county.fips = paste0(as.character(STATE), as.character(COUNTY) ) )

## translate the polygons into a standard dataframe, for use with ggmap
GA.df = fortify(GA, region='county.fips') %>%
  left_join( select(TFR.df, county.fips, iTFR, Bayes, sdBayes), by=c('id'='county.fips'))

tmp = data.frame( coordinates(GA), county.fips = GA$county.fips) %>%
       left_join( select(TFR.df, county.name, county.fips, Bayes), by='county.fips') %>%
       rename(long=X1, lat=X2) %>%
       mutate(Bayes = round(Bayes,1))
                  
G = ggplot(data=GA.df, aes(x=long, y=lat, group=id, fill=Bayes)) + 
  geom_polygon(color='lightgrey') + 
  scale_fill_gradient('TFR',low='white', high='dodgerblue') + 
  theme_nothing(legend=TRUE) +
  theme(title=element_text(size=14, face='bold', hjust=0.5)) +
  #  labs(title='Posterior Median TFR\nGA Counties 2010') +
  geom_text(data=tmp, aes(x=long, y=lat, label=Bayes), 
            color='black', inherit.aes = FALSE, size=3,
            check_overlap = TRUE) +
#  coord_equal() +
  geom_segment(x=-84.5, y=30.4, xend=-83.46, yend=30.4, lwd=1.5) +
  geom_text(x=-84, y=30.5, label='100 km')

print(G)

## add metro ATL outline

ATL = GA[ GA$NAME %in% c('Cobb', 'Fulton', 'DeKalb'),]
ATL = gUnaryUnion(ATL)
ATL.df = fortify(ATL)

UNIV = GA[ GA$NAME %in% c('Clarke','Bulloch'),]
UNIV.df = fortify(UNIV)

JAIL = GA[ GA$NAME %in% c('Pulaski'),]
JAIL.df = fortify(JAIL)

BASE = GA[ GA$NAME %in% c('Chattahoochee'),]
BASE.df = fortify(BASE)

highlight.df = rbind( ATL.df, UNIV.df, JAIL.df, BASE.df)

G + geom_polygon(data=highlight.df, aes(x=long,y=lat, group=id), 
                 inherit.aes = FALSE, fill=NA, color='black',
                 lwd=1) +
    annotate("text", y=34.9, x=-82.65, label='Metro Atlanta', size=5, adj=0) +
    annotate("segment", y=34.9, x=-82.7, yend=33.95, xend=-84.39, lwd=1.2,
             alpha=.30) +

    annotate("text", y=34.2, x=-82.3, label='Universities', size=5,adj=0) +
    annotate("segment", y=34, x=-83.2,yend=34.2, xend=-82.4, lwd=1.2,
           alpha=.30) +
    annotate("segment", y=32.6, x=-81.75,yend=34.05, xend=-81.75, lwd=1.2,
         alpha=.30) +
  
    annotate("text", y=31.5, x=-85.9, label='Womens\nPrison', size=5,adj=0) +
   annotate("segment", y=32.25, x=-83.6,yend=31.3, xend=-85.3, lwd=1.2,
         alpha=.30) +
  
   annotate("text", y=32.5, x=-85.9, label='Army\nBase', size=5,adj=0) +
   annotate("segment", y=32.3, x=-85.4,yend=32.3, xend=-85.05, lwd=1.2,
         alpha=.30) 

ggsave(file='GA County TFR map.eps', device=cairo_ps)
ggsave(file='GA County TFR map.pdf', device=cairo_pdf)


Q        = data.frame( summary(fit,'TFR', probs=c(.10,.50,.90))[[1]][ , c('10%', '50%', '90%')])
Q$iTFR   = TFR.df$iTFR
names(Q) = c('Q10', 'Q50', 'Q90', 'iTFR')
Q$index  = rank(Q$Q50)

ggplot(data=Q, aes(x=Q50, y=index)) + 
  xlim(1,3.5) + xlab('County TFR') + ylab('') +
  geom_segment(aes(x=Q10, xend=Q90, y=index, yend=index), col='grey') +
  geom_point(size=1) +
  geom_point(aes(x=iTFR, y=index), shape='x',size=3, col='red') +
  theme_bw() 

ggsave(file='GA counties TFR - 80pct intervals.eps', device=cairo_ps,
       width=10, height=7)
ggsave(file='GA counties TFR - 80pct intervals.pdf', device=cairo_pdf,
       width=10, height=7)

get_elapsed_time(fit)

#################

pdf(file='Stan GA 2010 counties misc plots.pdf') 

    TFR.df = TFR.df %>%
      mutate(pct = (Bayes-iTFR)/iTFR,  
             P   = (W25+W30)/F1549)
    
    # % diff Bayes - iTFR
    ggplot(data=TFR.df,  aes(pct)) + 
      geom_density(adjust=1.5) +
      geom_vline(xintercept = 0) +
      theme_bw()
    
    # fraction of women 25-34 vs (Bayes-TFR)
    ggplot(data=TFR.df,  aes(x=P, y=pct)) + 
      xlab('Fraction of Women 25-34') + ylab('% diff Bayes - iTFR') +
      geom_point() +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 2/7) +
      theme_bw()
    
    a = seq(15,45,5)
    
    # + adjustments: Bayes > iTFR
    isel = tail(order(TFR.df$pct), 5)   # indices for 5 biggest upward adjustments iTFR -> Bayes
    for (i in isel) {
      this.label = paste('Positive Adj:', TFR.df$county.name[i],'County', '\niTFR=',round(TFR.df$iTFR[i],2),'Bayes=',round(TFR.df$Bayes[i],2))
      W = TFR.df[i, paste0('W',a)]
      plot(a+2.5, W , type='o', ylim=range(0,W), main=this.label)
      points(c(25,30)+2.5, W[c('W25','W30')] , col='blue', pch=16, cex=1.5)
      abline(h=sum(W)/7, lty=2)
    }
    
    # - adjustments: Bayes < iTFR
    isel = head(order(TFR.df$pct), 5)   # indices for 5 biggest upward adjustments iTFR -> Bayes
    for (i in isel) {
      this.label = paste('Negative Adj:', TFR.df$county.name[i],'County', '\niTFR=',round(TFR.df$iTFR[i],2),'Bayes=',round(TFR.df$Bayes[i],2))
      W = TFR.df[i, paste0('W',a)]
      plot(a+2.5, W , type='o', ylim=range(0,W), main=this.label)
      points(c(25,30)+2.5, W[c('W25','W30')] , col='red', pch=16, cex=1.5)
      abline(h=sum(W)/7, lty=2)
    }
    
    # almost no adjustments: Bayes = iTFR
    isel = which( abs(TFR.df$pct) < .01)   # indices
    for (i in isel) {
      this.label = paste('Almost Zero Adj:', TFR.df$county.name[i],'County', '\niTFR=',round(TFR.df$iTFR[i],2),'Bayes=',round(TFR.df$Bayes[i],2))
      W = TFR.df[i, paste0('W',a)]
      plot(a+2.5, W , type='o', ylim=range(0,W), main=this.label)
      points(c(25,30)+2.5, W[c('W25','W30')] , col='black', pch=16, cex=1.5)
      abline(h=sum(W)/7, lty=2)
    }


dev.off()


######### Bayes - iTFR adjustments as a function of %Women 25-34

county.names = TFR.df$county.name

W = as.matrix( select(TFR.df, contains('W')))
dimnames(W) = list( county.names, seq(15,45,5))

delta = TFR.df$Bayes - TFR.df$iTFR 

P2534 = rowSums(W[,c('25','30')]) / rowSums(W)

tmp = cbind( select(TFR.df, county.name, county.fips, iTFR, Bayes), 
             delta=delta,
             P2534=P2534)

G = ggplot(data=tmp, aes(x=P2534, y= 7*Bayes/iTFR)) +
  xlab('Women 25-34 / Women 15-49') + ylab('Bayes Median / Child-Woman Ratio') +
  geom_point(size=2.5) +  
  geom_hline(yintercept = 7, lty=2) +
  geom_vline(xintercept = 2/7, lty=2) +
  geom_text(x=2/7, y= 5.75, label='2/7' ,size=5) +
  theme_bw() +
  geom_text(x=.155, y=9.5, label='Fayette', adj=-1, size=5) +
  geom_text(x=.375, y=5.75, label='Chattahoochee', adj=1, size=5) +
  scale_x_continuous(breaks=seq(.10,.50,.05), minor_breaks = NULL) +
  scale_y_continuous(breaks=0:10, minor_breaks = NULL) 

print(G)

##
miniplot = function(df) {
  ggplot(data=df,  aes(x=a, y=W)) + ylim(0,max(W)) + xlim(10,48) +
    geom_line(lwd=1, color='darkgrey') +
    geom_point(size=1, color='black') + 
    scale_x_continuous(breaks=c(20,30,40), labels=c('20','30','40')) +
    theme_bw() +
    theme(axis.line = element_line(colour = 'darkgrey', size=0.5),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size=6, face='bold'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
}
##

W = TFR.df %>% 
  filter(county.name=='Fayette') %>%
  select(starts_with('W'))

mini.df = data.frame(a=seq(15,45,5), W=as.numeric(W))

mini_Fayette = miniplot(mini.df)

W = TFR.df %>% 
  filter(county.name=='Chattahoochee') %>%
  select(starts_with('W'))

mini.df = data.frame(a=seq(15,45,5), W=as.numeric(W))

mini_Chattahoochee = miniplot(mini.df)

vp_C <- viewport(width = 0.1, height = 0.1, x = .69,y=.116)
vp_F <- viewport(width = 0.1, height = 0.1, x = .28,y=.93)

print(G)
print(mini_Chattahoochee, vp=vp_C)
print(mini_Fayette, vp=vp_F)

cairo_pdf(file='GA counties Bayes iTFR difference.pdf', height=8, width=8)
  print(G)
  print(mini_Chattahoochee, vp=vp_C)
  print(mini_Fayette, vp=vp_F)
dev.off()

cairo_ps(file='GA counties Bayes iTFR difference.eps', height=8, width=8)
  print(G)
  print(mini_Chattahoochee, vp=vp_C)
  print(mini_Fayette, vp=vp_F)
dev.off()



############# compare to NCHS estimates for big counties ##############

nchs = read.csv('NCHS fertility data.csv', stringsAsFactors = FALSE) %>%
        mutate(county.fips = as.character(county.fips))

df = left_join(nchs, select(TFR.df, -(W15:W45)), by='county.fips')

ggplot(df, aes(x=TFR, y=Bayes, label=county.name)) +
    geom_point(color='orangered', alpha=.50, size=4) +
    xlim(1.3, 2.4) + ylim(1.3, 2.4) +
    xlab('NCHS (gold standard) TFR') + ylab('Bayes Posterior Median and 80% interval') +
    geom_abline(intercept=0,slope=1) +
    geom_segment(aes(x=TFR, xend=TFR, y=Q10, yend=Q90), 
                 color='orangered', lwd=1.5, alpha=.50) +
    geom_text(size=3) + 
    theme_bw()

ggsave(file='NCHS vs Bayes intervals.eps', device = cairo_ps)
ggsave(file='NCHS vs Bayes intervals.pdf', device = cairo_pdf)


###### examine the implications of the NCHS level for the number of
###### children expected in a selected county

this.county = 'Fulton'

for (this.county in df$county.name) {
  
i.fit        = which(TFR.df$county.name == this.county) 
i.nchs       = which(df$county.name == this.county) 

## get the q5 draws for this county (prior and posterior should be extremely similar)
q5_sim = as.matrix(fit, paste0('q5[', i.fit, ']'))

## a crude and conservative estimate for L5 in all these simulations:
##  5*(1-q5) < L5 < 5.00,  so   (1-q5) < s0 < 1,  and we'll just use the
##  minimum possible 0-4 s0 value (=> max possible 1/s0 ==> max possible 
##  mortality effect observed on C and C/W)

s0_sim = 1-q5_sim

## simulate fertility age patterns from the prior

g_sim = stanDataList$m[1,] + 
         t(stanDataList$X) %*% matrix(rnorm(2*length(q5_sim)), nrow=2)

phi_sim = apply(g_sim, 2, function(x) exp(x) / sum(exp(x)))

tmp = rbind(0, phi_sim, 0)
rownames(tmp) = seq(10,50,5)

pa_sim = apply(tmp, 2, function(x) 1/2 * (head(x,-1) + tail(x,-1)))[-8,]
rownames(pa_sim) = seq(15,45,5)

gold_std_TFR = nchs[i.nchs, 'TFR']

W = as.numeric( select(TFR.df[i.fit,], starts_with('W')))

gold_std_expected_C = gold_std_TFR * 1/s0_sim * as.vector(W %*% pa_sim)

C_sim = data.frame(C=rpois( length(gold_std_expected_C), gold_std_expected_C))

C_obs = TFR.df[i.fit, 'C']
  
G = ggplot(data=C_sim, aes(x=C)) +
     geom_density(adjust=1.5, fill='dodgerblue',color='blue', alpha=.50) +
     xlab('Number of 0-4 year olds at NCHS level of TFR') +
     geom_vline(xintercept = C_obs, lwd=2, lty=2) +
     ggtitle(paste(this.county,'County, NCHS TFR=',
                   round(gold_std_TFR,2))) +
     theme_bw() 

print(G)

}

