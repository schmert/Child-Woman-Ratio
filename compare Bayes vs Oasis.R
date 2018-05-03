### Figures from the Stan model output for GA county TFRs 

rm(list=ls())
graphics.off()
windows(record=TRUE)

set.seed(490118)

library(rstan)
library(dplyr)

stan_file = 'Stan GA 2010 counties fit Fri 11Aug17 1410.RData'
load(stan_file)

ncounty = stanDataList$n

fit.summary = summary(fit,'TFR', probs=c(10,50,90)/100)$summary

TFR.df$Bayes   = fit.summary[,'50%']
TFR.df$sdBayes = fit.summary[,'sd']
TFR.df$Q10     = fit.summary[,'10%']
TFR.df$Q90     = fit.summary[,'90%']

TFR.df$county.fips = as.character(TFR.df$county.fips)

nchs = read.csv('nchs 2006-2010 tfr.csv', stringsAsFactors = FALSE) %>%
    mutate(county.fips = as.character(fips)) %>%
    rename(TFR_NCHS = TFR)


oasis = read.csv('BirthRates2006-2010.csv', skip=4,na.strings='.',
                 stringsAsFactors = FALSE) %>%
          rename(TFR=TFR1544) %>%
          mutate(county.fips = as.character(county.fips))

df = left_join(oasis, select(TFR.df, -(W15:W45)), by='county.fips')

df = left_join(df, select(nchs, county.fips, TFR_NCHS), by='county.fips')

# jettison 12 NA and create pop size factor
qq = quantile(df$F1549, c(.90), na.rm=TRUE)

df = filter(df, is.finite(TFR)) %>%
       mutate(sizecat = cut(F1549, breaks=c(-Inf,15000,+Inf),
                            labels=c('Small Counties\nLT 15000 Women 15-49',
                                     'Large Counties\nGT 15000 Women 15-49')),
              annotation = ifelse(county.name %in% c('Fulton','Gwinnett','DeKalb','Cobb'), 
                                  county.name,''))


G = ggplot(df, aes(x=TFR, y=Bayes, label=annotation,size=F1549)) +
      scale_x_continuous(breaks=seq(1.2, 3.0, 0.4)) +
      scale_y_continuous(breaks=seq(1.2, 3.0, 0.4)) + 
      labs( x='Vital Statistics TFR (2006-2010 avg)',
            y ='Bayes Posterior Median and 80% interval',
            size='Female Pop 15-49') +
      geom_abline(intercept=0,slope=1) +
      geom_segment(aes(x=TFR, xend=TFR, y=Q10, yend=Q90), 
                   color='navy', lwd=0.6, alpha=.50) +
      geom_point(color='navy', fill='navy',alpha=.50) +
      geom_text(size=4, aes(x=TFR+.05), hjust=0) +
      theme_bw() +
      theme(axis.title = element_text(face='bold',size=12),
            axis.text = element_text(face='bold',size=11))


G = G + 
     facet_wrap(~sizecat) + 
     theme(strip.text=element_text(face='bold',size=12))


cairo_pdf(file = 'compare-Bayes-vs-Oasis.pdf', width=12, height=7)
   print(G)
 dev.off()
 
cairo_ps(file = 'compare-Bayes-vs-Oasis.eps', width=12, height=7)
   print(G)
 dev.off()  

