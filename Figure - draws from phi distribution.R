#######################################################################
# make a plot with typical age patterns of fertility by simulating
# from gamma distribution
#######################################################################

graphics.off()
if (.Platform$OS.type == 'windows') windows(record=TRUE)
rm(list=ls() )

library(ggplot2)
library(dplyr)

IDB5 = read.csv('http://schmert.net/calibrated-spline/REBEP/IDB%205fx.csv', 
                header=TRUE, stringsAsFactors = FALSE, row.names=1)

asfr = as.matrix( IDB5[,-1] )
asfr[asfr==0] = 1e-2  #  per 1000... so 1e-5 per woman

# Nx7 matrix: proportions of fertility in age groups 15..45 
P     = sweep(asfr, 1, rowSums(asfr), '/') 

# Nx7 matrix: proportions in log-ratio form:
#   for each schedule/row, the 6 coefs are log( p[2]/p[1] ... log(p[7]/p[1]))
gamma  = t( apply(P, 1, function(x) log(x/x[1]) ) )[,paste0('f',seq(20,45,5))]  
mgamma = colMeans(gamma)
Vgamma = var(gamma)

small.table = rbind( round(mgamma,2), round(sqrt(diag(Vgamma)),2))
rownames(small.table) = c('mean','sd')

list('gamma distrib'=small.table,
    'gamma corr matrix'=round( cor(gamma), 2))

CH = t(chol(Vgamma))  # (CH)(CH)' = Vgamma

nsim = 25

zztarget = qchisq(p=.50, df=6) #median of a chisq(6) distrib

z = matrix(rnorm(length(mgamma)*nsim), ncol=nsim)
z = sweep( z, 2, sqrt(colSums(z^2) / zztarget), '/')


sim.gamma = mgamma + CH %*% z

## add a first column that corresponds to the mean gamma vector
sim.gamma = cbind( mgamma, sim.gamma)

sim.phi = prop.table( rbind(0,1,exp(sim.gamma)), margin=2) 
rownames(sim.phi) = paste0('F',seq(10,45,5))

ok = apply(sim.phi, 2, function(x) { ((x[2] < .25) & (x[8] < .05)) })
sim.phi = sim.phi[,ok]

# make a 'long' data frame for plotting
sim.data = data.frame( simnum  = as.vector(col(sim.phi)), 
                        a      = seq(10,45,5)[row(sim.phi)], 
                        phi    = as.vector(sim.phi) )
head(sim.data)

G = ggplot( data=sim.data, aes(x=a+2.5, y=phi, group=simnum, color=as.factor(simnum))) +
      geom_line(lwd=0.60) + geom_point(size=2,shape=1) +
      scale_x_continuous('Age',limits=c(10,50), expand=c(0,0.1)) +
      scale_y_continuous('Phi',limits=c(0,.40), expand=c(0,0.005)) + 
      theme_bw() +
      theme(legend.position='none', 
            text=element_text(size=15, face='bold') ) +
      geom_line(data=filter(sim.data,simnum==1), color='black', lwd=2, alpha=.70)

print(G)

ggsave( filename='Figure 2 - draws from phi distribution.eps', device=cairo_ps,
         width=8, height=5, path='./')








