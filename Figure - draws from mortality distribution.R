#######################################################################
# illustrate draws from the (q5, k) distribution
#######################################################################

graphics.off()
if (.Platform$OS.type == 'windows') windows(record=TRUE)
rm(list=ls() )

set.seed(6448)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)


  ## This is the learnBayes::beta.select function by Jim Albert
  ##  each argument is a list with $p and $x components for quantile prob. and variable value
  beta.params = function (quantile1, quantile2) 
  {
    betaprior1 = function(K, x, p) {
      m.lo = 0
      m.hi = 1
      flag = 0
      while (flag == 0) {
        m0 = (m.lo + m.hi)/2
        p0 = pbeta(x, K * m0, K * (1 - m0))
        if (p0 < p) 
          m.hi = m0
        else m.lo = m0
        if (abs(p0 - p) < 1e-04) 
          flag = 1
      }
      return(m0)
    }
    p1 = quantile1$p
    x1 = quantile1$x
    p2 = quantile2$p
    x2 = quantile2$x
    logK = seq(-3, 8, length = 100)
    K = exp(logK)
    m = sapply(K, betaprior1, x1, p1)
    prob2 = pbeta(x2, K * m, K * (1 - m))
    ind = ((prob2 > 0) & (prob2 < 1))
    app = approx(prob2[ind], logK[ind], p2)
    K0 = exp(app$y)
    m0 = betaprior1(K0, x1, p1)
    return(round(K0 * c(m0, (1 - m0)), 2))
  }
  
  ## identify a and b such that when q5 ~ beta(a,b) there is a 90% probability that q5 is
  ##  between one-half and twice the estimated q5 for this population

  q5hat = .025
  
  q5.ab = beta.params( list(x=q5hat/2,p=.05), list(x=2*q5hat,p=.95))
  
  
  #------------- MORTALITY model ------------------
  
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
             110,-0.2212,-0.0028,-0.0027,     0,-0.2966,-0.0229,-0.0041,     0')
  
  
    nsim = 30
    sim.q5 = rbeta(nsim, q5.ab[1], q5.ab[2] )
    
    sim.k = rnorm(nsim, mean=0, sd=1)
    
    sim = data.frame( q5=c(q5hat,sim.q5), k = c(0,sim.k))
    
    x = c(0,1, seq(5,110,5))
    
    for (v in c('af','bf','cf','vf')) assign(v, wilmoth[[v]])  # local variables af..vf
    
    logmx = apply(sim, 1, function(p) {
                  h = log(p[1])
                  k = p[2]
                  return( af + bf*h + cf*h^2 + vf*k)
                  })
    
    rownames(logmx) =x 
    
    logmx['1',] = log( -1/4 * (exp(logmx['0',]) + log(1-sim$q5)))
    
    
    matplot(x, logmx, type='o', lty=1)
    
    ## cut off at age 50
    logmx = logmx[x<50,]
    n = c(1,4,rep(5,9))
    
    mx = exp(logmx)
    Hx = apply(mx*n, 2, cumsum)
    
    lx = apply(Hx,2, function(H) exp(-c(0,H)))
    rownames(lx) = c(0,1,seq(5,50,5))
    
    bigL = function(lx) {
      c( 1/2*(lx[1] + lx[2]) + 4/2*(lx[2]+lx[3]),
         5/2*(lx[3:11] + lx[4:12])
         )
    }
    
    Lx = apply(lx,2,bigL)
    
    matplot(seq(0,45,5), Lx, type='o', pch=1, ylim=range(Lx,5),lty=1,
              xlab='a', ylab='La')

    
    qx = 1-lx
    x  = c(0,1,seq(5,50,5))
    rownames(qx) = x
    
    matplot(x, qx, type='o', lty=1)
    
    matplot(x, qx, type='o', lty=1, xlim=c(0,5.5), ylim=range(0,qx[x<10,]),pch=1)
    points( 5, q5hat, pch=16, cex=2)
    
    
## a two-panel ggplot: sampled L0...L45 values and sampled q(x) from (h,k) prior
    
    # make a 'long' data frame for plotting
    sim.df = data.frame( simnum  = as.vector(col(Lx)), 
                           a      = seq(0,45,5)[row(Lx)], 
                           La    = as.vector(Lx) )
    
    G1 = ggplot( data=sim.df, aes(x=a, y=La, group=simnum, color=as.factor(simnum))) +
      geom_line(lwd=0.60) + geom_point(size=2,shape=1) +
      scale_x_continuous('a',limits=c(0,50)) +
      scale_y_continuous(expression(L[a]),limits=range(sim.df$La,5), expand=c(0,0.005)) + 
      theme_bw() +
      theme(legend.position='none', 
            text=element_text(size=15, face='bold') ) +
      geom_line(data=filter(sim.df,simnum==1), color='black', lwd=2, alpha=.70) +
      ggtitle('Life Table Person-Years')
    
    
    print(G1)
    
    
    # make a 'long' data frame for plotting
    sim.df = data.frame( simnum  = as.vector(col(qx)), 
                         x      = c(0,1,seq(5,50,5))[row(qx)], 
                         qx    = as.vector(qx) ) %>%
             filter(x<15)
    
    
    G2 = ggplot( data=sim.df, aes(x=x, y=qx, group=simnum, color=as.factor(simnum))) +
      geom_line(lwd=0.60) + geom_point(size=2,shape=1) +
      scale_x_continuous('x',limits=c(0,10),breaks=c(0,1,5,10)) +
      scale_y_continuous('q(x)',limits=range(0,sim.df$qx), expand=c(0,0.001)) + 
      theme_bw() +
      theme(legend.position='none', 
            text=element_text(size=15, face='bold') ) +
      geom_line(data=filter(sim.df,simnum==1), color='black', lwd=2, alpha=.70) +
      geom_point( x=5, y=q5hat, shape=16, size=6, col='black') +
      ggtitle('Childhood q(x)')
    
    print(G2)
    
    
G = grid.arrange(G1, G2, ncol=2)

ggsave( filename='Figure - draws from mortality distribution.eps', 
        plot=G,device=cairo_ps,
        width=10, height=5, path='./')

    