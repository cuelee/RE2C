##########################
##
##    RE2C 
##    New RE2 more powerful than RE2
##   
#####  Buhm Han & Cue Lee

### 2/18/15: v1: launch
### 2/29/15: v2: Speed up using eigen decomp. Separated calculating low/high thresholds. (Assume we always input these)
### 2/20/15: v3: I try to make p-value continuous. (not stacked at 1) 
### 9/4/15: v4: with null test function.. I will try empirical distribution one more time.
### 9/8/15: v5: fast implementation of using empirical table.
### 8/20/16: v6: set Shet.high at Stat1
### 11/18/16: v7: account correlation of shared controls
### 12/9/16: v8: update the default tlow values under GWASthreshold of 5*10^-8
### 1/11/17: v9: change the p-value estimation to sum_{x in SFE bin} ( Prob(Shet > max(S-x, Thres(x)))*Prob(x) )
### 3/20/17: v10: update to check dependency of mvtnorm package before run RE2C

### RE2C calculate the parameter tlow; To estimate tlow, we need 1. Han and Eskin pvalue table, 2.number of studies, 3.target GWAS threshold (default 5*10^-8)
## I recommend to calculate tlow for every available study number which is ranged from 2 to 50 (later 100) and run RE2h
## provided default tlow with conventional GWASthreshold of  5*10^-8

## RE2C uses R. package names <mtvnorm>
suppressWarnings(pkgTest("mvtnorm"))

## RE2h main
RE2C <- function(beta, stders, sigma, stat2_cdf = FALSE, FEt.lows = FALSE, NR, tau2.zero.prob) {
   # n study numbers
   # stopifnot (length(beta) == length(stders))
   n=length(beta)
   vars <- stders ** 2
   ws <- 1/vars 
   sigmainv <- solve(sigma)
   beta_rv <- matrix(beta,nrow=1,ncol=n)
   ones <- matrix(rep(1,n),nrow=1,ncol=n)
   sumW <- sum(ws)
   sumW2 <- sum(ws ** 2)
   meanBeta <- (ones %*% sigmainv %*% t(beta_rv)) / (ones %*% sigmainv %*% t(ones))
   Q <- (beta_rv - rep(meanBeta,n))** 2 %*% ws
   meanW <- mean(ws)
   Sw2 <- 1/(n-1) * (sumW2 - n*meanW*meanW)
   U = (n-1)*(meanW - Sw2/(n*meanW))
   tau2 <- max( (Q-(n-1))/U, 0 )
   
   ##-----------------------------------------------
   ## Eigen-decomposition optimization (EMMA- style)
   ##-----------------------------------------------
   K <- sigma
   eig.L <- eigen(K,symmetric=TRUE)
   L.values <- eig.L$values
   L.vectors <- eig.L$vectors
   S <- diag(n)-matrix(1/n,n,n)
   eig.R <- eigen (S%*%K%*%S,symmetric=TRUE)
   R.values <- eig.R$values[1:(n-1)]
   R.vectors <- eig.R$vectors[,1:(n-1)]
   etas <- crossprod(R.vectors,t(beta_rv))
   etasqs <- etas^2
   xis <- L.values
   lambdas <- R.values
   
   mle.tau2 <- NR(tau2,n,xis,etasqs,lambdas)
   Hinv <- solve(sigma+mle.tau2*diag(n))
   mle.beta <- (ones %*% Hinv %*% t(beta_rv)) / (ones %*% Hinv %*% t(ones))
   mle.ll <- -LL.fun(mle.tau2,n,xis,etasqs,lambdas)
   
   null.ll = -ll.fun(c(0,0),beta_rv,sigma)
   lrt.mle <- -2*(null.ll-mle.ll)
   
   stat1 = -2*(null.ll+ll.fun(c(meanBeta,0),beta_rv,sigma))
   stat2 = -2*(-ll.fun(c(meanBeta,0),beta_rv,sigma)+ll.fun(c(mle.beta,mle.tau2),beta_rv,sigma))
   
   p.RE2_cond <- tau2.zero.prob[n-1]*pchisq(lrt.mle,1,lower.tail=F)+(1-tau2.zero.prob[n-1])*pchisq(lrt.mle,2,lower.tail=F)
   
   ##----------------------------------------------
   ## Given two stats, let's calculate p-value
   ##----------------------------------------------
   stat=stat1+stat2
   approx_stat1=min(stat1,49.974)
   if(c(FEt.lows)[floor(approx_stat1*20)+1]>stat2) {
     return(
       list(
         c(stat1,stat2,p.RE2_cond,1)
         ,c(Q)
       )
     )
     }
   if(stat>50){
      p.RE2C <- RE2C_ext(stat,stat2_cdf=stat2_cdf,tmax=FEt.lows[1000],n,tau2.zero.prob=tau2.zero.prob)
      return(
        list(
          c(stat1,stat2,p.RE2_cond,p.RE2C)
          ,c(Q)
        )
      )
   }
   
   p.RE2C <- ((1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*(apply(cbind(FEt.lows,(rep(stat,1000)-FEs)),1,max)))+1])%*%FEprobs + (1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*FEt.lows[1000]+1)]*pchisq(50, 1, lower.tail=F)
   return(
     list(
       c(stat1,stat2,p.RE2_cond,p.RE2C)
       ,c(Q)
     )
   )
}

RE2C_ext <- function(stat,stat2_cdf,tmax,n,tau2.zero.prob){
   extra<-stat-50
   modFEs<-seq(extra+0.25,stat,0.25)
   modFEprobs<-cal_FEprobs(modFEs)
   p.RE2C=(((1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*(apply(cbind(rep(tmax,200),(rep(50,200)-FEs_ext)),1,max)))+1])%*%modFEprobs + (1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*tmax+1)]*pchisq(stat, 1, lower.tail=F))
   return(p.RE2C)
}
