# print("default_use : LS(beta,stders,cor)")

## Fixed effect method (weighted mean average)
LS_R <- function(beta, stders, cor, ...) {
   ## conventional FE approach
   
   x <- beta
   V <- diag(stders) %*% cor %*% diag(stders)
   Vinv <- solve(V)
   ones <- matrix(rep(1,length(beta)),nrow=1)
   
   newx <- (ones %*% Vinv %*% x) / (ones %*% Vinv %*% t(ones))
   newv <- 1 / (ones %*% Vinv %*% t(ones))
   newstd <- sqrt(newv)
   
   newp <- pchisq(newx*newx/newv, 1, lower.tail=F)
   return(c(round(newx,2),round(newv,2),newp))
}
LS <- cmpfun(LS_R ,options=list(optimize=1,suppressAll=T,suppressUndefined=T))