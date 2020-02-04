

library(splines)
library(genlasso)
library(wavethresh)

root <- "C:/Users/rbrow/OneDrive/Documents/Public Github/Sparse-Spline-Penalties/"
source(paste0(root,"spline_fitting_functions.R"))


# x <- seq(0,5,length.out=500)
# y <- sapply(x,bathtub1)
# plot(x,y,type="l")
# 
# x <- seq(0,5,length.out=500)
# y <- sapply(x,doppler1)
# plot(x,y,type="l")



oneiter <- function(truefunc,x,nknots)
{
  ##########################################################
  # Run one iteration of the simulation for for a particular outcome-knot combination
  ########################################################## 
  
  # truefunc <- dp1
  # x <- x
  # nknots <- 50
  
  y <- adderror(truefunc)
  fit.r <- tryCatch( fit.spline(y=y, x=x,nknots=nknots,degree=3,basis="tpoly",smooth.type="ridge.mlm"),
                    error=function(e) { return(list(fitted=rep(NA,150)))} )
  fit.l <- fit.spline(y=y, x=x,nknots=nknots,degree=3,basis="tpoly",smooth.type="lasso")
  fit.f <- fit.spline(y=y, x=x,nknots=nknots,degree=3,basis="tpoly",smooth.type="fusion")
  
  fitted.matrix <- rbind(fit.r$fitted,fit.l$fitted,fit.f$fitted)
  dimnames(fitted.matrix) <- list( method=c("mlm","lasso","fusion"),x=NULL)
  
  return(fitted.matrix)
}


onesetting <- function(truefunc,x,nknots,Nsim) # Nsim iterations for a particular outcome-knot combo
{ 
  out <- replicate(Nsim,oneiter(truefunc,x,nknots),simplify="array")  
  names(dimnames(out))[3] <- "simiter"
  dimnames(out)[[3]] <- 1:Nsim
  return(out)
}


########################################


Nsim <- 5
n <- 150
x <- seq(0,5,length.out=n)
dp1 <- sapply(x,doppler1)
bt1 <- sapply(x,bathtub1)
nknots <- c(15,50,NA)

####### timing tests ##############
# test <- onesetting(dp1,x,50,3)
# 
# system.time( onesetting(dp1,x,NA,5)  )
# system.time( onesetting(dp1,x,50,5)  )
# system.time( onesetting(dp1,x,15,5)  )
# 
# system.time( onesetting(bt1,x,NA,5)  )
# system.time( onesetting(bt1,x,50,5)  )
# system.time( onesetting(bt1,x,15,5)  )


## for 15 knots 1000 sims takes ~50 minutes
## for 50 knots 1000 sims takes ~65 minutes
## for knot @ each x, 1000 sims takes ~110 minutes 

## so, we'll split into 4 separate processes running in parallel
##   1. 15 & 50 knots, bathtub
##   2. 15 & 50 knots, doppler
##   3. knot @ each x, bathtub
##   4. knot @ each x, doppler
## THEN, simulation should be done in ~ 2 hours for 1000 sims

########################################

Nsim <- 1000


# first core
set.seed(515)
bt_eachx <- onesetting(bt1,x,NA,Nsim)
save(bt_eachx,file=paste0(root,"Simulation Results/bt_eachx.rdata"))


# second core
set.seed(515)
dp_eachx <- onesetting(dp1,x,NA,Nsim)
save(dp_eachx,file=paste0(root,"Simulation Results/dp_eachx.rdata"))


# third core
set.seed(515)
bt_15kn <- onesetting(bt1,x,15,Nsim)
save(bt_15kn,file=paste0(root,"Simulation Results/bt_15kn.rdata"))

set.seed(515)
bt_50kn <- onesetting(bt1,x,50,Nsim)
save(bt_50kn,file=paste0(root,"Simulation Results/bt_50kn.rdata"))


# fourth core
set.seed(515)
dp_15kn <- onesetting(dp1,x,15,Nsim)
save(dp_15kn,file=paste0(root,"Simulation Results/dp_15kn.rdata"))

set.seed(515)
dp_50kn <- onesetting(dp1,x,50,Nsim)
save(dp_50kn,file=paste0(root,"Simulation Results/dp_50kn.rdata"))


