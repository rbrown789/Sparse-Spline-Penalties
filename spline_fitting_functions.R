
# library(splines)
# library(genlasso)
# library(wavethresh)
# 
# # root <- "C:/Users/rbrow_000/Google Drive/Rich Par Models/Final Project/"
# root <- "C:/Users/roland/Google Drive/Rich Par Models/Final Project/"


# function for linear interpolation of y values given target x
# and surrounding y and x values
lin.int.y <- function(targx,xlo,xhi,ylo,yhi)
{
  slope <- (yhi-ylo)/(xhi-xlo)
  yval <- slope*(targx-xlo)+ylo
  return(yval)
}


# function for linear interpolation of x value given target y
# and surrounding y and x values
lin.int.x <- function(targy,xlo,xhi,ylo,yhi)
{
  slope <- (yhi-ylo)/(xhi-xlo)
  xval <- (1/slope)*(targy-ylo)+xlo
  return(xval)
}


basisgen <- function(xvals,nknots,degree,type="tpoly")
{
  #####################################################################
  # Function to generate the spline basis.  Returns a list with the design matrix,
  # number of knots, knot locations, degree of the polynomial basis, and the x vector. 
  # 
  # Arguments:
  #  - xvals: the x vector
  #  - nknots: number of knots, either numeric, or NA.  If numeric, basis is 
  #            generated with nknots equally spaced knots.  If NA, a knot is placed
  #            at each unique x-value.
  #  
  #  - degree: degree of the basis
  #  - type: either "poly" for truncated polynomial or "bspl" for a b-spline basis
  #
  ######################################################################
  
  # xvals <- gm$Yearstd
  
  # initilaize the output matrix
  X <- NULL
  
  ##############################
  
  ### find the knot locations
  
  # if a specific number were specified
  if(is.numeric(nknots)) { klocs <- seq(min(xvals),max(xvals),length.out=nknots+2)[-c(1,nknots+2)]
  
  # if a knot at every X-value was specified
  } else if(is.na(nknots)) { 
   
    klocs <- xvals[-length(xvals)] }
  
  ###############################
  
  ### for a b-spline basis
  if(type=="bspl") { X <- bs(xvals,df=NULL,knots=klocs,degree=degree,intercept=T)
  
  ### for a truncated polynomial basis
  } else if (type=="tpoly") {
  
    # generate the initial columns
    for(i in 0:degree) { X <- cbind(X,xvals^i) }
    
    # generate the columns corresponding to the knots
    for(i in klocs) { X <- cbind(X,ifelse( (xvals-i)<0,0,(xvals-i)^degree)) }
  
  }
  
  ################################
  return(list(X=X,nknots=length(klocs),knotlocs=klocs,degree=degree,basis=type,xvals=xvals))
}


##########################################################################################################

Dmatgen <- function(basisgenobj,gamma)
{
  #####################################################
  # Function to generate the penalty matrix D for the sparse fused lasso to be applied
  #  to a spline basis. This consists of a difference matrix D1 to implement the pure fusion
  #  penalty stacked on top of a scaled identity matrix gamma*I, to implement the standard
  #  lasso penalty.  For truncated polynomial basis, no penalty on the first degree+1 columns
  #
  #  Also returns the lasso and fusion penalty matrices.  Designed to only penalize the knot coefficients.
  #  
  #  Arguments;
  #    - basisgenobj: object returned from basis gen
  #    - gamma: the weight parameter that weights the fusion and lasso penalties
  #
  ######################################################
  
  # basisgenobj <- basisgen(gm$Yearstd,nknots=4,degree=2,type="tpoly")
  # gamma <- 0.5
  
  deg <- basisgenobj$degree
  npencol <- basisgenobj$nknots
  
  Dlasso <- cbind(matrix(0,ncol=deg+1,nrow=npencol), diag(1,npencol))
  Dfusion <- cbind(matrix(0,ncol=deg+1,nrow=npencol-1),getD1d(npencol))
  Dfinal <- rbind(Dfusion,gamma*Dlasso)

  return(list(Dlasso=Dlasso,Dfusion=Dfusion,Dfinal=Dfinal))
}
  


# calculate the log determinant
logdet <- function(mat){  return( as.numeric( determinant(mat,logarithm=T)[[1]] ) )}


rl  <- function( vpars,y,X,Z )
{
  ##########################################################################
  # Function to build the restricted likelihood function for the penalized spline
  #  
  # Arguments:
  #   - vpars: sigma2e and lambda in mlm parameterization of pen. spline, note sigma2s = sigma2e/lambda^2
  #   - y: the outcome vector
  #   - X: the fixed effect design matrix
  #   - Z: the random effect design matrix
  ###########################################################################
  
  
  # vpars <- c(1,7)
  
  sigma2e <- vpars[1]
  lambda <- vpars[2]
  
  R <- sigma2e*diag(nrow(X))
  G <- (sigma2e/lambda^2)*diag(ncol(Z))
  
  V <- Z%*%G%*%t(Z) + R
  Vinv <- solve(V)
  xtvinvx <- t(X)%*%Vinv%*%X
  
  out <- -0.5*(  logdet(V) + logdet(xtvinvx) + t(y)%*%( Vinv - Vinv%*%X%*%solve(xtvinvx)%*%t(X)%*%Vinv )%*%y )
  return(as.numeric(-out))
}



aic.genlasso <- function(genlassoobj,smooth.type="lasso")
{
  ###############################################################
  #  For a generalized lasso fit, calculate sigma, degrees of freedom, AICc, 
  #   and Cp for each kink in the solution path
  #
  #  Arguments: 
  #    - genlossoobj: output from genlasso() function, only using the lasso or fusion 
  #        penalty matrices from Dmatgen()
  #    - smooth.type: either 'lasso' or 'fusion' to indicate which penalty was used
  #
  ###############################################################
  
  # genlassoobj <- fit.glasso
  n <- length(genlassoobj$y)
  y <- genlassoobj$y
  
  if(smooth.type=="lasso") {
    
    genlassoobj$beta <- round(genlassoobj$beta,7)
    genlassoobj$df <- apply(genlassoobj$beta,2,function(x) { sum(abs(x)>0) } ) # df = No. of non-zero coefficients
    names(genlassoobj$df) <- NULL
  
  } else if(smooth.type=="fusion") {
    
    genlassoobj$beta <- round(genlassoobj$beta,5)
    genlassoobj$df <- apply(genlassoobj$beta,2,function(x) { length(unique(x))} ) # df = No. of fused groups of coefficicents
    names(genlassoobj$df) <- NULL
  }
    
  # calculate sigma, cp, and aicc
  sigma <- cp <- aicc <- rep(NA,length(genlassoobj$lambda))
  for(i in 1:length(genlassoobj$lambda))
  {
    df <- genlassoobj$df[i]
    sse <- sum((y-genlassoobj$fit[,i])^2)
    sig.tmp <- sqrt(sse / (n-df-1) ) # from https://cran.r-project.org/web/packages/selectiveInference/selectiveInference.pdf
    sigma[i] <- ifelse( is.finite( sig.tmp ),sig.tmp,NA)
    cp[i] <- (1/n)*(sse+2*df*sigma[i]^2) # from https://en.wikipedia.org/wiki/Mallows%27s_Cp
    aicc.tmp <- sse + (2*df*(df+1))/(n-df-1)  # from http://www2.stat.duke.edu/~rcs46/lectures_2015/08-reg2/03Models_v2.pdf
    aicc[i] <- ifelse(!is.na(sigma[i]),aicc.tmp,NA)
  }
  
  genlassoobj$sigma <- sigma
  genlassoobj$aicc <- aicc
  genlassoobj$cp <- cp
  
  return(genlassoobj)
}


fit.spline <- function(y,x,nknots,degree,basis="tpoly",smooth.type="ridge.mlm")
{
  ########################################################
  # Function to fit MLM version of or lasso- or fusion- penalized splines.  For MLM,
  #  model selection is automatic within REML.  For fusion and lasso penalties,
  #  we use AICc to choose the ideal model.  
  # 
  #   - y: the outcome vector
  #   - x: the predictor vector
  #   - nknots: either a positive number K for K equally spaced knots, or NA for a knot at each x-value
  #   - degree: degree of the polynomial basis
  #   - smooth.type: 'ridge.mlm' for MLM implementation of penalized splines, 'lasso' for the lasso-penalized spline,
  #                  'fusion' for the fusion penalized spline
  #
  ########################################################
  
  # x <- gm$Yearstd; y <- gm$tempstd
  # nknots <- NA
  # degree <- 2
  # type <- "tpoly"
  # smooth.type <- "fusion"
  
  
  # generate the X matrix
  Xlist <- basisgen(x,nknots=nknots,degree=degree,type=basis)
  
  ###########################
  
  # define the matrices
  C <- Xlist$X
  X <- Xlist$X[,1:(Xlist$degree+1)]
  Z <- Xlist$X[,(Xlist$degree+2):ncol(Xlist$X)]
  
  ##############################
  
  ##### if the MLM version of splines is requested #####
  if(smooth.type=="ridge.mlm") {
    
    # maximize the restricted likelihood with optim and calculate the Rhat and Ghat matrices
    rlmax <- optim(c(1,5),rl,y=y,X=X,Z=Z,method="L-BFGS-B",lower=c(0.0001,0.0001),upper=c(Inf,Inf))
    parshat <- rlmax$par
    sigma2ehat <- parshat[1]
    
    Rhat <- sigma2ehat*diag(nrow(X))
    Ghat <- (sigma2ehat/parshat[2]^2)*diag(ncol(Z))
    
    ################################
    
    # point estimates of the mean structure and fitted values given estimated variance components
    d <- ncol(C)-ncol(Z)
    blah <- cbind( matrix(0,ncol=d,nrow=ncol(C)), rbind( matrix(0,ncol=ncol(Ghat),nrow=d),solve(Ghat) ))
    Rhatinv <- solve(Rhat)
    
    # point estimates and degrees of freedom
    int <- solve( t(C)%*%Rhatinv%*%C + blah  )%*%t(C)%*%Rhatinv
    df <- sum(diag(C%*%int))
    est <- as.numeric( int%*%y )
    
    # fitted values
    yhat <- as.numeric( C%*%est)
  
  ##############################################################################
  ##############################################################################
  # smooth.type <- "fusion"
  
  ### if it's a lasso or fusion penalty #####
  } else if(smooth.type %in% c("lasso","fusion"))  {
    
    # generate the penalty matrix
    Dmatlist <- Dmatgen(Xlist,gamma=1)
    if(smooth.type=="lasso") { Dmat <- Dmatlist$Dlasso } else if(smooth.type=="fusion") { Dmat <- Dmatlist$Dfusion }
    
    # fit the model 
    fit.glasso <- suppressWarnings( genlasso(y,C,Dmat) )
    fit.glasso <- suppressWarnings( aic.genlasso(fit.glasso,smooth.type=smooth.type) )
    
    # grab the model with lowest aicc
    ind <- which.min(fit.glasso$aicc)
    yhat <- fit.glasso$fit[,ind]
    est <- fit.glasso$beta[,ind]
    df <- fit.glasso$df[ind]
    d <- Xlist$degree+1
  }
  
  ##############################################################################
  ##############################################################################
  
  out.list <- list(fitted=yhat,beta.ur=est[1:d],beta.smooth=est[ (d+1):length(est)],df.fit=df,smooth.type=smooth.type,y=y,xvals=Xlist$xvals)
  # plot(out.list$x,out.list$y);lines(out.list$x,out.list$fitted)
  
  return(out.list)
}



#######################################################################################################3
#######################################################################################################3

# functions for fitting spline of desired degrees of freedom

rl.sigmaeonly <- function( sigma2e,lambda,y,X,Z)
{

  R <- sigma2e*diag(nrow(X))
  G <- (sigma2e/lambda^2)*diag(ncol(Z))
  
  V <- Z%*%G%*%t(Z) + R
  Vinv <- solve(V)
  xtvinvx <- t(X)%*%Vinv%*%X
  
  out <- -0.5*(  logdet(V) + logdet(xtvinvx) + t(y)%*%( Vinv - Vinv%*%X%*%solve(xtvinvx)%*%t(X)%*%Vinv )%*%y )
  return(as.numeric(-out))
  
}


df_from_lambda <- function(lambda,C,X,Z)
{
  #########################################################
  # Function to calculate degrees of freedom from lambda in 
  #  MLM penalized spline
  #
  ########################################################
  
  Ip <- diag(ncol(Z))
  pmat <-  cbind( matrix(0,ncol=ncol(X),nrow=ncol(C)), rbind( matrix(0,ncol=ncol(Ip),nrow=ncol(X)),Ip))
  out <- sum( diag( C %*% solve( t(C)%*%C + (lambda^2)*pmat  )%*%t(C) ) )
  return(out)
}





fit.spline.bydf <- function(y,x,nknots,degree,basis="tpoly",smooth.type="ridge.mlm",df.des)
{
  ########################################################
  # Function to fit MLM version of or lasso- or fusion- penalized splines.  Here we choose the
  #  model by desired degrees of freedom.  
  # 
  #   - y: the outcome vector
  #   - x: the predictor vector
  #   - nknots: either a positive number K for K equally spaced knots, or NA for a knot at each x-value
  #   - degree: degree of the polynomial basis
  #   - smooth.type: 'ridge.mlm' for MLM implementation of penalized splines, 'lasso' for the lasso-penalized spline,
  #                  'fusion' for the fusion penalized spline
  #
  ########################################################
  
  # x <- gm$Yearstd; y <- gm$tempstd
  # nknots <- 50
  # degree <- 2
  # basis <- "tpoly"
  # smooth.type <- "fusion"
  
  
  # generate the X matrix
  Xlist <- basisgen(x,nknots=nknots,degree=degree,type=basis)
  
  ###########################
  
  # define the matrices
  C <- Xlist$X
  X <- Xlist$X[,1:(Xlist$degree+1)]
  Z <- Xlist$X[,(Xlist$degree+2):ncol(Xlist$X)]
  
  ##############################
  
  ##### if the MLM version of splines is requested #####
  if(smooth.type=="ridge.mlm") {
    
    # df.des <- 12
    
    # find the correct value of lambda, this is done by calculating df for a grid of lambdas
    # then using linear interpolation on the grid to find the lambda for desired df
    lambdas <- c( seq(0.0001,0.1,length.out=25),seq(0.11,1,length.out=25),seq(2,100,length.out=25))
    dfs <- sapply(lambdas,df_from_lambda,C,X,Z)
    
    df.diff <- dfs-df.des
    ind.lo <- which.min( df.diff[df.diff > 0] )
    ind.hi <- ind.lo+1
    xlo <- lambdas[ind.lo]
    xhi <- lambdas[ind.hi]
    ylo <- dfs[ind.lo]
    yhi <- dfs[ind.hi]
    
    lambda <- lin.int.x(df.des,xlo,xhi,ylo,yhi)
    
    #####
    # now, fit maximize restricted likelihood with known lambda and unknown sigma2e
    rlmax <- optim(1,rl.sigmaeonly,lambda=lambda,y=y,X=X,Z=Z,method="L-BFGS-B",lower=0.0001,upper=Inf)
    sigma2ehat <- rlmax$par
    
    Rhat <- sigma2ehat*diag(nrow(X))
    Ghat <- (sigma2ehat/lambda^2)*diag(ncol(Z))
    
    ################################
    
    # point estimates of the mean structure and fitted values given estimated variance components
    d <- ncol(C)-ncol(Z)
    blah <- cbind( matrix(0,ncol=d,nrow=ncol(C)), rbind( matrix(0,ncol=ncol(Ghat),nrow=d),solve(Ghat) ))
    Rhatinv <- solve(Rhat)
    
    # point estimates and degrees of freedom
    int <- solve( t(C)%*%Rhatinv%*%C + blah  )%*%t(C)%*%Rhatinv
    df <- sum(diag(C%*%int));df;lambda
    est <- as.numeric( int%*%y )
    
    # fitted values
    yhat <- as.numeric( C%*%est)
  
  ##############################################################################
  ##############################################################################
  # smooth.type <- "fusion"
    
  ### if it's a lasso or fusion penalty #####
  } else if(smooth.type %in% c("lasso","fusion"))  {
    
    # generate the penalty matrix
    Dmatlist <- Dmatgen(Xlist,gamma=1)
    if(smooth.type=="lasso") { Dmat <- Dmatlist$Dlasso } else if(smooth.type=="fusion") { Dmat <- Dmatlist$Dfusion }
    
    # fit the model 
    fit.glasso <- suppressWarnings( genlasso(y,C,Dmat) )
    fit.glasso <- suppressWarnings( aic.genlasso(fit.glasso,smooth.type=smooth.type) )
    
    # grab the model with desired df and lowest aicc
    inds <- which( fit.glasso$df-df.des==0)
    yhatsnew <- fit.glasso$fit[,inds]
    estnew <- fit.glasso$beta[,inds]
    aiccnew <- fit.glasso$aicc[inds]
    
    if(is.matrix(yhatsnew)) { yhat <- yhatsnew[,which.min(aiccnew)] } else { yhat <- yhatsnew }
    if(is.matrix(estnew)) { est <- estnew[,which.min(aiccnew)] } else { est <- estnew}
    df <- df.des
    d <- Xlist$degree+1
  }
  
  ##############################################################################
  ##############################################################################
  
  out.list <- list(fitted=yhat,beta.ur=est[1:d],beta.smooth=est[ (d+1):length(est)],df.fit=df,smooth.type=smooth.type,y=y,xvals=Xlist$xvals)
  # plot(out.list$x,out.list$y);lines(out.list$x,out.list$fitted)
  
  return(out.list)
  
}









###################################################################################################################

# Testing on global mean surface temperature data

# gm <- read.csv(paste0(root,"Global_mean_surface_temperature_data.csv"),header=T)
# 
# # center and standardize the Year variable
# gm$Yearcent <- gm$Year-mean(gm$Year)
# gm$Yearstd <- (gm$Year-mean(gm$Year))/sd(gm$Year)
# gm$tempstd <- (gm$temp.dev-mean(gm$temp.dev))/sd(gm$temp.dev)
# 
# 
# plot(gm$Yearstd,gm$tempstd)
# test <- fit.spline(y=gm$tempstd, x=gm$Yearstd,nknots=NA,degree=2,basis="tpoly",smooth.type="ridge.mlm")
# lines(test$x,test$fitted)
# test <- fit.spline(y=gm$tempstd, x=gm$Yearstd,nknots=NA,degree=2,basis="tpoly",smooth.type="lasso")
# lines(test$x,test$fitted,col="red")
# test <- fit.spline(y=gm$tempstd, x=gm$Yearstd,nknots=NA,degree=2,basis="tpoly",smooth.type="fusion")
# lines(test$x,test$fitted,col="blue")

##############################################################################################################

### Code to generate the true functions for the simulation study ######

# bathtub function
bathtub1 <- function(x)
{
  if(x<1 | x>=4 ) { y <- 1 
  } else if(x >=1 & x<1.4) { y <- exp(-10*(x-1)) 
  } else if(x>1.4 & x<=3.6) { y <- exp(-10*0.4) 
  } else if(x>3.6 & x<4) { y <- exp(-10*(-x+4)) }
  return(y)
}


# doppler style function
doppler1 <- function(x)
{
  if(x>=3.8){ y <- 1.07-0.15*x
  } else if(x<3.8 ){ y <- doppler( (x+1.2)/5)+0.5 }
  return(y)
}

# add gaussian error
adderror <- function(truevals,sigma=0.15) { return(truevals + rnorm(length(truevals),mean=0,sd=sigma) ) }


###################################################################################################################






