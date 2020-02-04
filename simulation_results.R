
library(abind)
library(splines)
library(genlasso)
library(wavethresh)


root <- "C:/Users/rbrow/OneDrive/Documents/Public Github/Sparse-Spline-Penalties/"
source(paste0(root,"spline_fitting_functions.R"))

load(paste0(root,"Simulation Results/bt_15kn.rdata"))
load(paste0(root,"Simulation Results/bt_50kn.rdata"))
load(paste0(root,"Simulation Results/bt_eachx.rdata"))
load(paste0(root,"Simulation Results/dp_15kn.rdata"))
load(paste0(root,"Simulation Results/dp_50kn.rdata"))
load(paste0(root,"Simulation Results/dp_eachx.rdata"))


# combine into single 5-d array and remove the saller arrays
bt <- abind(bt_15kn,bt_50kn,bt_eachx,along=4)
dp <- abind(dp_15kn,dp_50kn,dp_eachx,along=4)
sdat <- abind(bt,dp,along=5)
names(dimnames(sdat))[1] <- "method"
names(dimnames(sdat))[2] <- "x"
names(dimnames(sdat))[3] <- "simiter"
dimnames(sdat)[[4]] <- c("15kn","50kn","eachx");names(dimnames(sdat))[4] <- "knots"
dimnames(sdat)[[5]] <- c("bt1","dp1"); names(dimnames(sdat))[5] <- "fxn"
rm(list=c("bt","bt_15kn","bt_50kn","bt_eachx","dp","dp_15kn","dp_50kn","dp_eachx"))

# reorder the array
sdat <- aperm(sdat,c(1,4,5,3,2))


#################################################################################################
#################################################################################################
#################################################################################################
export <- F


## True Functions and Evaluation Locations
xplot <- seq(0,5,length.out=1000)
dpplot <- sapply(xplot,doppler1)
btplot <- sapply(xplot,bathtub1)
xlocs.bt <- c(0.25,1,1.25,2.5,3.75,4,4.75)
xlocs.dp <- c(0.17,0.47,0.9,1.2,1.56,2.05,2.7,4.25)


if(export) { pdf(paste0(root,"Slides/Final Talk/truefunx_withpts.pdf"),width=8,height=6)}
plot("n",ylab="",xlab="",main="True Functions",xlim=c(0,5),ylim=c(0,1))

lines(xplot,btplot,lwd=1.5,col="yellowgreen")
ylocs.bt <- sapply(xlocs.bt,bathtub1)
segments(xlocs.bt,ylocs.bt-0.015,xlocs.bt,ylocs.bt+0.015,lwd=2)

lines(xplot,dpplot,lwd=1.5,col="deepskyblue2")
ylocs.dp <- sapply(xlocs.dp,doppler1)
segments(xlocs.dp,ylocs.dp-0.015,xlocs.dp,ylocs.dp+0.015,lwd=2)

legend("bottomright",c("Bathtub","Doppler"),lwd=1.5,col=c("yellowgreen","deepskyblue2"),bty="n" )
if(export){dev.off()}


## true functions without evaluation locations
if(export) { pdf(paste0(root,"Slides/Final Talk/truefunx.pdf"),width=8,height=6)}
plot("n",ylab="",xlab="",main="True Functions",xlim=c(0,5),ylim=c(0,1))

lines(xplot,btplot,lwd=1.5,col="yellowgreen")
ylocs.bt <- sapply(xlocs.bt,bathtub1)
# segments(xlocs.bt,ylocs.bt-0.015,xlocs.bt,ylocs.bt+0.015,lwd=2)

lines(xplot,dpplot,lwd=1.5,col="deepskyblue2")
ylocs.dp <- sapply(xlocs.dp,doppler1)
# segments(xlocs.dp,ylocs.dp-0.015,xlocs.dp,ylocs.dp+0.015,lwd=2)

legend("bottomright",c("Bathtub","Doppler"),lwd=1.5,col=c("yellowgreen","deepskyblue2"),bty="n" )
if(export){dev.off()}


######################################

## Separate Plots ##
if(export) { pdf(paste0(root,"Slides/Final Talk/bathtub_ticks.pdf"),width=6,height=6)}
plot("n",ylab="",xlab="",main="Bathtub",xlim=c(0,5),ylim=c(0,1))
lines(xplot,btplot,lwd=1,col="yellowgreen")
ylocs.bt <- sapply(xlocs.bt,bathtub1)
segments(xlocs.bt,ylocs.bt-0.015,xlocs.bt,ylocs.bt+0.015,lwd=2)
if(export) { dev.off() }

if(export) { pdf(paste0(root,"Slides/Final Talk/doppler_ticks.pdf"),width=6,height=6)}
plot("n",ylab="",xlab="",main="Doppler",xlim=c(0,5),ylim=c(0,1))
lines(xplot,dpplot,lwd=1,col="deepskyblue2")
ylocs.dp <- sapply(xlocs.dp,doppler1)
segments(xlocs.dp,ylocs.dp-0.015,xlocs.dp,ylocs.dp+0.015,lwd=2)
if(export) { dev.off()}

## Bathtub Plots with simulated data ##
if(export) { pdf(paste0(root,"Slides/Final Talk/bt_withpts.pdf"),width=8,height=6)}
set.seed(123)
plot("n",ylab="",xlab="",main="Bathtub",xlim=c(0,5),ylim=c(-0.3,1.3))
x <- seq(0,5,length.out=150)
y <- adderror( sapply(x,bathtub1))
points(x,y,cex=0.8,pch=20)
lines(xplot,btplot,lwd=1.5,col="yellowgreen")

if(export){dev.off()}

###########################################


### Plots with differing degrees of freedom ####

set.seed(1234)
x <- seq(0,5,length.out=150)
true <- sapply(x,doppler1)
y <- adderror( sapply(x,doppler1))


if(export) { png(paste0(root,"Slides/Final Talk/dp_bydf/trueonly.png"),width=8,height=6,units="in",res=1000) } 
plot("n",ylab="",xlab="",main="Doppler: True Function",xlim=c(0,5),ylim=c(-0.3,1.3))
points(x,y,cex=0.8,pch=20)
lines(x,true,lwd=2,col="grey80")
if(export) {dev.off()}

for ( i in c(6,8,10,12,16) )
{
  if(export) { png(paste0(root,"Slides/Final Talk/dp_bydf/",i,".png"),width=8,height=6,units="in",res=1000) } 
  plot("n",ylab="",xlab="",main=paste0("Doppler: DF = ",i),xlim=c(0,5),ylim=c(-0.3,1.3))
  points(x,y,cex=0.8,pch=20)
  yhatmlm <- fit.spline.bydf(y,x,50,3,basis="tpoly",smooth.type="ridge.mlm",df.des=i)$fitted
  yhatlas <- fit.spline.bydf(y,x,50,3,basis="tpoly",smooth.type="lasso",df.des=i)$fitted
  yhatfus <- fit.spline.bydf(y,x,50,3,basis="tpoly",smooth.type="fusion",df.des=i)$fitted
  
  lines(x,true,lwd=2,col="grey80")
  lines(x,yhatmlm,lwd=1.5,col="orange")
  lines(x,yhatlas,lwd=1.5,col="yellowgreen")
  lines(x,yhatfus,lwd=1.5,col="deepskyblue2")

  legend( "topright",c("True","MLM","Lasso","Fusion"),col=c("grey80","orange","yellowgreen","deepskyblue"),
          lwd=c(2,1.5,1.5,1.5),bty="n")
  if(export) {dev.off()}
}


if(export) { png(paste0(root,"Slides/Final Talk/dp_bydf/selected.png"),width=8,height=6,units="in",res=1000) } 
if(export) { pdf(paste0(root,"Slides/Final Talk/dp_bydf/selected.pdf"),width=8,height=6) } 
plot("n",ylab="",xlab="",main=paste0("Doppler: Selected Model"),xlim=c(0,5),ylim=c(-0.3,1.3))
points(x,y,cex=0.8,pch=20)


mlm <- fit.spline(y,x,50,3,basis="tpoly",smooth.type="ridge.mlm")
las <- fit.spline(y,x,50,3,basis="tpoly",smooth.type="lasso")
fus <- fit.spline(y,x,50,3,basis="tpoly",smooth.type="fusion")

yhatmlm <- mlm$fitted
yhatlas <- las$fitted
yhatfus <- fus$fitted

lines(x,true,lwd=2,col="grey80")
lines(x,yhatmlm,lwd=1.5,col="orange")
lines(x,yhatlas,lwd=1.5,col="yellowgreen")
lines(x,yhatfus,lwd=1.5,col="deepskyblue2")

legend( "topright",c("True",paste0( c("MLM","Lasso","Fusion"),": DF = ",c(round(mlm$df.fit),las$df.fit,fus$df.fit))),
                           col=c("grey80","orange","yellowgreen","deepskyblue"),lwd=c(2,1.5,1.5,1.5),bty="n")


if(export) {dev.off()}





### Plots with differing degrees of freedom ####

set.seed(124)
x <- seq(0,5,length.out=150)
true <- sapply(x,bathtub1)
y <- adderror( sapply(x,bathtub1))


if(export) { pdf(paste0(root,"Slides/Final Talk/dp_bydf/dp_bydf.pdf"),width=8,height=6) } 
plot("n",ylab="",xlab="",main="Bathtub: True Function",xlim=c(0,5),ylim=c(-0.3,1.3))
points(x,y,cex=0.8,pch=20)
lines(x,true,lwd=2,col="grey80")
# if(export) {dev.off()}

for ( i in c(6,8,10,12,16) )
{
  # if(export) { png(paste0(root,"Slides/Final Talk/dp_bydf/",i,".png"),width=8,height=6,units="in",res=1000) } 
  plot("n",ylab="",xlab="",main=paste0("Bathtub: DF = ",i),xlim=c(0,5),ylim=c(-0.3,1.3))
  points(x,y,cex=0.8,pch=20)
  yhatmlm <- fit.spline.bydf(y,x,50,3,basis="tpoly",smooth.type="ridge.mlm",df.des=i)$fitted
  yhatlas <- fit.spline.bydf(y,x,50,3,basis="tpoly",smooth.type="lasso",df.des=i)$fitted
  yhatfus <- fit.spline.bydf(y,x,50,3,basis="tpoly",smooth.type="fusion",df.des=i)$fitted
  
  lines(x,true,lwd=2,col="grey80")
  lines(x,yhatmlm,lwd=1.5,col="orange")
  lines(x,yhatlas,lwd=1.5,col="yellowgreen")
  lines(x,yhatfus,lwd=1.5,col="deepskyblue2")
  
  legend( "bottomright",c("True","MLM","Lasso","Fusion"),col=c("grey80","orange","yellowgreen","deepskyblue"),
          lwd=c(2,1.5,1.5,1.5),bty="n")
  # if(export) {dev.off()}
}


# if(export) { png(paste0(root,"Slides/Final Talk/dp_bydf/selected.png"),width=8,height=6,units="in",res=1000) } 
if(export) { pdf(paste0(root,"Slides/Final Talk/dp_bydf/selected.pdf"),width=8,height=6) }
plot("n",ylab="",xlab="",main=paste0("Bathtub: Selected Model"),xlim=c(0,5),ylim=c(-0.3,1.3))
points(x,y,cex=0.8,pch=20)


mlm <- fit.spline(y,x,50,3,basis="tpoly",smooth.type="ridge.mlm")
las <- fit.spline(y,x,NA,3,basis="tpoly",smooth.type="lasso")
fus <- fit.spline(y,x,50,3,basis="tpoly",smooth.type="fusion")

yhatmlm <- mlm$fitted
yhatlas <- las$fitted
yhatfus <- fus$fitted

lines(x,true,lwd=2,col="grey80")
lines(x,yhatmlm,lwd=1.5,col="orange")
lines(x,yhatlas,lwd=1.5,col="yellowgreen")
lines(x,yhatfus,lwd=1.5,col="deepskyblue2")

legend( "bottomright",c("True",paste0( c("MLM","Lasso","Fusion"),": DF = ",c(round(mlm$df.fit),las$df.fit,fus$df.fit))),
        col=c("grey80","orange","yellowgreen","deepskyblue"),lwd=c(2,1.5,1.5,1.5),bty="n")


if(export) {dev.off()}





#################################################################################################
#################################################################################################
#################################################################################################


yhat.calc <- function(fitvec,xloc,xvec)
{
  ###############################################################
  # Given a particular x-value, the x locations for the fitted values,
  #  and the fitted values, calculate the value of the predicted function
  #  at the x-value.  Finds the x locations and correspondingn fitted values
  #  on either side of the desired x value, and gets yhat via linear interpolation
  #  
  #  Arguments:
  #    - xloc: the x-value for which a predicted yhat is desired
  #    - xvec: the x vector used in the penalized spline
  #    - fitvec:  the fitted values corresponding to the x vector
  #
  ###############################################################
  
  # xloc <- 0.25
  # xvec <- x
  # fitvec <- sdat["lasso","50kn","bt1",1,]

  dvec <- xvec-xloc
  
  if(0 %in% dvec){ yhat <- fitvec[which(dvec==0)]
  
  } else {
  
  ind.lo <- which.max( dvec[dvec < 0] )
  ind.hi <- ind.lo+1
  xlo <- xvec[ind.lo]
  xhi <- xvec[ind.hi]
  ylo <- fitvec[ind.lo]
  yhi <- fitvec[ind.hi]
  
  yhat <- lin.int.y(targx=xloc,xlo,xhi,ylo,yhi)
  
  }
  
  return(yhat)
}


#################################################################################################
#################################################################################################
#################################################################################################


x <- seq(0,5,length.out=150)


calc.bvm <- function(targx,method,nknots,fxn,xvec,sdat)
{
  # targx <- 0.25
  # method <- "mlm"
  # nknots <- "50kn"
  # fxn <- "bt1"
  # xvec <- x
  
  # calculate the true y-value
  if(fxn=="bt1") {truey <- bathtub1(targx)} else if(fxn=="dp1") {truey <- doppler1(targx)}
  
  # grab the fitted values for the 1000 simulations from the specified settings
  tmp <- sdat[method,nknots,fxn,,]
  
  # calcualte the predicted values, bias, variance, and mse
  yhats <- apply(tmp,1,yhat.calc,xloc=targx,xvec=xvec)
  bias <- mean(yhats-truey,na.rm=T)
  vari <- var(yhats,na.rm=T)
  mse <- bias^2 + vari
  
  return(c(bias,vari,mse))
  
}  


### Generate dataframe of all scenarios and x locations to calculate bias, variance, and mse
method <- c("mlm","lasso","fusion")
knots <- c("15kn","50kn","eachx")

sitmat.bt <- expand.grid(method=method,knots=knots)
sitmat.bt$fxn <- "bt1"
sitmat.bt <- merge(sitmat.bt,data.frame(xloc=xlocs.bt,fxn="bt1"),all.x=T)

sitmat.dp <- expand.grid(method=method,knots=knots)
sitmat.dp$fxn <- "dp1"
sitmat.dp <- merge(sitmat.dp,data.frame(xloc=xlocs.dp,fxn="dp1"),all.x=T)

sitmat <- rbind(sitmat.bt,sitmat.dp)





### calculate bias, variance and mse for all scenarios and x locations
outmat <- NULL
for(i in 1:nrow(sitmat))
{
  out <- calc.bvm(targx=sitmat$xloc[i],method=sitmat$method[i],nknots=sitmat$knots[i],fxn=sitmat$fxn[i],xvec=x,sdat=sdat)
  outmat <- rbind(outmat,out)
}

sitmat <- cbind(sitmat,data.frame(outmat))
names(sitmat)[5:7] <- c("bias","var","mse")


# calculate mse, summed over all original 150 x locations and other regions
scen.mat <- unique(sitmat[,1:3])
mses.overall <- mses.gt25 <- mses.le25 <- mses.17to33 <-  mses.08to155 <- c()
for(i in 1:nrow(scen.mat))
{
  blah <- sapply(x, calc.bvm,method=scen.mat$method[i],nknots=scen.mat$knots[i],fxn=scen.mat$fxn[i],x,sdat)
  mse.sum.overall <- rowSums(blah)[3]
  mse.sum.gt25 <- rowSums(blah[,x>2.5])[3]
  mse.sum.le25 <- rowSums(blah[,x<=2.5])[3]
  mse.sum.0.8to1.55 <- rowSums(blah[,x<=1.55 & x>=0.8])[3]
  mse.sum.1.7to3.3 <- rowSums(blah[,x<=3.3 & x>=1.7])[3]
  
  
  mses.overall <- c(mses.overall,mse.sum.overall)
  mses.gt25 <- c(mses.gt25,mse.sum.gt25)
  mses.le25 <- c(mses.le25,mse.sum.le25)
  mses.08to155 <- c(mses.08to155,mse.sum.0.8to1.55)
  mses.17to33 <- c(mses.17to33,mse.sum.1.7to3.3)
}

scen.mat$mse.overall <- mses.overall
scen.mat$mse.gt25 <- mses.gt25
scen.mat$mse.le25 <- mses.le25
scen.mat$mse.08to155 <- mses.08to155
scen.mat$mse.17to33 <- mses.17to33




## cumulative MSE tables in paper
scen.mat[scen.mat$fxn=="dp1",c(2,3,6,5,4)]
scen.mat[scen.mat$fxn=="bt1",c(2,3,7,8,4)]




## Bias-Variance-MSE tables in paper
sitmat[sitmat$fxn=="dp1" & sitmat$xloc==xlocs.dp[1],]
sitmat[sitmat$fxn=="dp1" & sitmat$xloc==xlocs.dp[8],]

sitmat[sitmat$fxn=="bt1" & sitmat$xloc==xlocs.bt[2],]
sitmat[sitmat$fxn=="bt1" & sitmat$xloc==xlocs.bt[3],]
sitmat[sitmat$fxn=="bt1" & sitmat$xloc==xlocs.bt[4],]











