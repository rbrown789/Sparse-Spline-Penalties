
library(splines)
library(genlasso)
library(wavethresh)

root <- "C:/Users/rbrow/OneDrive/Documents/Public Github/Sparse-Spline-Penalties/"

source(paste0(root,"spline_fitting_functions.R"))
gm <- read.csv(paste0(root,"gmst_data.csv"),header=T)

#################################################################################################
# center and standardize the Year variable
gm$Yearcent <- gm$Year-mean(gm$Year)
gm$Yearstd <- (gm$Year-mean(gm$Year))/sd(gm$Year)
gm$tempstd <- (gm$temp.dev-mean(gm$temp.dev))/sd(gm$temp.dev)



###### Fits to the Global Mean Surface Temperature Data #############

tp30 <- basisgen(gm$Yearstd,nknots=NA,degree=2,type="tpoly")
Dtp30 <- Dmatgen(tp30,gamma=1)


testlasso <- genlasso(y=gm$tempstd,X=tp30$X,D=Dtp30$Dlasso)
testfusion <- genlasso(y=gm$tempstd,X=tp30$X,D=Dtp30$Dfusion)


fuseind <- 50
lassind <- 50


betafuse <- round( testfusion$beta[,fuseind],7); length(unique(betafuse))
betalasso <- round( testlasso$beta[,lassind],7); sum(abs(betalasso) > 0)


export <- T
export <- F

if(export){ pdf(paste0(root,"wip1plot2.pdf"),width=9,height=6) }
par(mfrow = c(2,2),
    oma = c(4,4,0,2) + 0.1,
    mar = c(0,0,1.5,1.5) + 0.1,
    xpd=FALSE)

cex <- 0.65; pch <- 20; xlim <- c(-1.7,1.7); ylim2 <- c(-2.3,2.3)

plot(gm$Yearstd,gm$tempstd,axes=F,xlab="",ylab="Standardized Temperature",main="Lasso Penalty: df=11",
     xlim=xlim,pch=pch,cex=cex)
box();axis(2)
points(gm$Yearstd,tp30$X%*%betalasso,type="l",col="yellowgreen")


plot(gm$Yearstd,gm$tempstd,axes=F,xlab="",ylab="",main="Fusion Penalty: df=11",
     xlim=xlim,pch=pch,cex=cex)
box();axis(4)
points(gm$Yearstd,tp30$X%*%betafuse,type="l",col="deepskyblue2")


plot("n",xlab="Standarized Years",ylab="Coefficient Value", xlim=xlim,ylim=ylim2)
abline(h=0,col="grey60")
points(tp30$knotlocs,testlasso$beta[,lassind][-(1:(tp30$degree+1) )],pch=pch,cex=cex)

plot("n",xlab="Standarized Years",ylab="",axes=F, xlim=xlim,ylim=ylim2)
box();axis(4);axis(1)
abline(h=0,col="grey60")
points(tp30$knotlocs,testfusion$beta[,fuseind][-(1:(tp30$degree+1) )],pch=pch,cex=cex)

title(xlab = "Standardized Years",
      ylab = "    Coefficient Value                        Standardized Temperature       ",
      outer = TRUE, line = 2,cex.lab=1)


if(export){dev.off()}

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################






