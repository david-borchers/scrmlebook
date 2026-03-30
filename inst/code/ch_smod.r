library(secr)
library(fields) # needed for secr.surface
library(maptools) # needed for secr.surface
library(sp)
library(raster)
library(scrmlebook)

#setwd("/Users/dlb/MeetingsTalks/SCR Workshop St Andrews 2016/slides/")
#setwd("/Users/dlb/MeetingsTalks/ISEC_Seattle_2016/Workshop/")
#load("/Users/dlb/Dropbox/kruger (1)/workshop-fits.RData")
#load("/Users/dlb/Dropbox/kruger (1)/popn.RData")
#kruger.popn=popn
#save(kruger.popn,file="./scrmlebook/data/kruger-popn.RData")


data("ERdetfun-data")
data("ERdetfun-fits")
data("kruger-popn")

# crop kruger boundary to focus plots on area with data
cpoly = crop(attributes(kruger.mask)$polygon,extent(4079702,4144000,-3328069,-3257132))
attributes(kruger.mask)$polygon = cpoly

habitat.model=list(D~habitat.cov,lambda0~1,sigma~1)
habitat.fit=secr.fit(kruger.capt,model=habitat.model,mask=kruger.mask,detectfn="HHN")

Ds = predictDsurface(habitat.fit) # get density surface from fitted model in fit
a = attributes(kruger.mask)$area # get cell area
E.N = sum(attributes(Ds)$area*covariates(Ds)$D.0) # calculate "integral"
region.N(habitat.fit) # this calculates the integral for you (instead of doing line above)


pdf("./keepfigure/KrugerDs.pdf",h=5,w=10)
par(mfrow=c(1,2))
Dplot=plot(Ds,key=FALSE,bty="n",xaxt="n",yaxt="n",xlab="Easting",ylab="Northing",col=tim.colors(40),asp=1)
points(kruger.cams,pch="+",col="white")
persp(Dplot,phi=65,theta=3,expand=0.75,r=1,d=2,xlab="Easting",ylab="Northing",zlab="Density")
dev.off()

## create population
#set.seed(12345)
#popn = sim.popn("D.0",Ds,model2D="IHP",poly=attributes(kruger.mask)$polygon,Nsessions=1)

ylim=c(min(kruger.popn$y),max(kruger.popn$y));ylim[2]=ylim[2]+0.05*abs(diff(ylim))
xlim=c(min(kruger.popn$x),max(kruger.popn$x));xlim[2]=xlim[2]+0.05*abs(diff(xlim))
asp=abs(diff(ylim))/abs(diff(xlim))
sex=as.matrix(covariates(kruger.popn))
pch=rep(19,length(sex))
pch[sex=="F"]=17
pdf("./keepfigure/KrugerPop.pdf",h=6*asp*0.8,w=6)
#par(mar=c(4,4,3,1))
plot(kruger.popn$x,kruger.popn$y,pch=19,cex=1.25,xlim=xlim,ylim=ylim,
     bty="n",xaxt="n",yaxt="n",xlab="Easting",ylab="Northing",asp=1)
plot(attr(Ds,"polygon"),xlim=xlim,add=TRUE)
points(kruger.cams,pch="+",col="darkgray")
dev.off()


pdf("./keepfigure/KrugerHabitat.pdf",h=5,w=10)
par(mfrow=c(1,2))
Dplot=plot(Ds,covariate="habitat.cov",contour=FALSE,bty="n",xaxt="n",yaxt="n",xlab="Easting",ylab="Northing",col=terrain.colors(40))
dev.off()


# This is a quick-fix. Really need the data on detections from Ben
# --------------------
ch=sim.capthist(kruger.cams,kruger.popn,detectfn="HHN",renumber=FALSE,seed=12345,noccasions=1,
                detectpar=list(lambda0=exp(coefficients(habitat.fit)["lambda0",1]),
                               sigma=exp(coefficients(habitat.fit)["sigma",1])))
seen = dimnames(ch)[[1]]
kruger.dets = kruger.popn[seen,]
pdf("./keepfigure/KrugerThinnedPop.pdf",h=6*asp*0.8,w=6)
plot(kruger.popn$x,kruger.popn$y,pch=1,cex=1.25,xlim=xlim,ylim=ylim,
     bty="n",xaxt="n",yaxt="n",xlab="Easting",ylab="Northing",asp=1)
points(kruger.dets$x,kruger.dets$y,pch=19,cex=1.25)
plot(attr(Ds,"polygon"),xlim=xlim,add=TRUE)
points(kruger.cams,pch="+",col="darkgray")
dev.off()


# Thinning plots:
# ==============


hrhn <- function(d, pars){
  pars[1]*exp(-d^2/(2*pars[2]^2))
}
lambda0.1=0.5
sigma.1=5450
lambda0.2=2
sigma.2=850
pars1=c(lambda0.1,sigma.1)
pars2=c(lambda0.2,sigma.2)
#mask=mix.fit$mask ##
ncams=dim(kruger.cams)[1]
d=matrix(rep(NA,dim(kruger.mask)[1]*ncams),nrow=ncams)
for(i in 1:ncams) {
  d[i,]=sqrt((kruger.cams$x[i]-kruger.mask$x)^2 + (kruger.cams$y[i]-kruger.mask$y)^2)
}
h=hrhn(d,pars1)
H=apply(h,2,sum)
H[is.na(H)]=0

# Thinning 
# --------

# thinning surfaces
dethaz=data.frame(x=kruger.mask$x,y=kruger.mask$y,z=H)
p=q=dethaz
p$z=1-exp(-dethaz$z)
q$z=exp(-dethaz$z)
# Thinned density surfaces
Dp=Dq=p
Dp$z=p$z*covariates(Ds)$D.0
Dq$z=q$z*covariates(Ds)$D.0
# boundary polygon
bound=attributes(Ds)$polygon

xrange=range(kruger.mask$x);xrange
yrange=range(kruger.mask$y);yrange
window.x.range=c(xrange[1],xrange[1]+3*diff(xrange))
window.y.range=c(yrange[1],yrange[1]+2*diff(yrange)*1.2)
yadd1=diff(yrange)*1.3
xadd1=diff(xrange)
xadd3=2*diff(xrange)

plot1=Ds
plot1$x=plot1$x+xadd1
plot1$y=plot1$y+yadd1
cams1=kruger.cams
cams1$x=cams1$x+xadd1
cams1$y=cams1$y+yadd1
bound1=bound
bound1@polygons[[1]]@Polygons[[1]]@coords[,1]=bound1@polygons[[1]]@Polygons[[1]]@coords[,1]+xadd1
bound1@polygons[[1]]@Polygons[[1]]@coords[,2]=bound1@polygons[[1]]@Polygons[[1]]@coords[,2]+yadd1


plot3=Dq
plot3$x=plot3$x+xadd3
cams3=kruger.cams
cams3$x=cams3$x+xadd3
bound3=bound
bound3@polygons[[1]]@Polygons[[1]]@coords[,1]=bound3@polygons[[1]]@Polygons[[1]]@coords[,1]+xadd3

zlim=range(c(0,covariates(Ds)$D.0))

quartz(h=8,w=7)
#par(mar=c(4,2,2,8))
# Density plot
Dplot=plot(plot1,main=expression(D(s)),xaxt="n",yaxt="n",bty="n",xlab="",ylab="",asp=1,addpoly=FALSE,
           col=tim.colors(30),xlim=window.x.range,ylim=window.y.range,zlim=zlim,key=FALSE)
plot(bound1,add=TRUE)
points(cams1,pch="+",col="white")
# Two thinned plots
pest=prep4image(Dp,main=expression(D(bold(s))*p(bold(s))),xaxt="n",yaxt="n",bty="n",xlab=expression(D(bold(s))*p(bold(s))),ylab="",asp=1,
                ,col=tim.colors(30),key=FALSE,add=TRUE,zlim=zlim)
plot(bound,add=TRUE)
points(kruger.cams,pch="+",col="white")

pest=prep4image(plot3,main=expression(D(bold(s))*(1-p(bold(s)))),xaxt="n",yaxt="n",bty="n",xlab=expression(D(bold(s))*(1-p(bold(s)))),ylab="",asp=1,
                ,col=tim.colors(30),key=FALSE,add=TRUE,zlim=zlim)
plot(bound3,add=TRUE)
points(cams3,pch="+",col="white")



#pdf("krugerthinning.pdf",h=9,w=4)
par(mar=c(2,2,2,2),mfrow=c(3,1))
Dplot=plot(Ds,main=expression(D(s)),xaxt="n",yaxt="n",bty="n",xlab="",ylab="",asp=1,addpoly=FALSE)
points(kruger.cams,pch="+",col="white")


pest=prep4image(p,main=expression(p(bold(s))),xaxt="n",yaxt="n",bty="n",xlab="",ylab="",asp=1)
points(kruger.cams,pch="+",col="white")

pest=prep4image(q,main=expression(1-p(bold(s))),xaxt="n",yaxt="n",bty="n",xlab="",ylab="",asp=1)
points(kruger.cams,pch="+",col="white")

zlim=range(Dp$z,Dq$z)
pest=prep4image(Dp,main=expression(D(bold(s))*p(bold(s))),xaxt="n",yaxt="n",bty="n",xlab="",ylab="",asp=1)
points(kruger.cams,pch="+",col="white")

pest=prep4image(Dq,main=expression(D(bold(s))*(1-p(bold(s)))),xaxt="n",yaxt="n",bty="n",xlab="",ylab="",asp=1)
points(kruger.cams,pch="+",col="white")


dev.off()


# INVERSE Thinning
# ----------------
pdf("kruger0thinning.pdf",h=9,w=4)
par(mar=c(2,2,2,2),mfrow=c(3,1))
Dest=prep4image(Ddf,main=expression(D(x)),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))

dethaz=data.frame(x=mask$x,y=mask$y,z=H)
p0=dethaz
p0$z=exp(-dethaz$z)
p0est=prep4image(p0,main=expression(1-p(x)),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))

Dp0=p0
Dp0$z=Dp0$z*Ddf$z
pest=prep4image(Dp0,main=expression(D(x)*(1-p(x))),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))

dev.off()




#===================== Stuff below here is old and as yet unused for book ======================


plot(popn$x,popn$y,pch=19,cex=0.5,asp=1,xlim=xlim,ylim=ylim)

Dfit=predictDsurface(mix.fit)
plot(Dfit,border=0,legend=FALSE)
krugerCH=mix.fit$capthist
krugerCAMS=traps(mix.fit$capthist)
krugerMASK=mix.fit$mask
plot(krugerCAMS,add=TRUE)

#save(krugerCAMS,krugerCH,krugerMASK,prep4image,secr.surface,file="/Users/dlb/MeetingsTalks/ISEC_Seattle_2016/Workshop/RforParticipants/Lectures1&2.rda")
#rm(list=c("krugerCAMS","krugerCH","krugerMASK","prep4image","secr.surface"))

load("/Users/dlb/MeetingsTalks/ISEC_Seattle_2016/Workshop/RforParticipants/Lectures1&2.rda")

Ddf=data.frame(x=Dfit$x,y=Dfit$y,z=covariates(Dfit)$D.0)
covdf=data.frame(x=Dfit$x,y=Dfit$y,z=covariates(Dfit)$habitat.cov)
Dest=prep4image(Ddf)
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCH,tracks=TRUE,add=TRUE,rad=100000)
plot(krugerCAMS,add=TRUE)

habcov=prep4image(covdf,col=terrain.colors(40))
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCH,tracks=TRUE,add=TRUE,rad=100000)
plot(krugerCAMS,add=TRUE)

persp(Dest,phi=65,theta=3,expand=0.75,r=1,d=2,xlab="Easting",ylab="Northing",zlab="Density")
persp(habcov,phi=65,theta=3,expand=0.75,r=1,d=2,xlab="Easting",ylab="Northing",zlab="Habitat")


summary(krugerCH)
covariates(krugerCH)
plot(krugerCH,tracks=TRUE,rad=1000,type="petal",gridlines=FALSE)
plot(krugerCAMS,add=TRUE)

dy=0.05*abs(diff(range(popn$y)))
cex=1
pdf("popn.pdf")
par(mar=c(0,0,0,3))
plot(popn$x,popn$y,ylim=c(min(popn$y)-dy,max(popn$y)+dy),
     xaxt="n",yaxt="n",bty="n",xlab="",ylab="",cex=cex)
plot(attributes(Dfit)$polygon,add=TRUE)
dev.off()

pdf("density.pdf")
par(mar=c(0,0,0,2))
Dest=prep4image(Ddf,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",cex=cex)
plot(attributes(Dfit)$polygon,add=TRUE)
points(popn$x,popn$y,cex=cex)
dev.off()


hrhn <- function(d, pars){
  pars[1]*exp(-d^2/(2*pars[2]^2))
}
lambda0.1=0.5
sigma.1=5450
lambda0.2=2
sigma.2=850
pars1=c(lambda0.1,sigma.1)
pars2=c(lambda0.2,sigma.2)
mask=mix.fit$mask
ncams=dim(krugerCAMS)[1]
d=matrix(rep(NA,dim(mask)[1]*ncams),nrow=ncams)
for(i in 1:ncams) {
  d[i,]=sqrt((krugerCAMS$x[i]-mask$x)^2 + (krugerCAMS$y[i]-mask$y)^2)
}
h=hrhn(d,pars1)
H=apply(h,2,sum)
H[is.na(H)]=0

# Thinning plots
# --------------
pdf("krugerthinning.pdf",h=9,w=4)
par(mar=c(2,2,2,2),mfrow=c(3,1))
Dest=prep4image(Ddf,main=expression(D(x)),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))

dethaz=data.frame(x=mask$x,y=mask$y,z=H)
#Hest=prep4image(dethaz,main="Detection rate",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
#plot(attributes(Dfit)$polygon,add=TRUE)
#plot(krugerCAMS,add=TRUE,detpar=list(col="black"))

#En=dethaz
#En$z=En$z*Ddf$z
#pest=prep4image(En,main="Detection probability",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
#plot(attributes(Dfit)$polygon,add=TRUE)
#plot(krugerCAMS,add=TRUE,detpar=list(col="black"))

p=dethaz
p$z=1-exp(-dethaz$z)
pest=prep4image(p,main=expression(p(x)),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))

Dp=p
Dp$z=Dp$z*Ddf$z
pest=prep4image(Dp,main=expression(D(x)*p(x)),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))

dev.off()


# INVERSE Thinning plots
# -----------------------
pdf("kruger0thinning.pdf",h=9,w=4)
par(mar=c(2,2,2,2),mfrow=c(3,1))
Dest=prep4image(Ddf,main=expression(D(x)),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))

dethaz=data.frame(x=mask$x,y=mask$y,z=H)
p0=dethaz
p0$z=exp(-dethaz$z)
p0est=prep4image(p0,main=expression(1-p(x)),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))

Dp0=p0
Dp0$z=Dp0$z*Ddf$z
pest=prep4image(Dp0,main=expression(D(x)*(1-p(x))),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))

dev.off()

# plot for mask slide
mask=mix.fit$mask
dy=0.05*abs(diff(range(popn$y)))
pdf("krugermask.pdf")
plot(krugerCAMS$x,krugerCAMS$y,pch="+",col="red",ylim=c(min(popn$y)-dy,max(popn$y)+dy),
     xlim=range(mask$x),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(mask,add=TRUE,dots=FALSE,meshcol="gray",col=gray(level=0.99))
points(krugerCAMS$x,krugerCAMS$y,pch="+",col="red")
plot(attributes(Dfit)$polygon,add=TRUE)
dev.off()


# plots of cameras, etc:
pdf("krugerCAMS.pdf")
par(mar=c(0,0,2,0))
plot(krugerCAMS,border=0,gridlines=FALSE,detpar=list(cex=1.5))
dev.off()

pdf("krugerCH.pdf")
par(mar=c(0,0,2,0))
plot(krugerCH,border=0,gridlines=FALSE,detpar=list(cex=1.5),cappar=list(cex=4),labpar=list(cex=1.25,col="white"),tracks=TRUE,ncap=TRUE)
dev.off()


pdf("krugerMASK.pdf")
par(mar=c(0,0,0,0))
plot(krugerMASK,border=0,dots=FALSE,meshcol="gray",col=gray(level=0.99))
plot(krugerCAMS,add=TRUE)
dev.off()


# Basic fitting code
detectfn="HHN"
model=list(D~habitat.cov,lambda0~1,sigma~1)
fit=secr.fit(krugerCH,model=model,mask=krugerMASK,detectfn=detectfn)
fit
model0=list(D~1,lambda0~1,sigma~1)
fit0=secr.fit(krugerCH,model=model0,mask=krugerMASK,detectfn=detectfn)
fit0
AIC(fit,fit0)

# plot of densities etc:
xlab="Easting";ylab="Northing"
pdf("krugerHimage.pdf")
habcov=secr.surface(mix.fit,what="habitat.cov",xlab=xlab,ylab=ylab,col=terrain.colors(64),main="Habitat")
dev.off()
pdf("krugerHpersp.pdf")
mar=c(0,0,0,0)
persp(habcov,phi=65,theta=3,expand=0.75,r=0,d=2,xlab=xlab,ylab=ylab,zlab="Habitat")
dev.off()

pdf("krugerDimage.pdf")
Dsurf=secr.surface(fit,xlab=xlab,ylab=ylab,main="Density")
dev.off()
pdf("krugerDpersp.pdf")
persp(Dsurf,phi=65,theta=3,expand=0.75,r=1,d=2,xlab=xlab,ylab=ylab,zlab="Density")
dev.off()

# Confidence intervals, etc.
D=secr.surface(fit,scale=10^4,plot=FALSE)
D.cv=secr.surface(fit,"D.cv",scale=10^4,plot=FALSE)
D.lcl=secr.surface(fit,"D.lcl",scale=10^4,plot=FALSE)
D.ucl=secr.surface(fit,"D.ucl",scale=10^4,plot=FALSE)
zlim=range(na.omit(c(D$z,D.cv$z,D.lcl$z,D.ucl$z)))
pdf("krugerDCIs.pdf",h=10,w=10)
par(mfrow=c(2,2))
D=secr.surface(fit,scale=10^4,zlim=zlim,xlab=xlab,ylab=ylab,main=expression(hat(D)))
D.cv=secr.surface(fit,"D.cv",scale=10^4,xlab=xlab,ylab=ylab,main=expression(hat(cv)(hat(D))))
D.lcl=secr.surface(fit,"D.lcl",scale=10^4,zlim=zlim,xlab=xlab,ylab=ylab,main=expression(hat(D)[lower]))
D.ucl=secr.surface(fit,"D.ucl",scale=10^4,zlim=zlim,xlab=xlab,ylab=ylab,main=expression(hat(D)[uppper]))
dev.off()

# Density plot with selected region
pdf("krugerDimageBox.pdf")
Dsurf=secr.surface(fit,xlab=xlab,ylab=ylab,main="Density")
rect(4080000,-3328069,4110000,-3290000,lwd=4,col="yellow",density=0)
dev.off()

# Abundance and mask subsetting:
region.N(fit)
keep=(krugerMASK$x<4110000 & krugerMASK$y< -3290000)
smallmask=subset(krugerMASK, subset=keep)
region.N(fit,region=smallmask)

pdf("krugerDsubimage.pdf")
Dsubsurf=secr.surface(fit,region=smallmask,xlab=xlab,ylab=ylab,main="Density")
dev.off()
pdf("krugerDsubpersp.pdf")
persp(Dsubsurf,phi=65,theta=3,expand=0.75,r=1,d=2,xlab=xlab,ylab=ylab,zlab="Density")
dev.off()

# number captures in subarea
camkeep=(krugerCAMS$x<4110000 & krugerCAMS$y< -3290000)
subCH=krugerCH[,,camkeep]
length(which(apply(subCH,1,sum)>0))




# Splitting plots
# --------------
pdf("split.pdf",h=4,w=20)
par(mar=c(2,2,2,2),mfrow=c(1,5))
Dest=secr.surface(mix.fit,main=expression(D(x)),xaxt="n",yaxt="n",bty="n",xlab="",ylab="",key=FALSE)
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))

plot(0,0,pch="=",cex=10,xlim=c(-1,1),ylim=c(-1,1),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")

dethaz=data.frame(x=mask$x,y=mask$y,z=H)

p=dethaz
p$z=1-exp(-dethaz$z)

Dp=p
Dp$z=Dp$z*Ddf$z
#pdf("seen.pdf",h=6,w=6)
pest=prep4image(Dp,main=expression(D(x)*p(x)),xaxt="n",yaxt="n",bty="n",xlab="",ylab="",key=FALSE)
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))

plot(0,0,pch="+",cex=10,xlim=c(-1,1),ylim=c(-1,1),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")

p0=dethaz
p0$z=exp(-dethaz$z)
Dp0=p0
Dp0$z=Dp0$z*Ddf$z
#pdf("unseen.pdf",h=6,w=6)
pest=prep4image(Dp0,main=expression(D(x)*(1-p(x))),xaxt="n",yaxt="n",bty="n",xlab="",ylab="",key=FALSE)
plot(attributes(Dfit)$polygon,add=TRUE)
plot(krugerCAMS,add=TRUE,detpar=list(col="white"))
dev.off()





# After Seattle Workshop:
# =======================
load("survey-objects.RData")
kruger.cams=addCovariates(traps(kruger.capt),kruger.mask)
names(covariates(kruger.cams))
covariates(kruger.cams)$habitat.cov
traps(kruger.capt)=kruger.cams
traps(kruger.capt.bin)=kruger.cams
save(kruger.cams,kruger.capt,kruger.capt.bin,kruger.mask,file="survey-objects.RData")

detfn="HHN"
er.model=list(D~1,lambda0~1,sigma~1)
er.fit=secr.fit(kruger.capt,model=er.model,mask=kruger.mask,detectfn=detfn)
er.fit

detfn="HN"
df.model=list(D~1,g0~1,sigma~1)
df.fit=secr.fit(kruger.capt,model=df.model,mask=kruger.mask,detectfn=detfn)
df.fit

plot(er.fit,xval=seq(0,20000,length=200))
plot(df.fit,xval=seq(0,20000,length=200),add=TRUE,col="blue")
esa(er.fit)
esa(df.fit)
