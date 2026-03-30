#load("/Volumes/1TB/AIMSbackup/Meetings & Talks/AIMS Feb 2016/data/bcams.rda")
#load("/Volumes/1TB/AIMSbackup/Meetings & Talks/AIMS Feb 2016/data/bchist.rda")
#load("/Volumes/1TB/AIMSbackup/Meetings & Talks/AIMS Feb 2016/data/bmask.rda")
##setwd("/Users/dlb/MeetingsTalks/ISEC_Seattle_2016/Workshop/")
#save(list=c("bcams","bchist","bmask"),file="/Users/dlb/MeetingsTalks/ISEC_Seattle_2016/Workshop/bcams_bchist_bmask.rda")
load("bcams_bchist_bmask.rda")

library(secr)
library(mgcv)

standardize=function(x,centre=mean(x),scale=sd(x)) (x-centre)/scale

plot(bchist,border=0,tracks=TRUE)
names(covariates(bmask))
plot(bmask,covariate="elevation")
elev.mean=mean(covariates(bmask)$elevation)
elev.sd=sd(covariates(bmask)$elevation)
covariates(bmask)$elevation=standardize(covariates(bmask)$elevation)
bcams=traps(bchist)
names(covariates(bcams))
covariates(traps(bchist))$elevation=(covariates(traps(bchist))$elevation-elev.mean)/elev.sd
hist(covariates(bmask)$elevation)
points(covariates(traps(bchist))$elevation,covariates(traps(bchist))$elevation*0,pch="|",col="red")

# aside: extract some traps for use in slides
plot(bmask,covariate="elevation",legend=FALSE,border=0)
text(bcams$x,bcams$y,as.character(1:length(bcams$x)),cex=0.5)
keep=c(20,26,35,40,42)
kcams=bcams[keep,]
xrange=range(kcams$x)
yrange=range(kcams$y)
dx=diff(xrange)
dy=diff(yrange)
xlim=c(min(kcams$x)-dx/2,max(kcams$x)+dx/2)
ylim=c(min(kcams$y)-dy/2,max(kcams$y)+dy/2)
keep=which(xlim[1]<=bmask$x & bmask$x<=xlim[2] & ylim[1]<=bmask$y & bmask$y<=ylim[2])
kmask=bmask[keep,]
class(kmask)=class(bmask)
covariates(kmask)=covariates(bmask)[keep,]
plot(kmask,covariate="elevation",legend=FALSE,border=0)
points(kcams$x,kcams$y)

hframe=data.frame(x=kmask$x,y=kmask$y,z=covariates(kmask)$elevation)
gfit=gam(z~s(x,y,k=30),family=gaussian,data=hframe)
plot(gfit)
nxy=64
xs=rep(seq(min(hframe$x),max(hframe$x),length=nxy),nxy)
ys=rep(seq(min(hframe$y),max(hframe$y),length=nxy),rep(nxy,nxy))
pdat=data.frame(x=xs,y=ys)
ele=predict(gfit,type="response",newdata=pdat)
imdat=data.frame(x=xs,y=ys,z=ele)


pdf("locnhetero.pdf",h=3,w=5.5)
par(mar=c(0,0,0,3))
mp=prep4image(imdat,contour=FALSE,col=terrain.colors(25))
points(kcams$x,kcams$y,col="white",cex=3.5,pch=19)
points(kcams$x,kcams$y,col="red",cex=3.5,lwd=2)
text(kcams$x,kcams$y,c(4,3,5,2,1),col="black",cex=1)
xi=c(-1,-21.8)
points(xi[1],xi[2],pch=19,cex=1.5)
text(xi[1],xi[2],expression(s[1]),pos=1)
xj=c(-3.5,-16.5)
points(xj[1],xj[2],pch=19,cex=1.5)
text(xj[1],xj[2],expression(s[2]),pos=1)
dev.off()



pdf("habitat.pdf",h=5)
mp=prep4image(imdat,contour=FALSE,col=terrain.colors(25))
points(kcams$x,kcams$y,col="white",cex=3.5,pch=19)
points(kcams$x,kcams$y,col="red",cex=3.5,lwd=2)
text(kcams$x,kcams$y,c(4,3,5,2,1),col="black",cex=1)
xi=c(-1,-21.8)
points(xi[1],xi[2],pch=19,cex=2)
text(5,-21,"Activity centre",pos=3)
arrows(3,-21,xi[1]+0.5,xi[2]+0.1,length=0.1)
dev.off()

pdf("usage.pdf",h=5)
gdist=(sqrt((imdat$x-xi[1])^2+(imdat$y-xi[2])^2)/(max(ele)-ele+0.5))
use=(max(gdist)-gdist)^10
udat=data.frame(x=xs,y=ys,z=use/max(use)*12)
mp=prep4image(udat)
points(kcams$x,kcams$y,col="white",cex=3.5,pch=19)
points(kcams$x,kcams$y,col="red",cex=3.5,lwd=2)
text(kcams$x,kcams$y,c(4,2,5,3,1),col="black",cex=1)
points(xi[1],xi[2],pch=19,cex=2)
text(5,-21,"Activity centre",pos=3)
arrows(3,-21,xi[1]+0.5,xi[2]+0.1,length=0.1)
dev.off()



## Logit, cloglog, and log links
## ===================
cloglog=function(p) log(-log(1-p))
invcloglog=function(lp) 1-exp(-exp(lp))

x=seq(-10,10,length=200)
g0.1 = plogis(x/2)
g0.2 = plogis(-x/2)

lambda.1=exp(x/3)
lambda.2=exp(-x/3)

g0.1.cll=invcloglog(x/3-0.5)
g0.2.cll=invcloglog(-x/3-0.5)

pdf("./keepfigure/linkfuns.pdf",h=3,w=8)
par(mfrow=c(1,2))
plot(x,g0.1,type="l",lwd=2,xlab=expression(eta==beta[0]+beta[1]*x),ylab=expression(g[0]))
lines(x,g0.2,lty=2,lwd=2)
lines(x,g0.1.cll,lty=1,lwd=2,col="gray")
lines(x,g0.2.cll,lty=2,lwd=2,col="gray")
plot(x,lambda.1,type="l",lwd=2,xlab=expression(eta==beta[0]+beta[1]*x),ylab=expression(lambda[0]))
lines(x,lambda.2,lty=2,lwd=2)
dev.off()










# Stuff I think is unused:
sex2.fit=secr.fit(bchist,model=list(g0~h2,sigma~h2),hcov="sex",mask=bmask)
sexg0.fit=secr.fit(bchist,model=list(g0~h2),hcov="sex",mask=bmask) # Hessian does not converge
sexsig.fit=secr.fit(bchist,model=list(sigma~h2),hcov="sex",mask=bmask)

eleg0.fit=secr.fit(bchist,model=list(g0~elevation),mask=bmask)
eleg0.fit
elesigma.fit=secr.fit(bchist,model=list(sigma~elevation),mask=bmask)
elesigma.fit
ele2.fit=secr.fit(bchist,model=list(g0~elevation,sigma~elevation),mask=bmask) # Hessian does not converge
#ele2.fit

# Convergence failure for all these:
sex2ele2.fit=secr.fit(bchist,model=list(g0~elevation+h2,sigma~elevation+h2),hcov="sex",mask=bmask)
sex2eleg0.fit=secr.fit(bchist,model=list(g0~elevation+h2,sigma~h2),hcov="sex",mask=bmask)
sex2elesigma.fit=secr.fit(bchist,model=list(g0~h2,sigma~elevation+h2),hcov="sex",mask=bmask)
sexD2.fit=secr.fit(bchist,model=list(D~elevation,g0~h2,sigma~h2),hcov="sex",mask=bmask)
sigmah2.eleDsigma.fit=secr.fit(bchist,model=list(D~elevation,g0~elevation,sigma~h2),hcov="sex",mask=bmask)

AIC(eleg0.fit,elesigma.fit,ele2.fit)
AIC(sex2.fit,sexg0.fit)

# DNC:
hrsex2.fit=secr.fit(bchist,model=list(g0~h2,sigma~h2,z~1),detectfn="HR",hcov="sex",mask=bmask)
sex2a.fit=secr.fit(bchist,model=list(a0~h2,sigma~h2),hcov="sex",mask=bmask)
sexa.fit=secr.fit(bchist,model=list(a0~h2,sigma~1),hcov="sex",mask=bmask)
sexsigma.a.fit=secr.fit(bchist,model=list(a0~1,sigma~h2),hcov="sex",mask=bmask)

sex2.fit
cov2cor(sex2.fit$beta.vcv)
cov2cor(sexg0.fit$beta.vcv)




