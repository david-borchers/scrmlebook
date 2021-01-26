library(secr)
library(scrmlebook)
library(raster)
library(gdistance)
library(igraph)
library(sp)

# OLD Figure (see below for new figure)
# =====================================
#ny=4; nx=4 # dimensions of mask
## set costs (NA is "hole" - nothing there & can't go there):
#z=c(log(2),log(2),log(2),log(2),NA,log(2),log(2),NA,0,NA,log(2),0,0,0,0,0) 
#z=c(1,1,1,1,NA,1,1,NA,0,NA,1,0,0,0,0,0) 
#beta.0 = 0
#alpha.z = log(2)
#lp = beta.0 + alpha.z*z
#rmesh=data.frame(x=rep(1:nx,ny),y=rep(1:ny,rep(nx,ny)),lp=lp) # make data frame with coords and costs
#
#rmask=read.mask(data=rmesh,columns="lp")  # make mask with linear predictor of conductance as covariate
#
#pdf("./keepfigure/igraphplot.pdf",h=6,w=12)
#par(mfrow=c(1,2))
#plotcovariate(rmask,covariate="lp",col=c("lightgray","gray"),what="image") # look at the covariate surface
#text(rmesh$x,rmesh$y,labels=c(13:16,9:12,5:8,1:4),cex=1.25)
#
#layout=cbind(c(1:4,1,3,4,2,3,1:4),c(4,4,4,4,3,3,3,2,2,1,1,1,1))
#vcol=c(rep("lightgray",5),"gray","lightgray",rep("gray",6))
#ig=make_igraph(rmask,"lp",cfun=mean,conductance="FALSE")
#layout=layout[order(as.numeric(names(V(ig)))),]
#vcol=vcol[order(as.numeric(names(V(ig))))]
#plot(ig, edge.label=signif(E(ig)$weight, 2), edge.label.cex=0.9,layout=layout,vertex.label.cex=1,vertex.size=25,vertex.color=vcol)
#dev.off()
#
#pdf("./keepfigure/lcpath.pdf",h=6,w=6)
#from=c(4,1); to=c(2,4) # corrdinates of start and end points (Note: need not be exactly on mask coords)
#dist = plot_lcpath(from,to,rmask,lwd=2,linecol="black",col=c("lightgray","gray"),key=FALSE) # calculate and plot the least cost path, returning cost
#dev.off()

# Example with very simple 4x4 grid (and buffer for commuteDistance)
# ------------------------------------------------------------------
values=matrix(c(NA,NA,NA,NA,NA,NA,
                NA,1,1,1,1,NA,
                NA,1,3,3,3,NA,
                NA,1,3,1,NA,NA,
                NA,1,3,NA,1,NA,
                NA,NA,NA,NA,NA,1),nrow=6,byrow=TRUE)
rx = raster(values, xmn=-0.5, xmx=5.5, ymn=-0.5, ymx=5.5)
temp = Polygons(list(poly1=Polygon(data.frame(x=c(0.5,0.5,4.5,4.5),y=c(0.5,4.5,4.5,0.5)))),ID=1)
boundary = SpatialPolygons(list(temp))
plot(rx,col=c("lightgray","gray"))
plot(boundary,add=TRUE)
region = mask(rx,boundary)
xy = coordinates(region)[!is.na(values(region)),]
invalues = values(region)[!is.na(values(region))]
region.sp = SpatialPixelsDataFrame(xy,data.frame(values=invalues))

rmask = read.mask(data=coordinates(rx))
covariates(rmask)$lp = as.vector(values)
# reduce to just the area within the boundary:
inside = (0<rmask$x & rmask$x<5) & (0<rmask$y & rmask$y<5)
rmask = subset(rmask,inside)

pdf("./keepfigure/igraphplot.pdf",h=6,w=12)
par(mfrow=c(1,2))
plot(region.sp,col=c("lightgray","gray"),what="image")
text(rmask$x,rmask$y,labels=1:16,cex=1.25)

ig=make_igraph(rmask,"lp",cfun=mean,conductance="FALSE")
layout=cbind(c(1:4,1:4,1:3,1:2,4),c(4,4,4,4,3,3,3,3,2,2,2,1,1,1))
layout=layout[order(as.numeric(names(V(ig)))),]
vcol=c(rep("lightgray",5),rep("darkgray",3),rep(c("lightgray","darkgray","lightgray"),2))
vcol=vcol[order(as.numeric(names(V(ig))))]
plot(ig, edge.label=signif(E(ig)$weight, 2), edge.label.cex=0.9,layout=layout,vertex.label.cex=1,vertex.size=25,vertex.color=vcol)
dev.off()


tr3<-transition(rx,transitionFunction=function(x) 1/mean(x),directions=8,symm=symm)
costDist = costDistance(tr3,xy)
rSPDist = rSPDistance(tr3,xy,xy,theta=0)
rSP5Dist = rSPDistance(tr3,xy,xy,theta=5)
commDist = commuteDistance(tr3,xy)

xfrom = 4
yfrom = 4
from = which(xy[,1]==xfrom & xy[,2]==yfrom)
from.Edist = sqrt((xy[from,1]-xy[,1])^2 + (xy[from,2]-xy[,2])^2)
Edist = SpatialPixelsDataFrame(xy,data.frame(dist=from.Edist/max(from.Edist)))
from.costdist = as.matrix(costDist)[,from]
costdist = SpatialPixelsDataFrame(xy,data.frame(dist=from.costdist/max(from.costdist)))
from.rSPdist = as.matrix(rSPDist)[,from]
rSPdist = SpatialPixelsDataFrame(xy,data.frame(dist=from.rSPdist/max(from.rSPdist)))
from.rSP5dist = as.matrix(rSP5Dist)[,from]
rSP5dist = SpatialPixelsDataFrame(xy,data.frame(dist=from.rSP5dist/max(from.rSP5dist)))
from.commdist = as.matrix(commDist)[,from]
commdist = SpatialPixelsDataFrame(xy,data.frame(dist=from.commdist/max(from.commdist)))

to1 = c(1,1)
to2 = c(4,1)
lcp1 = shortestPath(tr3,c(xfrom,yfrom),c(1,1),"SpatialLines")
lcp2 = shortestPath(tr3,c(xfrom,yfrom),c(4,1),"SpatialLines")

pdf("./keepfigure/distancemetrics.pdf",h=6,w=6)
layout(matrix(c(1,2,3,4),nrow=2,byrow=TRUE))
plot(region.sp,main="Resistance",what="image",col=c("lightgray","darkgray"))
plot(boundary,add=TRUE)
plot(Edist,main="Euclidian Distance",what="image",col=gray.colors(40,start=0.9,end=0.1))
plot(boundary,add=TRUE)
points(xy[from,1],xy[from,2],col="green",pch=19)
points(to1[1],to1[2],col="red",pch=19)
segments(xy[from,1],xy[from,2],to1[1],to1[2],col="white")
#plot(lcp2,add=TRUE,col="white")
plot(costdist,main="Least-Cost Distance",what="image",col=gray.colors(40,start=0.9,end=0.1))
plot(boundary,add=TRUE)
points(xy[from,1],xy[from,2],col="green",pch=19)
points(xy[from,1],xy[from,2],col="green",pch=19)
points(to1[1],to1[2],col="red",pch=19)
#points(to2[1],to2[2],col="red",pch=19)
plot(lcp1,add=TRUE,col="white")
#plot(lcp2,add=TRUE,col="white")
#plot(rSP5dist,main="Directed Resistance Distance (5)",what="image",col=parula(12))
#plot(boundary,add=TRUE)
#plot(rSPdist,main="Directed Resistance Distance",what="image",col=gray.colors(40,start=0.9,end=0.1))
#plot(boundary,add=TRUE)
#points(xy[from,1],xy[from,2],col="green",pch=19)
plot(commdist,main="Resistance Distance",what="image",col=gray.colors(40,start=0.9,end=0.1))
plot(boundary,add=TRUE)
points(xy[from,1],xy[from,2],col="green",pch=19)
#plot(commdist,main="Resistance Distance",what="scale",col=gray.colors(40,start=0.9,end=0.1),scale.shrink=0.1)
dev.off()

