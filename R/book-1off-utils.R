#' @title Plots overall detection probability surface p..(s) on mesh
#'
#' @description
#'  Plots overall detection probability surface p..(s) on mesh. Currenly only plots for a 
#'  single individual and only multi-catch traps and count detectors with a single session 
#'  implemented.
#'  
#'  Requires function calc.p.logPis() and plotcovariate().
#'  
#' @param capthist Object of class  '\code{capthist}' from package \code{secr}. Currently only 
#' multi-catch traps and count detectors with a single session implemented.
#' @param mesh Object of class  '\code{capthist}' from package \code{secr}.
#' @param pars Vector of parameters, in the format required by function \code{scr.negloglik}.
#' @param dist K by M Matrix distances from each of the K detectors to each of the M points on 
#' the mesh.
#' @param contour TRUE if you want a contour plotted on top of the image plot.
#' @param key TRUE if you want a key (legend) added showing probability levels.
#' @param ... Other arguments for \code{plotcovariate}.
#' 
#' @export plot.p..
plot.p..=function(capthist, mesh, pars, dist=NULL, contour=TRUE, key=TRUE, ...) {
  M=dim(mesh)[1] # number of mesh points
  a = rep(attributes(mesh)$area,M) # area of each mesh cell (vector of length M)
  
  D = rep(exp(pars["D"]),M) # Density at each mesh point (vector of length M)
  
  # Calculate p..(s) and log(P_i(s)) at each mesh point
  scr = calc.p.logPis(capthist, mesh, pars, dist=dist)
  
  covariates(mesh)$p..s = p..s
  plotcovariate(mesh,covariate="p..s",main="p..s",contour=contour, key=key, ...)
}



#' @title Image plot contour of activity centre location.
#'
#' @description
#'  Image plot contour of activity centre location. Currenly only plots for a single individual
#'  and only multi-catch traps and count detectors with a single session implemented.
#'  
#'  Requires function calc.p.logPis() and plotcovariate().
#'  
#' @param i A scalar or vector indices of individuals to plot - corresponding to the order they 
#' appear in the capture history.
#' @param capthist Object of class  '\code{capthist}' from package \code{secr}. Currently only 
#' multi-catch traps and count detectors with a single session implemented.
#' @param mesh Object of class  '\code{capthist}' from package \code{secr}.
#' @param pars Vector of parameters, in the format required by function \code{scr.negloglik}.
#' @param dist K by M Matrix distances from each of the K detectors to each of the M points on 
#' the mesh.
#' @param contour TRUE if you want a contour plotted on top of the image plot.
#' @param key TRUE if you want a key (legend) added showing probability levels.
#' @param ... Other arguments for \code{plotcovariate}.
#' 
#' @export plot.Pi
plot.Pi=function(i,capthist, mesh, pars, dist=NULL, contour=TRUE, key=TRUE, ...) {
  M=dim(mesh)[1] # number of mesh points
  a = rep(attributes(mesh)$area,M) # area of each mesh cell (vector of length M)
  
  dets=traps(capthist) # detectors
  
  D = rep(exp(pars["D"]),M) # Density at each mesh point (vector of length M)
  
  # Calculate p..(s) and log(P_i(s)) at each mesh point
  scr = calc.p.logPis(capthist, mesh, pars, dist=dist)
  
  covariates(mesh)$Pi.s = scr$Pi.s[i,]/sum(scr$Pi.s[i,])
  plotcovariate(mesh,covariate="Pi.s",main=expression(f(s[i]*"|"*capthist[i])),contour=contour, key=key, ...)
  plot(dets,add=TRUE,detpar=list(pch=19,col="black"))
  if(dim(ch)[2]==1) freq = apply(apply(ch,c(2,3),sum),2,sum)
  else freq=apply(capthist[i,,],2,sum)
  detected=which(freq>0)
  text(dets$x[detected],dets$y[detected],labels=freq[detected],cex=0.75,col="white")
}



#' @title Encounter rate fit plot
#'
#' @description
#'  Plots mean encounter rate by distance interval, with estimated encounter rate function
#'  overlaid. NB: Distances calculated from estimated activity centre mode. Don't over-interpret
#'  the plot - the activity centre locations may have very high uncertainty attached to them, 
#'  and hence the distances may be very uncertain!
#'  
#'  The function only works with traps of class '\code{count}'.
#'  
#'  Uses function '\code{histline}'.
#'  
#' @param scrfit Object of class  '\code{secr}' from package \code{secr}. The function only 
#' works with traps of class '\code{count}'.
#' @param dmax Maximum distance to plot.
#' @param binwidth width of each distance interval to plot.
#' @param lwd Line width for encounter rate function.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param hcol Colour of histogram.
#' @param hlineonly If TRUE, uses \code{\link{lines}} to draw histogram lines on current plot; 
#' else uses \code{\link{plot}} to draw lines on new plot.
#' @param houtline If TRUE, draws only the outline (profile) of the histogram; else draws each 
#' complete bar.
#' @param hfill TRUE, fills bards with \code{hcol}, else leaves them empty.
#' 
#' @export plot.er.fit
plot.er.fit=function(scrfit,dmax=NULL,binwidth=NULL,lwd=2,xlab="x",ylab="y",
                     hcol="red",hlineonly=FALSE,houtline=FALSE,hfill=TRUE) {
  
  if (detector(traps(scrfit$capthist))!="count") stop("This function only works with count detectors.")
  
  ch = scrfit$capthist
  n=summary(ch)$counts["n","Total"] # number of individuals
  loc=data.frame(x=rep(0,n),y=rep(0,n)) # data frame for locations of individuals
  for(i in 1:n) loc[i,] = fxi.mode(scrfit,i=i) # add estimated modes of activity centess to data frame
  dists = distances(loc,traps(ch)) # distances between indiv modes and detectors
  if(is.null(dmax)) dmax = max(dists)
  if(is.null(binwidth)) binwidth = dmax/10
  capts = apply(ch,c(1,3),sum) # sum over occasions
  dvec = as.vector(dists) # turn into a vector of distances
  cvec = as.vector(capts) # turn into a vector of counts
  keep = dvec<=dmax # look only within given max dist
  
  breaks=seq(0,dmax,binwidth)
  nb=length(breaks)-1
  mean.n=rep(0,nb)
  for(i in 1:nb) {
    keep=which(breaks[i]<=dvec & dvec<breaks[i+1])
    if(length(keep)>0) mean.n[i]=mean(cvec[keep])
  }
  
  # gymnastics to avoid plotting while getting plot values in pl
  ff <- tempfile()
  png(filename=ff)
  pl=plot(scrfit,xval=seq(0,dmax,length=200),lwd=2)
  dev.off()
  unlink(ff)
  
  # now do the plot we want
  histline(mean.n,breaks,xlab=xlab,ylab=ylab,col=hcol,fill=hfill,lineonly=hlineonly,outline=houtline,
           ylim=c(0,max(c(mean.n,pl$y))))
  plot(scrfit,xval=seq(0,dmax,length=200),lwd=lwd,add=TRUE)
}





rotatexy = function (df, degrees) 
{
  rotatefn <- function(xy) {
    x <- xy[1] - centrexy[1]
    y <- xy[2] - centrexy[2]
    x2 <- x * cos(theta) + y * sin(theta) + centrexy[1]
    y2 <- -x * sin(theta) + y * cos(theta) + centrexy[2]
    c(x2, y2)
  }
  centrexy <- c(mean(df$x), mean(df$y))
  theta <- 2 * pi * degrees/360
  df2 <- data.frame(t(apply(df, 1, rotatefn)))
  names(df2) <- c("x", "y")
  
  return(df2)
}


prep4image = function(dat, plot = TRUE, contour = TRUE, key = TRUE, ...){
  
  # convert factor dat$z to integer:
  zfactor=FALSE
  if(is.factor(dat$z)) {
    zfactor=TRUE
    fac=dat$z
    facnames=levels(fac)
    nlevels=length(facnames)
    dat$z=rep(0,length(fac))
    got=rep(FALSE,nlevels)
    for(i in 1:nlevels){
      j=which(fac==facnames[i])
      if(length(j)>0) got[i]=TRUE
      dat$z[j]=i
    }
    facnames=facnames[got] # remove factor names not in mask
  }
  dat = as.matrix(dat)
  
  x = sort(unique(dat[,"x"]))
  y = sort(unique(dat[,"y"]))
  
  z = matrix(NA, nrow = length(x), ncol = length(y))
  
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      m = which(dat[,"x"] == x[i] & dat[,"y"] == y[j]) ; m
      z[i,j] = if(length(m) == 0) NA else dat[,"z"][m]
    }
  }
  
  if(plot){
    if(key){
      image.plot(x, y, z, ...)
    }else{
      image(x, y, z, ...)
    }
    if(contour) contour(x, y, z, add = TRUE)
  }
  
  outlist=list(x = x, y = y, z = z)
  if(zfactor) attributes(outlist)$facnames=facnames
  
  invisible(outlist)
  
}




plot.eland.detprob = function(fit,mask=NULL,occ="all",...) {
  cams = traps(fit$capthist)
  ch = fit$capthist
  if(is.null(mask)) mask = fit$mask
  
  betas = coefficients(fit) # beta paameters
  
  # Make linear predictor and calculate g0
  g0s = which(substr(row.names(betas),1,2)=="g0")
  beta.g0 = betas$beta[g0s]
  X.g0 = model.matrix(fit$model$g0,data=elandmask)
  masklp.g0 = X.g0%*%beta.g0
  maskg0.hat=invlogit(masklp.g0)[,1]
  
  # Make linear predictor and calculate sigma
  sigmas = which(substr(row.names(betas),1,5)=="sigma")
  beta.sigma = betas$beta[sigmas]
  X.sigma = model.matrix(fit$model$sigma,data=elandmask)
  masklp.sigma = X.sigma%*%beta.sigma
  masksigma.hat=exp(masklp.sigma)[,1]
  
  #function to calculate Halfnormal dectection function:
  pdet = function(d,g0,sigma,nocc) {
    npar = length(g0)
    ncam = dim(d)[1]
    p = matrix(rep(NA,npar*ncam),nrow=ncam)
    for(i in 1:ncam) p[i,] = 1 - (1-g0*exp(-d[i,]^2/(2*sigma^2)))^nocc
    return(p)
  }
  # Function to combine across independen detect probs
  combdet = function(ps) 1 - prod(1-ps)
  
  ncam = dim(cams)[1]
  nmask = dim(elandmask)[1]
  dists = matrix(rep(NA,nmask*ncam),nrow=ncam)
  for(i in 1:ncam) {
    dists[i,] = distancetotrap(mask,cams[i,])
  }
  nocc = dim(ch)[2]
  #  pest = rbind(pest,pestot)
  
  if(occ=="all") {
    pest = pdet(dists,maskg0.hat,masksigma.hat,nocc=nocc)
    pestot = apply(pest,2,combdet)
  } else if(occ=="single") {
    pest = pdet(dists,maskg0.hat,masksigma.hat,nocc=1)
    pestot = apply(pest,2,combdet)
  } else stop("Invalid occ passed")
  
  covariates(mask)$p = pestot
  plotcovariate(mask,"p",...)
  
}


