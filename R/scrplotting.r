#' @title Plots density as function of a covariate
#' 
#' @param fit An object of class "secr".
#' @param cov A character variable of length 1 giving name of covariate to use 
#' for the x-axis. This must be the name of a covariate on the mask
#' @param mask If not NULL, is object of class "mask". If NULL, takes from fit.
#' @param per.km If TRUE, converts density to number per square km, else per hectare.
#' @param add.traps TRUE if you want the locations of traps to appear as ticks on the x-axis.
#' @param sess Number of session to use (if NULL, first session is used)
#' @param logD If TRUE, log(Density+1) is plotted instead of Density
#' @param Dlim vector containing min and max Density to plot (y-axis is restircted to this range)
#' 
#' @export
plotDline = function(fit,cov,mask=NULL,per.km=TRUE,add.traps=TRUE,sess=NULL,logD=FALSE,Dlim=NULL) {

  if(!is.null(Dlim) & (length(Dlim)!=2)) stop("Dlim must be a vector of length 2.")

    scale = 1
  unitlabel = "(/ha)"
  if(per.km) {
    scale = 100^2
    unitlabel = "(/sq. km)"
  }
  component = "D" # this was when I thought I'd use the function for things other than D!
  # Get variables in model
  vars = all.vars(fit$model[[component]]) # get names of variables in model
  if(!is.element(cov,vars)) stop(paste(cov," is not a variable in the model.",sep=""))
  if(is.null(mask)) mesh = fit$mask
  else mesh = mask
  if(is.list(fit$mask) & !is.list(mesh)) 
    stop("When fit$mask is a list, the mask that is passed must be a list of the same length")
  if(inherits(fit$mask,"list")) {
    for(i in 1:length(fit$mask)) mesh[[i]] = make.covmesh(mesh[[i]],vars,cov)
  } else mesh = make.covmesh(mesh,vars,cov)
  # and predict, with CI
  Dpred = predictDsurface(fit,mask=mesh,se.D=TRUE,cl.D=TRUE)
  
  if(inherits(fit$mask,"list")) {
    if(is.null(sess)) {
      warning("Argument 'sess'' is NULL; first session being plotted.")
      sess = 1
    } 
    if(sess>length(fit$mask)) stop("Argument 'sess' > number sessions.")
    Dpred = Dpred[[sess]]
    traps = traps(fit$capthist[[sess]])
    mesh = mesh[[sess]]
  } else {
    traps = traps(fit$capthist)
  }
  
  Dcov = covariates(Dpred)
  Dframe = data.frame(x=Dcov[,cov],Dhat=Dcov$D.0,D.lcl=Dcov$lcl.0,D.ucl=Dcov$ucl.0)
  Dframe[,2:4] = Dframe[,2:4]*scale
  
  if(logD) {
    if(is.null(Dlim)) {
      Dlim = log(range(Dframe$D.lcl,Dframe$D.ucl)+1)
      ylim = range(0,Dlim)
    } else ylim = log(Dlim+1)
    plot(Dframe$x,log(Dframe$Dhat+1),type="l",ylim=ylim,main="",xlab=cov,ylab=paste("log(Density+1)",unitlabel))
    CIdf = data.frame(Dframe[,1],log(Dframe[,3]+1),log(Dframe[4]+1))
    shadedCI(CIdf,col=adjustcolor("gray", alpha.f=0.4),border=NA)
  } else {
    if(is.null(Dlim)) {
      Dlim = range(Dframe$D.lcl,Dframe$D.ucl)
      ylim = range(0,Dlim)
    } else ylim = Dlim
    plot(Dframe$x,Dframe$Dhat,type="l",ylim=ylim,main="",xlab=cov,ylab=paste("Density",unitlabel))
    shadedCI(Dframe[,c(1,3,4)],col=adjustcolor("gray", alpha.f=0.4),border=NA)
  }
  
  trapcov = NULL
  
  if(add.traps) {
    # get indices of mesh points closest to traps:
    if(inherits(fit$mask,"list"))
      pseudotraps = read.traps(data=data.frame(x=fit$mask[[sess]]$x,y=fit$mask[[sess]]$y))
    else 
      pseudotraps = read.traps(data=data.frame(x=fit$mask$x,y=fit$mask$y))
    covind = nearesttrap(data.frame(x=traps$x,y=traps$y),pseudotraps)
    # get corresponding covariate values
    if(cov=="x") trapcov = mesh$x[covind] # mesh covairate nearest to each trap
    else if(cov=="y") trapcov = mesh$y[covind] # mesh covairate nearest to each trap
    else trapcov = covariates(mesh)[covind,cov] # mesh covairate nearest to each trap
    points(trapcov,rep(0,length(trapcov)),pch="|")
  }
  
  invisible(list(Density=Dframe,traps=trapcov))
  
}

# Function to make mesh with covariates in density model (vars) equal to their means, 
# except for the specified covariate (cov), which is set to equally-spaced points
# spanning its range.
# In the case of factors, the first value is used in place of the mean value, and
# if cov is a factor, the function stops with an error message. It is designed to
# deal only with numeric variables.
make.covmesh = function(mesh,vars,cov) {
  M = dim(mesh)[1]
  newmesh = mesh # output mesh
  covs = covariates(mesh) 
  # get model variables that are covariates on the mesh
  varcols = which(is.element(names(covs),vars)) 
  # keep only covariates that are in model (i.e. in vars)
  #  covs = covs[,varcols,drop=FALSE]
  covnames = names(covs) # get names of selected covariates
  ncovs = length(varcols)
  if(!is.element(cov,covnames) & !is.element(cov,"x") & !is.element(cov,"y")) 
    stop(paste(cov," is not a variable on the mesh.",sep=""))
  # Deal with mesh covariates
  if(length(varcols)>0) {
    #    covmeans = covs[1,]
    #    names(covmeans) = covnames
    for(i in 1:ncovs) {
      if(is.factor(covs[,varcols[i]]) | is.character(covs[,varcols[i]])) {
        warning(paste("Arbitrarily chosen first level (",covs[1,varcols[i]],") of factor ",covnames[varcols[i]]," for plotting.",sep=""))
        covariates(newmesh)[,varcols[i]] = as.factor(covs[,varcols[i]])
        covariates(newmesh)[,varcols[i]] = covariates(newmesh)[1,varcols[i]]
        #        levels(covariates(newmesh)[,varcols[i]]) = levels(covs[,varcols[i]])
      } else {
        covariates(newmesh)[,varcols[i]] = mean(covs[,varcols[i]])
      }
    }
  }
  # x and y are not in mesh covariates, so need to handle differently:
  if(cov=="x") { 
    covrange = range(newmesh$x)
    newmesh$x = seq(covrange[1],covrange[2],length=M) # insert range of x's
    newmesh$y = rep(mean(newmesh$y),M) # insert mean values of y 
  } else if(cov=="y") { 
    covrange = range(newmesh$y)
    newmesh$y = seq(covrange[1],covrange[2],length=M) # insert range of y's
    newmesh$x = rep(mean(newmesh$x),M) # insert mean values of x 
  } else { # here if cov is a covariate on the mesh
    if(is.factor(covs[,cov]) | is.character(covs[,cov])) stop("Not yet programmed plotting for factors.")
    covrange = range(covs[,cov])
    # set to ordered values from min to max of cov
    covariates(newmesh)[,which(names(covariates(newmesh))==cov)] = seq(covrange[1],covrange[2],length=M)
    # insert mean values of x and y
    newmesh$x = rep(mean(newmesh$x),M)
    newmesh$y = rep(mean(newmesh$y),M)
  }
  
  return(newmesh)
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
plot.er.fit=function(scrfit,dmax,binwidth,lwd=2,xlab="x",ylab="y",
                     hcol="red",hlineonly=FALSE,houtline=FALSE,hfill=TRUE) {
  
  if (detector(traps(scrfit$capthist))!="count") stop("This function only works with count detectors.")
  
  ch = scrfit$capthist
  n=summary(ch)$counts["n","Total"] # number of individuals
  loc=data.frame(x=rep(0,n),y=rep(0,n)) # data frame for locations of individuals
  for(i in 1:n) loc[i,] = fxi.mode(er.fit,i=i) # add estimated modes of invivs to data frame
  dists = distances(loc,kruger.cams) # distances between indiv modes and detectors
  capts = apply(kruger.capt,c(1,3),sum) # sum over occasions
  dvec = as.vector(dists) # turn into a vector of distances
  cvec = as.vector(capts) # turn into a vector of counts
  keep = dvec<=dmax # look only within given max dist
  
  breaks=seq(0,dmax,binwidth)
  nb=length(breaks)-1
  mean.n=rep(0,nb)
  for(i in 1:nb) {
    keep=which(breaks[i]<=dvec & dvec<breaks[i+1])
    mean.n[i]=mean(cvec[keep])
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


#' @title Adds covariate to mask from raster object
#' @param ras is an object of class `RasterLayer'
#' @param mask is an object of class `mask'
#' @param maskCRS is an object of class "CRS" giving the coordinate reference 
#' system of the mask coordinates
#' @param covname optional argument specifying what the covariate added to the
#' mask is to be called.
#' @export
addRaster2Mask2 = function(ras,mask,maskCRS,covname=NULL) {
  # To transform mask coordinates onto ras CRS, need to make SpatialPoints object with mask CRS
  maskspdf = SpatialPoints(data.frame(x=mask$x,y=mask$y),proj4string=maskCRS)
  # Then transform onto ras CRS
  tfmaskspdf = spTransform(maskspdf,CRS(proj4string(ras)))
  cov = extract(ras,tfmaskspdf) # exctact covariate values at maaks points
  covariates(mask) = cbind(covariates(mask),cov)
  names(covariates(mask))[length(names(covariates(mask)))] = covname[1]
  return(mask)
}


#' @title Plots mask covariate using package sp functions
#' @param mask is an object of class `mask'
#' @param covariate is a character variable with the name of one of the covariates in mask (the one to plot)
#' @param ... other arguments to be passed to \code{sp} \code{plot} function.
#' @export
splotcovariate = function(mask,covariate,...) {
  if(!inherits(mask,"mask")) stop("mask must be of class `mask`.")
  require(secr)
  require(sp)
  require(pals)
  cnum = which(names(covariates(mask))==covariate)
  if(length(cnum)==0) stop("There is no such covariate on the mask.")
  if(length(cnum)>1) warning("There is more than one covariate of that name the mask; using the first.")
  if(is.character(covariates(mask)[,cnum[1]])) covariates(mask)[,cnum[1]] = as.factor(covariates(mask)[,cnum[1]])
  spdf = SpatialPixelsDataFrame(as.matrix(mask),data=covariates(mask))
  sp::plot(spdf[cnum[1]],...)
}

#' @title Plots image (and optionally contours) of mask covariate value
#' @param mask is an object of class `mask'
#' @param covariate is a character variable with the name of one of the covariates in mask (the one to plot)
#' @param contour is a logical, TRUE if want contour plots on image
#' @param ... other arguments to be passed to \code{prep4image}
#' @export
plotcovariate=function(mask, covariate, ...) {
  cnum=which(names(covariates(mask))==covariate)
  if(is.null(cnum)) stop(paste("No covariate(s) called",covariate))
  if(length(cnum)>1) warning("Can only plot one covariate at a time. First covariate being plotted.")
  dat=data.frame(x=mask$x,y=mask$y,z=covariates(mask)[[cnum]])
  prep4image(dat,...)
}


#' @title Plots image (and optionally contours) of density from Dsurface object
#' @param Dsurface is an object of class `Dsurface'
#' @param covariate The covariate to plot (density by default)
#' @param addpoly If TRUE, adds bounding polygon
#' @param contour is a logical, TRUE if want contour plots on image
#' @param scale Factor to multiply surface by for plotting
#' @param ... other arguments to be passed to \code{prep4image}
#' @export
plot.Dsurface=function(Dsurface,covariate="D.0",scale=1,contour=TRUE,addpoly=TRUE, ...) {
  dat=data.frame(x=Dsurface$x,y=Dsurface$y,z=covariates(Dsurface)[,covariate])
  plotdat=prep4image(dat,contour=contour,...)
  
  if(!is.null(attr(Dsurface,"poly.habitat"))) {
    if(attr(Dsurface,"poly.habitat") & addpoly) {
      sp:::plot(attributes(Dsurface)$polygon,add=TRUE)
    }
  }
  invisible(plotdat)
}


#' @title Plots density surface image
#' @param fit is an object of class `secr'
#' @param se.D If TRUE, returns std error of estimated density
#' @param cl.D If TRUE, returns confidence bounds for density
#' @param alpha Significance level
#' @param ... other arguments to be passed to \code{plot.SpatialPixelsDataFrame}
#' @export
splot.Dsurface=function(fit,se.D=FALSE,cl.D=FALSE,alpha=0.05, ...) {
  require(secr)
  if(!inherits(fit,"secr")) stop("fit must be an secr object.")
  Dhat = predictDsurface(fit,se.D=se.D,cl.D=cl.D,alpha=alpha,parameter="D")
  plotcovariate(Dhat,"D.0",...)
  invisible(plotdat)
}


#' @title Prepares data frame for plotting with image/contour/persp.
#'   
#' @description From an input data frame with columns x, y and z, this function 
#'   creates a list with elements x, y and z in a format suitable for passing to
#'   functions \code{\link{image}}, \code{\link{contour}} or 
#'   \code{\link{persp}}. The coordinates in \code{data} are assumed to come
#'   from a 2D grid of points.
#'   
#' @param data a data frame with columns x and y being Cartesian coordinates, 
#'   and z being the values of some variable at each coordinate.
#' @param plot if \code{TRUE} then an image plot will be drawn using 
#'   \code{\link{image.plot}}
#' @param contour if \code{TRUE} then contours will be added (only used when 
#'   \code{plot=TRUE})
#' @param key logical for whether or not to include key when \code{plot = TRUE} (\code{\link{image.plot}} is used when \code{key = TRUE}, \code{\link{image}} is used when \code{key = FALSE})
#' @param ... other arguments to pass to \code{\link{image}} or \code{\link{image.plot}} (only used 
#'   when \code{plot=TRUE})
#'   
#' @details Sorts z on values of x first, then y, then creates a matrix of 
#'   z-values from this. Returns a list with elements x (unique values of x, in 
#'   increasing order), y (unique values of y, in increasing order) and z 
#'   (matrix of z-values in appropriate order for image/contour/persp). 
#'   
#'   If the original z is a factor variabele, the z returned is a matrix of integers 
#'   between 1 and length(levels(z)) and the output list has an attributes called 
#'   ``facnames'' that is a character vector containing the levels as factor 
#'   variables, with z=1 corresponding to the first name, z=2 to the second, etc.
#' @export
prep4image = function(data, plot = TRUE, contour = FALSE, key = TRUE, ...){
  
  # convert factor data$z to integer:
  zfactor=FALSE
  if(is.factor(data$z)) {
    zfactor=TRUE
    fac=data$z
    facnames=levels(fac)
    nlevels=length(facnames)
    data$z=rep(0,length(fac))
    got=rep(FALSE,nlevels)
    for(i in 1:nlevels){
      j=which(fac==facnames[i])
      if(length(j)>0) got[i]=TRUE
      data$z[j]=i
    }
    facnames=facnames[got] # remove factor names not in mask
  }
  data = as.matrix(data)
  
  x = sort(unique(data[,"x"]))
  y = sort(unique(data[,"y"]))
  
  z = matrix(NA, nrow = length(x), ncol = length(y))
  
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      m = which(data[,"x"] == x[i] & data[,"y"] == y[j]) ; m
      z[i,j] = if(length(m) == 0) NA else data[,"z"][m]
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


#' @title add new covariate consisting of distances to points that have covariate==distance.to
#' @param mask is object of class "mask".
#' @param covariate is name of one of the elements of covariates(mask)
#' @param distance.to is the target value of covariate; distances to mask points with this value are calculated.
#' @param dname is the name you want to give to the new covariate
add.dist.to=function(mask,covariate,distance.to,dname=NULL,overwrite=FALSE){
  covs=names(covariates(mask))
  if(is.null(covs)) stop("No covariates in mask. Can't add distance to anything")
  ncov=length(covs)
  if(is.null(dname)) dname=paste("dist2",covariate,"=",distance.to,sep="")
  if(is.element(covariate,covs)) {
    if(dname %in% names(covariates(mask))) {
      if(overwrite) {
        warning(paste("Covariate called",dname,"already existed. It has been overwritten."))
        covi=which(names(covariates(mask))==dname)
      } else stop(paste("Covariate called",covariate,"already exists. Use overwrite=TRUE if you want to overwrite it."))
    } else {
      covi=ncov+1
    }
    cov=covariates(mask)[which(covs==covariate)]
    if(is.null(cov)) stop(paste("No covariate called",covariate,"in mask."))
    targets=which(cov==distance.to)
    dists=distances(mask,mask[targets,])
    d2=apply(dists,1,min)
    covariates(mask)[[covi]]=d2
    names(covariates(mask))[[covi]]=dname
    return(mask)
  } else stop(paste("No covariate called",covariate,"in mask."))
}



#' @title Plot conditional value of fiited density smooth
#'
#' @description  Calculates and plots the estimated density over the range of values of one spatial covariate, at
#' fixed values of all other spatial covariates. This function is pretty inflexible but also pretty simple. If you want
#' to do a more flexible plot, just take the few lines of code inside the function and cahnge some things.
#' 
#' @return Invisibly returns a list comprising:
#' \itemize{
#'  \item{preds}{ A data frame with these columns:
#'  \itemize{
#'  \item{x}{ Values of the chosen covariate at each of \code{nx} points spanning its range.}
#'  \item{Dhat}{ Estimated density at each value of \code{x}.}
#'  \item{se}{ Standard error of \code{Dhat}.}
#'  \item{lcl}{ Lower 95\% confidence bound for density at each value of \code{x}.}
#'  \item{ucl}{ Upper 95\% confidence bound for density at each value of \code{x}.}
#'  }}
#'  \item{detcovs}{ Values of the chosen covariate at each detector (in order of detectors in \code{traps(fit)})}
#'  \item{condvals}{ Named vector with values of each covariate being conditioned on.}
#'  \item{xvar}{ Name of the selected covariate.}
#'  }
#' 
#' @param fit An object of class \code{secr}.
#' @param covariate A characted giving the name of the selected covariate (which must be in \code{names(covariates(mask))}; see below).
#' @param otherlevels A named vector with the selected level/value of each of the other covariates used in \code{fit} 
#' (defaults to their mean). (NOTE: Use of this argument has not been properly tested - things might go wrong!)
#' @param mask An object of class \code{mask}, which contains the covariates 'covariate' and those for 'otherlevels'.
#' @param se TRUE if you want standard errors of density calculated and returned.
#' @param lcl TRUE if you want the lower 95\% bound of density calculated and returned.
#' @param se TRUE if you want the upper 95\% bound of density calculated and returned.
#' @param nx Number of values of 'covariate' to use.
#' @param Dscale A multiplier to scale density: the defaule Dscale=100^2 gives density in number per square km.
#' @param xlab Label for x-axis. Defaults to the \code{covariate} argument value.
#' 
#' @seealso \code{\link{condDensity}}
#' 
#' @examples 
#' \dontrun{
#' library(scrmlebook)
#' data("leopard-sim") # simulated Kruger Park data
#' # Fit overly complex density model - just for illustration:
#' fit = secr.fit(kruger.capt,model=list(D~x+s(y,k=3)+s(habitat.cov,k=3)),kruger.mask)
#' smest = condDensityplot(fit,covariate="y") 
#' # Plot conditional density in dist.to.water dimension
#' smest = condDensityplot(fit,covariate="habitat.cov") 
#' }
#' @export
#' 
condDensityplot = function(fit,covariate,otherlevels=NULL,mask=NULL,se=TRUE,cl=TRUE,nx=200,Dscale=100^2,
                           xlab=covariate,ylab=expression(hat(D)),main="",
                           Dlcol="black",Dlty=1,Dlwd=1,
                           CIlcol="black",CIlty=2,CIlwd=1,
                           tkcol="black",tkcex=0.75,tkpch="|") {
  smest = condDensity(fit,covariate=covariate,otherlevels=otherlevels,mask=mask,se=se,cl=cl,nx=nx) # calculate conditional density
#  if(is.null(xlab)) xlab=covariate
#  plot(smest$preds$x,smest$preds$Dhat*Dscale,ylim=range(0,na.omit(smest$preds[,-1]))*Dscale,
#       xlab=xlab,ylab=expression(hat(D)),type="n")
  plot(smest$preds$x,smest$preds$Dhat*Dscale,ylim=range(0,na.omit(smest$preds[,-1]))*Dscale,
       xlab=xlab,ylab=ylab,type="n",main=main)
  lines(smest$preds$x,smest$preds$Dhat*Dscale,col=Dlcol,lty=Dlty,lwd=Dlwd)
  if(cl) {
    lines(smest$preds$x,smest$preds$lcl*Dscale,col=CIlcol,lty=CIlty,lwd=CIlwd)
    lines(smest$preds$x,smest$preds$ucl*Dscale,col=CIlcol,lty=CIlty,lwd=CIlwd)
  }
  # add locations of camera traps
  points(smest$detcovs,rep(0,length(smest$detcovs)),pch=tkpch,cex=tkcex,col=tkcol)
  
  invisible(smest)
}



#' @title Plots 3D encounter rate on mesh
#' @param capthist is object of class "capthist".
#' @param mask If not NULL, is object of class "mask". 
#' @param occasions Scalar or vector of occasions to plot (encounters and effot
#' will be combined over these occasions).
#' @param add.points TRUE if you want points on top of the "pins" showing encounter
#' rate.
#' @param add.text TRUE if you want the ecounter rate printed on top of the "pins".
#' @param add.maskedge TRUE if you want the mask edge to be drawn (sometimes can
#' make mask look better, often not).
#' @param add.mask TRUE if you want the mask plotted on the base.
#' @param add.traps TRUE if you want the traps plotted as points on the base.
#' @param maskcov NULL, or name of the mask covariate you want to shade the mask with.
#' @param maskcol Name of the colour scheme to use in shading the mask.
#' @param ncols Number of colours to use in the mask colour scheme
#' @param pointcolr Name of color scheme to use in plotting points on "pins" 
#' (NULL gives you paula color scheme).
#' @param ptsize Size of points on "pins".
#' 
#' @export
plotER3d = function(capthist, mask=NULL, occasions=1, add.points=TRUE, 
                    add.text=FALSE, add.maskedge=FALSE, add.mask=TRUE, 
                    add.traps=TRUE, maskcov=NULL,maskcol=parula, ncols=40,
                    pointcolr=NULL,ptsize=5,...) {
  require(rgl)
  require(pals)
  
  nsessions = length(attr(capthist,"session"))
  if(nsessions>1) stop("Only one session at a time: capthist must be a matrix, not a list.")
  if(max(occasions)>dim(capthist)[2]) stop("occasion bigger than number of occasions in cpature history.")
  if(min(occasions)<1) stop("occasion must be greater than zero.")
  
  occasions = as.integer(occasions)
  noccasions = length(occasions)
  
  dets = traps(capthist)
  ntraps = dim(dets)[1]
  asp= c(1,diff(range(dets$y))/diff(range(dets$x)),1)
  if(is.null(usage(dets))) effort = rep(1,ntraps)
  else effort = apply(usage(dets)[,occasions,drop=FALSE],1,sum)
  ndets = matrix(rep(0,noccasions*ntraps),nrow=noccasions)
  for(i in 1:noccasions) ndets[i,] = apply(capthist[,occasions[i],,drop=FALSE],3,"sum")
  ndets = apply(ndets,2,sum)
  er = ndets/effort
  er[er==Inf] = 0
  zlim=range(0,er)
  if(!is.null(mask)){
    xlim = range(mask$x)
    ylim = range(mask$y)
    plot3d(dets$x,dets$y,er, size=10,type="h",lwd=1,xlim=xlim,ylim=ylim,zlim=zlim,
           aspect=asp,...)
  } else {
    if(add.mask | add.maskedge | !is.null(maskcov)) warning("No mask specified, so can't plot mask.")
    plot3d(dets$x,dets$y,er, size=10,type="h",lwd=1,zlim=zlim,
           aspect=asp,...)
  }
  rgl.bbox(xlen = 0, ylen = 0, zlen = 0, color = 'white')
  if(add.points) {
    if(is.null(pointcolr)) pointcolr = parula(max(ndets[er>0])+1)[ndets[er>0]+1]
    else pointcolr = rep(pointcolr[1],length((ndets[er>0])))
    points3d(dets$x[er>0],dets$y[er>0],er[er>0],col=pointcolr,size=ptsize,add=TRUE)
  }
  if(add.text) {
    textcolr = parula(max(ndets)+1)[ndets+1]
    text3d(dets$x,dets$y,1.05*er, texts=signif(er,2),col=textcolr,add=TRUE, cex=0.75)
  }
  
  if(!is.null(mask)){
    if(!add.maskedge & !add.mask) add.mask=TRUE # if pass mask, plot it!
    # add mask edge if asked to
    if(add.maskedge) {
      edge = plotMaskEdge(mask,plt=FALSE,add=TRUE)
      xedge = edge[c(1,3),]
      yedge = edge[c(2,4),]
      segments3d(xedge,yedge,z=-0.01*diff(zlim))
    }
    
    # add mask if asked to
    if(add.mask) {
      if(!is.null(maskcov)) {
        usecols = maskcol(ncols)
        covrange=range(covariates(mask)[maskcov])
        dcov = diff(covrange)
        covind = ceiling((covariates(mask)[maskcov]-min(covrange))/dcov * ncols)[,1]
        covind[covind==0] = 1
        covcols = usecols[covind]
        points3d(mask$x,mask$y,z=-0.01*diff(zlim),col=covcols)
      } else points3d(mask$x,mask$y,z=-0.01*diff(zlim),col="gray")
    }
  }
  # add traps
  if(add.traps) points3d(dets$x,dets$y,z=0,pch=1,col="red")
  
}



#' @title Plot top proportion of density
#' @param Dmask is object of class "Dsurface".
#' @param pc is top percentage of the abundance to plot. Must be betwen 0 and 100. 
#' For example, if p=25, the distribution of the highest density regions that accounts for 
#' 25% of the total abundance is plotted.
#' @param boundary is a SpatialPolygonsDataFrame defining the boundary of the sruvey region.
#' @param binar is a variable, which if true, causes the top pc to be assigned a single 
#' density value of 1, rather than its actual density.
plotNpc = function(Dmask,pc,boundary=NULL,binary=FALSE,...) {
  if(!inherits(Dmask,"Dsurface")) stop ("Dmask must be of class Dsurface.")
  if(pc<=0 | pc>=100) stop("pc must be less than 100 and greater than 0")
  pc = (100-pc)/100
  D = covariates(Dmask)$D.0
  ord = order(D)
  F = cumsum(D[ord])/sum(D)
  keep = which(F>=pc)
  covariates(Dmask)$Dpc = rep(NA,length(covariates(Dmask)$D.0))
  covariates(Dmask)$Dpc[ord[keep]] = covariates(Dmask)$D.0[ord[keep]]
  if(binary) covariates(Dmask)$Dpc[ord[keep]] = 1
  if(!is.null(boundary)) {
    plot(boundary,add=TRUE)
    plotcovariate(Dmask,covariate="Dpc",...,add=TRUE)
  } else {
    plotcovariate(Dmask,covariate="Dpc",...)
  }
  #  if(!is.null(boundary)) plot(boundary,add=TRUE)
}


