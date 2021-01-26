#' @title Calculate detection function or encounter rate parameters
#'
#' @description
#'  Calculates detection function or expected encounter rate at a range of
#'  distances, given detection function or encounter rate (hazard) function
#'  parameters. This is a very slightly hacked version of the \link{secr} function
#'  \link{secr::detectfnplot}, modified so as not to do the plotting and always 
#'  return the function values.
#'  
#' @param detectfn integer code or character string for shape of detection 
#' function 0 = halfnormal etc. – see \link{secr::detectfn}. 
#' @param pars list, vector or matrix of parameter values.
#' @param pars.vcv variance-covariance matrix of pars
#' @param details list of ancillary parameters - see \link{secr::detectfn} for 
#' details
#' @param rgr logical; if TRUE a scaled curve r.g(r) is plotted instead of g(r)
#' @param hazard logical; if TRUE the hazard of detection (the expected encounter
#' function) is plotted instead of the probability.
#' @param dists distances at which to calculate the function.
#' 
#' @export plam
#' 
plam = function(detectfn, pars, dists, pars.vcv=NULL, details=NULL, rgr=FALSE, hazard=FALSE) {
  gline <- function(pars) {
    dfn <- secr:::getdfn(detectfn)
#    if (sigmatick) {
#      sigma <- pars[2]
#      y <- dfn(sigma, pars, details$cutval)
#      if (hazard) 
#        y <- -log(1 - y)
#      dy <- par()$usr[4]/20
#      #segments(sigma, y - dy, sigma, y + dy)
#    }
    y <- dfn(dists, pars, details$cutval)
    if (hazard) 
      y <- -log(1 - y)
    if (rgr) {
      y <- dists * y
      ymax <- par()$usr[4]
      #lines(dists, y * 0.8 * ymax/max(y), lty = 2, ...)
    }
    #else lines(dists, y, ...)
    data.frame(x = dists, y = y)
  }
  if (is.list(pars)) {
    if (is.list(pars[[1]])) 
      pars <- matrix(unlist(pars), nrow = length(pars), byrow = T)
    else pars <- unlist(pars)
  }
  if (!is.matrix(pars)) 
    pars <- matrix(pars, nrow = 1)
  if (is.character(detectfn)) 
    detectfn <- secr:::detectionfunctionnumber(detectfn)
  needp <- c(2, 3, 2, 3, 2, 3, 3, 3, 3, 2, 3, 3, 5, 5, 2, 3, 2, 3, 3, 3)[detectfn + 1]
  if (ncol(pars) != needp) 
    stop("require ", needp, " parameters for ", detectionfunctionname(detectfn), " detection function")
#  if (is.null(ylim)) {
#    if (detectfn %in% c(10, 11, 12, 13)) {
#      ylim <- c(0, 1)
#    }
#    else {
#      ylim <- c(0, max(pars[, 1]))
#    }
#  }
  #  if (!add) {
  #    if (is.null(xlab)) 
  #      xlab <- "Distance  (m)"
  #    if (is.null(ylab)) 
  #      ylab <- "Detection"
  #    plot(type = "n", 0, 0, xlim = range(dists), ylim = ylim, 
  #         xlab = xlab, ylab = ylab, ...)
  #  }
  as.data.frame(apply(pars, 1, gline))
}




#' @title Draws a shaded confidence region for a line
#'
#' @description
#'  Adds a shaded polygon covering the region between lower and upper confidence
#'  bounds, to an existing plot.
#'  
#' @param xyci a data frame or matrix with first element being x-coordinates, 
#' second being the corresponding lower (or upper) confidence bound values, 
#' third being the corresponding upper (or lower) confidence bound values. 
#' @param col an argument of \code{polygon} specifying shading colour.
#' @param border an argument of \code{polygon} specifying border colour
#' @param ... other arguments to \code{polygon}.
#' 
#' @export shadedCI
#' 
shadedCI = function(xyci,col=adjustcolor("gray", alpha.f=0.4),border=NA,...) {
  if(!inherits(xyci,"data.frame") & !inherits(xyci,"matrix")) 
    stop("xyci must be a matrix or data frame.")
  if(ncol(xyci)!=3) stop("xyci must have 3 columns.")
  ci.x <- c(xyci[,1], rev(xyci[,1]))
  ci.y <- c(xyci[,2], rev(xyci[,3]))
  polygon(x=ci.x, y=ci.y,col=col,border=border,...)
}



#' @title Replaces substring(s) in a string
#'
#' @description
#'  Replaces all instances of specified substrings in a string (character variable) 
#'  with specified other substrings.
#'  
#' @param string the string (character variable) in which replacing is to occur.
#' @param from a (vector of) character variables that are to be replaced.
#' @param to a (vector of) character variables of the same length as \code{from}
#' containing the strings that are to replace \code{from}.
#' 
#' @export streplace
#' 
streplace = function(string,from,to){
  n2rep = length(from)
  if(length(to)!=length(from)) stop("Lengths of `from` and `to` must be the same.")
  nch = nchar(string)
  outch = vector(mode="character",length=0)
  for(n in 1:n2rep) {
    nfrom = nchar(from[n])
    if(nfrom<nch) {
      i = 1
      while(i<=nch) {
        if(substr(string,i,i+nfrom-1)==from[n]) {
          outch = paste(outch,to[n],sep="")
          i = i+nfrom
        }else {
          outch = paste(outch,substr(string,i,i),sep="")
          i = i+1
        }
      }
      string = outch
      nch = nchar(string)
      outch = vector(mode="character",length=0)
    }
  }
  return(string)
}

#' @title Summarises density estimates from secr objects
#'
#' @description
#'  Makes a data frame summarising density estimates, ordered by AIC or AICc.
#'  
#' @details 
#' Lists the following in each column of the data frame, for each model, with models arranged in increasing order of AIC or AICc
#' (whichever was specified in the call). Density estimate, SE, lower- and upper 95\% confidence intervals are multiplied by 
#' argument \code{Dscale}
#' \itemize{
#'  \item{Dhat:}{ Density estimate.}
#'  \item{se:}{ Standard error of density estimate.}
#'  \item{lcl:}{ Lower 95\% confidence bound for density.}
#'  \item{ucl:}{ Upper 95\% confidence bound for density.}
#'  \item{cv:}{ Coefficient of variation of density estimate.}
#'  \item{npar:}{ Number of parameters in the model.}
#'  \item{logLik:}{ Log-Likelihood}
#'  \item{AIC:}{ AIC or AICc (whichever was specified in the call).}
#'  \item{dAIC:}{ Delta AIC or Delta AICc (whichever was specified in the call).}
#'  \item{AICwt:}{ AIC weight or AICc weight (whichever was specified in the call).}
#'  }
#'  Requires function \code{AIC} from package \code{secr}.
#'  
#' @param allargs Object of class '\code{secrlist}' or '\code{secr}' from package \code{secr}.
#' @param ... Other object of class '\code{secrlist}' or '\code{secr}'.
#' @param criterion Character variable being 'AIC' or 'AICc'. 
#' @param Dscale Numeric scaling factor for density (which is multiplied by this factor).
#' 
#' @export makeDtable
makeDtable = function(object, ..., criterion=c("AICc","AIC"),Dscale=1) {
  if(class(object)!="secr" & class(object)!="secrlist") stop("object must be of class `secr` or `secrlist`")
  allargs <- list(...)
  modelnames <- (c(as.character(match.call(expand.dots = FALSE)$object), 
                   as.character(match.call(expand.dots = FALSE)$...)))
  allargs <- secrlist(object, allargs)
  names(allargs) <- modelnames
  nests = length(allargs)
  aicvals = nDpars = rep(NA,nests)
  NAs = rep(NA,nests)
  Dhat = data.frame(Dhat=NAs, se=NAs, lcl=NAs, ucl=NAs, cv=NAs, npar=NAs, logLik=NAs, AIC=NAs, dAIC=NAs, AICwt=NAs, model=NAs,
                    name=NAs)
  for(i in 1:nests) {
    if(class(allargs)=="secr") secrobj = allargs else secrobj = allargs[[i]]
    aicobj = AIC(secrobj,criterion=criterion[1])
    aicvals[i] = as.numeric(aicobj[which(names(aicobj)==criterion[1])])
    nDpars[i] = length(which(substr(row.names(coefficients(secrobj)),start=1,stop =1)=="D"))
  }
  mod = order(aicvals)
#  if(class(allargs)=="secr") secrobj = allargs else secrobj = allargs
#  aics <- AIC(secrobj)
  aics <- AIC(allargs)
  for(i in 1:nests) {
    msum = summary(allargs[[i]]$mask)
    A = msum$cellarea*msum$nmaskpoints
    if(nDpars[i]==1) {
      Dhat[i,1:5] = derived(allargs[[i]])[2,c(1:4,7)]
    }else {
      Dhat[i,1:4] = region.N(allargs[[i]])[1,1:4]/A
      Dhat[i,5] = Dhat[i,2]/Dhat[i,1]
    }
  }
  Dhat[,1:4] = Dhat[,1:4]*Dscale
  Dhat$cv = Dhat$cv*100 # make a percentage
#  allaics = AIC(allargs)
  Dhat$npar = aics$npar
  Dhat$logLik = aics$logLik
  Dhat$AIC = aicvals[mod]
  Dhat$name = modelnames[mod]
  if(nests>1) {
    Dhat$dAIC = Dhat$AIC-min(Dhat$AIC)
    Dhat$AICwt = aics[,8]
  }
  Dhat$model = aics$model
  
  return(Dhat)
}




#' @title Transforms secr object parameters from 'beta' to natural scale
#'
#' @description
#'  Transforms 'beta' parameters, SEs and confidence bounds onto natural scale.
#'  
#' @param fit An \code{secr} fitted object.
#' 
#' @details This function adapts the unexported \code{secr} functions \code{secr:::untransform} and
#'  \code{secr:::se.untransform} to operate directly on an object of class \code{secr}.
#' 
#' @export beta2natural
beta2natural = function(fit) {
  pnames = fit$betanames
  beta = coefficients(fit)
  link = fit$link
  npar = length(pnames)
  est = beta*0
  names(est)[1] = "par"
  names(est)[2] = "SE.par"
  for(i in 1:npar) {
    est[i,1] = switch(link[[i]], identity=beta[i,1], log=exp(beta[i,1]), neglog=-exp(beta[i,1]), 
                      logit=invlogit(beta[i,1]), odds=invodds(beta[i,1]), sin=invsine(beta[i,1]))
    est[i,2] = switch(link[[i]], identity=beta[i,1], log=exp(beta[i,1])*sqrt(exp(beta[i,2]^2)-1), 
                      neglog=exp(beta[i,1])*sqrt(exp(beta[i,2]^2)-1), logit=invlogit(beta[i,1]) * 
                        (1-invlogit(beta[i,1]))*beta[i,2], sin = NA)
    est[i,3] = switch(link[[i]], identity=beta[i,3], log=exp(beta[i,3]), neglog=-exp(beta[i,3]), 
                      logit=invlogit(beta[i,3]), odds=invodds(beta[i,3]), sin=invsine(beta[i,3]))
    est[i,4] = switch(link[[i]], identity=beta[i,4], log=exp(beta[i,4]), neglog=-exp(beta[i,4]), 
                      logit=invlogit(beta[i,4]), odds=invodds(beta[i,4]), sin=invsine(beta[i,4]))
  }
  return(est)
}

#' @title Calculate distances between points
#'
#' @description
#'  Returns all the distances between two sets of points, in a matix with as many rows as there are points 
#'  in \code{from} and as many columns as there are points in \code{t}..
#'  
#' @param from A data frame or matrix with two columns, being the x- and y-coordinates of a set of points
#' @param to A data frame or matrix with two columns, being the x- and y-coordinates of a set of points
#' 
#' @export distances
distances = function(from,to) {
  x = rbind(as.matrix(from),as.matrix(to))
  ds = as.matrix(dist(x))[1:dim(from)[1],(dim(from)[1]+1):dim(x)[1]]
}


#' @title Draws histogram.
#'
#' @description
#'  Utility function to draw histograms with more options than \code{hist} allows.
#'  
#' @param height Height of histogram bars.
#' @param breaks Locations of boundaries of histogram bins (must be 1 longer than \code{height}).
#' @param lineonly If TRUE, uses \code{\link{lines}} to draw lines on current plot; else uses 
#' \code{\link{plot}} to draw lines on new plot.
#' @param outline If TRUE, draws only the outline (profile) of the histogram; else draws each 
#' complete bar.
#' @param fill If TRUE, uses polygon() to fill barsl in this case valid arguments to polygon() 
#' are passed via argument(s) "...". If fill==FALSE, valid arguments to plot() or lines() are 
#' passed via argument(s) "..."
#' @param ylim Range of y-axis.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param ... See aargument \code{fill}.
#' 
#' @export histline
histline = function(height,breaks,lineonly=FALSE,outline=FALSE,fill=FALSE,ylim=range(height),
                  xlab="x",ylab="y",...)
{
  n=length(height)
  if(length(breaks)!=(n+1)) stop("breaks must be 1 longer than height")
  if(outline) {
    y=c(0,rep(height,times=rep(2,n)),0)
    x=rep(breaks,times=rep(2,(n+1)))
  }   else {
    y=rep(0,4*n)
    x=rep(0,4*n+2)
    for(i in 1:n) {
      y[((i-1)*4+1):(i*4)]=c(0,rep(height[i],2),0)
      x[((i-1)*4+1):(i*4)]=c(rep(breaks[i],2),rep(breaks[i+1],2))
    }
    x=x[1:(4*n)]
  }
  if(lineonly) {
    if(!fill) lines(x,y,...)
    else polygon(x,y,...)
  } else {
    if(!fill) plot(x,y,type="l",ylim=ylim,xlab=xlab,ylab=ylab,...)
    else {
      plot(x,y,type="n",ylim=ylim,xlab=xlab,ylab=ylab)
      polygon(x,y,...)
    }
  }
}


