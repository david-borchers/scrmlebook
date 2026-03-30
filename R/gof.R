#' @title Calculates encounter or detection function.
#'
#' @description
#' Calculates encounter function or detection function, scaled by distance or not.
#'
#' @param dists a matrix of distances.
#' @param detectfn \code{secr} detection function type (number or character).
#' @param pars list like that returned by the \code{secr} function \code{detectpar}.
#' @param hazard If TRUE, returns encounter function, else detection function.
#' @param effort Trapping effort (multiplies hazard). Must be a vector of length
#' \code{nrow(dist)}.
#' @param rgr If TRUE multiplies encounter or detection function by distance.
#' 
#' @description Calculates encounter function or detection function at given \code{dists}.
#' 
#' @examples
#' detectors <- make.grid (nx = 6, ny = 8, detector = "multi")
#' detections <- sim.capthist (detectors, popn = list(D = 10, buffer = 100), 
#'                             detectpar = list(g0 = 0.2, sigma = 25))
#' ## fit & print null (constant parameter) model
#' model <- secr.fit (detections)
#' dists = matrix(c(0,10,20,20,10,0),nrow=2,byrow=2)
#' pars = detectpar(model)
#' detectfn = model$detectfn
#' predetectfn(dists,detectfn,pars)
#' predetectfn(dists,detectfn,pars,effort=c(1,2)))
#' predetectfn(dists,detectfn,pars,effort=c(1,2)),hazard=FALSE)
#' predetectfn(dists,detectfn,pars,hazard=FALSE)
#' 
#' @export
predetectfn = function(dists, detectfn, pars, effort=rep(1,nrow(dists)), hazard=TRUE, rgr=FALSE) {
  
  require(secr)
  calchaz <- function(dists) {
    dfn <- secr:::getdfn(detectfn)
    y <- dfn(dists, pars)
    y <- -log(1 - y)
    return(y)
  }
  
  if (is.character(detectfn)) detectfn <- secr:::detectionfunctionnumber(detectfn)
  if (is.list(pars)) {
    if (is.list(pars[[1]])) 
      pars <- matrix(unlist(pars), nrow = length(pars), byrow = T)
    else pars <- unlist(pars)
  }
  if (!is.matrix(pars)) pars <- matrix(pars, nrow = 1)
  
  haz = apply(dists, 1, calchaz)
  y = t(haz%*%diag(as.vector(effort)))
  if (!hazard) 
    y <- 1-exp(-y)
  if (rgr) {
    y <- dists * y
  }
  
  return(y)
}

#' @title Calculates Freeman-Tukey GoF statistic.
#'
#' @description
#' Calculates the Freeman-Tukey GoF statistic. If \code{object} is of class \code{secr}, then 
#' \code{En.capt} is ignored and the expected number of captures is calculated from \code{object} 
#' and the statistic is calculated from  \code{object$capthist}. If If \code{object} is of class 
#' \code{capthist} then \code{En.capt} must contain the expected number of captures and the 
#' statistic is calculated from this and \code{capthist}.
#'
#' @param object An object of class \code{secr} fitted model (containing an element \code{$mask}), 
#' or an object of class \code{capthist}.
#' @param En.capt The expected number of captures in total.
#' 
#' @examples
#' require(secr)
#' Freeman_Tukey(secrdemo.0)
#' 
#' @export
Freeman_Tukey = function(object,En.capt=NULL) { 
  if(class(object)=="secr") {
    ch = object$capthist
    En.capt = calcEn(object)
  }else if(class(object)=="capthist") {
    ch = object
    if(is.null(En.capt)) stop("Need to pass a non-null En.capt.")
  }else {
    stop("Invalid class of object. Class must be 'secr' or 'capthist'.")
  }
  n = apply(ch,1,sum) # observed number of captures of each detected animal
  FTstat = sum((n-En.capt)^2) # test statistic
  return(FTstat)
}


#' @title Calculates expected number of captures.
#'
#' @description
#' Calculates the expected number of captures from an \code{secr} object.
#'
#' @param fit An \code{secr} fitted model. (This must have an element \code{$mask}.)
#' 
#' @examples
#' require(secr)
#' calcEn(secrdemo.0)
#' 
#' @export
calcEn = function(fit) {
  nocc = dim(fit$capthist)[2]
  dets = traps(fit$capthist)
  mesh = predictDsurface(fit,fit$mask)
  D = covariates(mesh)$D.0
  dists = edist(dets,mesh)
  En.capt = 0
  for(occ in 1:nocc) {
    effort=rep(1,nrow(dists))
    if(!is.null(usage(dets))) effort = as.numeric(usage(dets)[,occ])
    lambdas = predetectfn(dists, detectfn=fit$detectfn, detectpar(fit), effort, hazard=TRUE, rgr=FALSE)
    Ens = apply(lambdas,2,sum) # E(number detectionsd at any trap | AC)
    p.s = 1 - exp(-Ens) # p(detected | AC)
    En = sum(Ens*D)/sum(D)
    D.capt = p.s*D/sum(p.s*D)
    En.capt = En.capt + sum(Ens*D.capt)/sum(D.capt) # expected number of captures of each detected animal
  }
  return(En.capt)
}

#' @title Calculates Deviance/df statistic.
#'
#' @description
#' Calculates the scaled deviance divided by the model degrees of freedom.
#'
#' @param fit An \code{secr} fitted model. (This must have an element \code{$mask}.)
#' @param ... Unused argument, here to make use in \code{gof.sim} easy (as some functions
#' passed to \code{gof.sim} have other arguments, and this dummy argument allows the same
#' syntax to be used for them and this.)
#' 
#' @examples
#' require(secr)
#' Devdf(secrdemo.0)
#' 
#' @export
Devdf = function(fit, ...) return(deviance(fit) / df.residual(fit)) 


#' @title Simulate sample of a GoF test statistic.
#'
#' @description
#' Repeatedly simulates a population and capture history using the parameters of a 
#' fitted \code{secr} model as truth, fits a model of the same sort to the data,
#' and calculates the test statistic (specified via\code{statfn}) for each new
#' fitted model.
#' 
#' Only works with full likelhood models (those that include a model for density), and 
#' does not work for all of them. In particular, it works for models without
#' individual effects, except that it does work for a model with sex as a 
#' partially-observed individual effect in which both sexes have a common density.
#' 
#' Has only been tested with single-occasion surveys.
#'
#' @param fit An \code{secr} fitted model. (This must have an element \code{$mask}.)
#' @param statfn A function that calculates a  GoF test statistic from a fitted model or 
#' \code{capthist} object.
#' @param nsim The number of simulations to do.
#' @param newfit If TRUE, refits model to every simulated dataset and recalculates expected numbers.
#' @param seed Random number seed (setting it gives repeatability of results).
#' @param ... Other arguments to pass to statfn.
#' 
#' @examples
#' require(secr)
#' sim.gof(secrdemo.0)
#' 
#' @export
sim.gof = function(fit, statfn=Devdf, nsim=99, newfit=FALSE, seed=NULL, ...) {
  
  if(!is.null(seed)) set.seed(seed)
  
  require(tcltk2) # for progress bar function
  # Set up progress bar, before looping:
  pb <- tkProgressBar(title=paste("Simulation Progress (Nsim=",nsim,")",sep=""), min=0, max=nsim, width=400)
  
  dets = traps(fit$capthist)
  mesh = predictDsurface(fit,fit$mask)
  D = covariates(mesh)$D.0
  simres = vector(length=nsim,mode="list")
  if(is.null(fit$hcov)) {
    for(i in 1:nsim) {
      simpop = sim.popn(D=D, core=mesh, buffer=0, Ndist="poisson", model2D="IHP")
      simch = sim.capthist(dets, popn=simpop, 
                           detectfn=fit$detectfn, detectpar=detectpar(fit),
                           noccasions=dim(fit$capthist)[2])
      if(newfit) simres[[i]] = secr.fit(simch,model=fit$model,mask=mesh,detectfn=fit$detectfn,trace=0)
      else simres[[i]] = statfn(simch, ...)
      
      # Progress, inside loop
      setTkProgressBar(pb, i, label=paste( round(i/nsim*100, 0),"% done"))
    }
  } else if(fit$hcov=="sex") { # SPECIFIC model with sex structure and same density surface for both sexes
    pM = plogis(fit$fit$par[which(fit$betanames=="pmix.h2M")]) # proportion males
    for(i in 1:nsim) {
      # simulate males and females separately
      female.popn = sim.popn(D=D*(1-pM), core=mesh, buffer=0, Ndist="poisson", model2D="IHP") # females
      male.popn = sim.popn(D=D*pM, core=mesh, buffer=0, Ndist="poisson", model2D="IHP") # females
      covariates(female.popn)$sex="F"
      covariates(male.popn)$sex="M"
      
      # set detection parameters for simulated capture histories
      dpars = detectpar(fit,byclass=TRUE)
      female = 1; male = 2 # arbitrary, but makes no difference
      sigma.male = dpars[[male]]$sigma; sigma.female = dpars[[female]]$sigma
      lambda0.male = dpars[[male]]$lambda0; lambda0.female = dpars[[female]]$lambda0
      g0.male = 1-exp(-lambda0.male); g0.female = 1-exp(-lambda0.female)
      
      # Make count dataset
      ch.male = sim.capthist(dets, male.popn, detectfn=fit$detectfn,
                             detectpar=list(lambda0=lambda0.male, sigma=sigma.male),
                             noccasions=dim(fit$capthist)[2])
      ch.female = sim.capthist(dets, female.popn, detectfn=fit$detectfn,
                               detectpar=list(lambda0=lambda0.female, sigma=sigma.female),
                               noccasions=dim(fit$capthist)[2])
      simch = rbind(ch.male,ch.female)
      if(newfit) simres[[i]] = secr.fit(simch,model=fit$model,mask=mesh,detectfn=fit$detectfn,trace=0)
      else simres[[i]] = statfn(simch,...)
      
      # Progress, inside loop
      setTkProgressBar(pb, i, label=paste( round(i/nsim*100, 0),"% done"))
    }
  } else stop("hcov not null but not `sex`; can't deal with this case at present")
  
  if(newfit) Tdbn = unlist(lapply(simres,statfn,...))
  else Tdbn = unlist(simres)
  
  # Close progress bar after looping
  close(pb)
  
  return(Tdbn)
}


#' @title Simulate sample of a GoF test statistic.
#'
#' @description
#' Repeatedly simulates a population and capture history using the parameters of a 
#' fitted \code{secr} model as truth, fits a model of the same sort to the data,
#' and calculates the test statistic (specified via\code{statfn}) for each new
#' fitted model.
#' 
#' Only works with full likelhood models (those that include a model for density), and 
#' does not work for all of them. In particular, it works for models without
#' individual effects, except that it does work for a model with sex as a 
#' partially-observed individual effect in which both sexes have a common density.
#' 
#' Has only been tested with single-occasion surveys.
#'
#' @param object An \code{secr} fitted model. (This must have an element \code{$mask}.)
#' @param statfn A function that calculates a  GoF test statistic from a fitted model or 
#' \code{capthist} object.
#' @param nsim The number of simulations to do.
#' @param newfit If TRUE, refits model to every simulated dataset and recalculates expected numbers.
#' @param plot If TRUE, produces a histogram of simulated statistics with observed statistic marked.
#' @param nbins Number of histogram bins in plot
#' @param seed Random number seed (setting it gives repeatability of results).
#' 
#' @examples
#' require(secr)
#' gofdbn1 = gof_simtest(secrdemo.0,nsim=99,seed=1)
#' gofdbn1$pvalue
#' 
#' @export
gof_simtest = function(object, statfn=Devdf, nsim=99, newfit=FALSE, plot=TRUE, nbins=20, seed=NULL) {

  fn = as.character(substitute(statfn))
  En.capt = NULL
  if(class(object)=="secr") {
    if(!newfit) {
      if(fn=="Devdf") {
        newfit=TRUE
        warning("Argument 'newfit' changed to TRUE because have to refit for deviance statistic.")
      }else En.capt = calcEn(object)
    }
    Tobs = statfn(object) # observed value
    Tdbn = sim.gof(object, statfn=statfn, nsim=nsim, newfit=newfit, seed=seed, En.capt=En.capt) # simulated distribution
  }else if(class(object)=="capthist") {
    En.capt = calcEn(object)
    Tobs = statfn(object,En.capt) # observed value
    Tdbn = sim.gof(object, statfn=statfn, nsim=nsim, newfit=newfit, seed=seed, En.capt=En.capt) # simulated distribution
  }else {
    stop("Invalid class of object. Class must be 'secr' or 'capthist'.")
  }
  
  Ts = c(Tdbn,Tobs)
  pval = 1 - rank(Ts,ties.method="average")[length(Ts)]/length(Ts)
  if(plot) {
    if(fn=="Devdf") xlab = expression(D/df)
    if(fn=="Freeman_Tukey") xlab = "Freeman Tukey statistic"
    ylim = range(Ts)
    breaks = seq(ylim[1],ylim[2],length=nbins+1)
    hist(Tdbn,breaks=breaks,xlab=xlab,main="")
    points(Tobs,0,pch="*",cex=3)
  }
  print(paste("p-value =",signif(pval,3)))
  invisible(list(observed=Tobs, simulated=Tdbn, pvalue=pval))
}
