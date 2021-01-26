#' @title Calculates encounter or detection function values and CIs
#'
#' @description
#'  Calculates an encounter function or a detection function at specified distances, 
#'  given an \code{secr} object, a data frame containing covariates used by the 
#'  object, and the distances.
#'  Returns a data frame of x (distance), y (detection function) and (optionally) 
#'  lower and upper confidence interval values, calculated either using the 
#'  method of the package \code{secr}, or by parameteric bootstrap..
#'  
#' @param x an object of class \code{secr}. 
#' @param newdata a data frame with a column containing appropriate covariate
#' values for each of the covariates in \code{x} (both density model and 
#' detection or encounter model covariates, even though denstiy covariates are 
#' not used for detection function calculation).
#' @param dists distances at which to calculate the detection function.
#' @param rgr TRUE if you want distance multiplied by detection function value
#' to be returned (as per point transects, for example).
#' @param ci If TRUE, confidence intervals are calculated and returned.
#' @param ci.method If 'secr' then the method of package \code{secr} (i.e. a kind
#' of Delta Method) is used to calculate confidence intervals, else if 'bootstrap' 
#' nonparametric bootstrap, resampling from a multivariate normal with mean 
#' equal to the estimated model parameters and variance-covariance equal to 
#' their estimated variance-covariance matrix, is used.
#' @param alpha Significance level for confidence intervals
#' @param what If 'p' the detection function is returned; if 'er' the encounter 
#' function is returned.
#' @param B number of bootstrap samples to take
#' 
#' @export erdetfun
#' 
erdetfun = function(x, newdata=NULL, dists=0:200, rgr=FALSE, ci=FALSE, ci.method="secr", alpha=0.05, what="p",B=999) {
  if(!is.element(ci.method,c("secr","bootstrap"))) stop("method must be `secr` or `bootstrap`.")
  if(!is.element(what,c("p","er"))) stop("method must be `p` (for detection probability) or `er` (for encouner rate).")
  # get detection function
  pdat = noplot.secr(x=x, newdata=newdata, xval=dists, rgr=rgr, limits=ci,alpha=alpha)
  npcols = dim(pdat)[2]
  if(ci.method=="bootstrap" & ci) {
    require(MASS) # needed for function mvrnorm
    require(tcltk2)
    bpars = mvrnorm(B, x$fit$par,x$beta.vcv) # bootstrap parameters
    ndists = length(dists)
    bpest = matrix(rep(NA,B*ndists),ncol=B) # bootstrapped det funs
    # Set up progress bar, before looping:
    pb <- tkProgressBar(title=paste("Bootstrap Progress (B=",B,")",sep=""), min=0, max=B, width=400)
    bx = x
    for(i in 1:B) {
      # Progress, inside loop
      setTkProgressBar(pb, i, label=paste( round(i/B*100, 0),"% done"))
      # put bootstrap pars in secr object bx
      bx$fit$par = bpars[i,]
      # estimate detection prob without CI, with bootrap pars
      bpest[,i] = noplot.secr(x=bx, newdata=newdata, xval=dists, rgr=rgr, limits=FALSE)$y
    }
    bci = t(apply(bpest,1,quantile,probs=c(alpha/2,1-alpha/2)))
    pdat[,3:4] = bci
    # Close progress bar after looping
    close(pb)
  }
  if(what=="er") pdat[,2:npcols] = -log(1-pdat[,2:npcols])
  return(pdat)
}


#' @title Calculates detection function values
#'
#' @description
#'  Calculates a detection function at specified distances, given an \code{secr}
#'  object, a data frame containing covariates used by the object, and the 
#'  distances
#'  This is a hacked version of the \link{secr} function \code{plot.secr} - 
#'  hence the name '\code{noplot.secr}', The hack is so as not to plot, and to 
#'  return a data frame of x (distance), y (detection function) and (optionally) 
#'  lower and upper confidence interval values.
#'  
#' @param x an object of class \code{secr}. 
#' @param newdata a data frame with a column containing appropriate covariate
#' values for each of the covariates in \code{x} (both density model and 
#' detection or encounter model covariates, even though denstiy covariates are 
#' not used for detection function calculation).
#' @param rgr TRUE if you want distance multiplied by detection function value
#' to be returned (as per point transects, for example).
#' @param alpha Significance level for confidence intervals
#' @param xval distances at which to calculate the detection function.
#' 
#' @export noplot.secr
#' 
noplot.secr = function (x, newdata=NULL, rgr=FALSE, limits=FALSE, alpha=0.05, xval=0:200) 
{
  gline <- function(predicted, rowi = 1, eps = 1e-10) {
    if (!is.data.frame(predicted)) {
      out <- list()
      for (i in 1:length(predicted)) out[[i]] <- gline(predicted[[i]],i)
      names(out) <- names(predicted)
      out
    }
    else {
      pars <- predicted[secr:::parnames(x$detectfn), "estimate"]
      pars[is.na(pars)] <- unlist(x$fixed)
      dfn <- secr:::getdfn(x$detectfn)
#      if (sigmatick) {
#        sigma <- pars[2]
#        y <- dfn(sigma, pars, x$details$cutval)
#        dy <- par()$usr[4]/20
 #     }
      y <- dfn(xval, pars, x$details$cutval)
      if (rgr) {
        y <- xval * y
        ymax <- par()$usr[4]
      }
      if (limits & !rgr) {
        grad <- matrix(nrow = length(xval), ncol = length(x$fit$par))
        if (is.null(newdata)) 
          newdata <- secr.make.newdata(x)
        parnamvec <- secr:::parnames(x$detectfn)
        if (!parnamvec[1] %in% c("g0", "lambda0", "beta0")) 
          stop("first detection parameter not g0 or lambda0")
        lkdfn <- function(beta, r) {
          real <- numeric(length(parnamvec))
          names(real) <- parnamvec
          for (rn in parnamvec) {
            par.rn <- x$parindx[[rn]]
            mat <- secr:::general.model.matrix(x$model[[rn]], 
                                               data = newdata[rowi, , drop = FALSE], gamsmth = x$smoothsetup[[rn]], 
                                               contrasts = x$details$contrasts)
            lp <- mat %*% matrix(beta[par.rn], ncol = 1)
            real[rn] <- secr:::untransform(lp, x$link[[rn]])
          }
          gr <- dfn(r, real, x$details$cutval)
          if (parnamvec[1] == "lambda0") 
            log(-log(1 - gr))
          else logit(gr)
        }
        for (i in 1:length(xval)) grad[i, ] <- secr:::gradient(pars = x$fit$par, 
                                                               fun = lkdfn, r = xval[i])
        vc <- vcov(x)
        gfn <- function(gg) {
          gg <- matrix(gg, nrow = 1)
          gg %*% vc %*% t(gg)
        }
        se <- apply(grad, 1, gfn)^0.5
        if (parnamvec[1] == "lambda0") {
          lcl <- ifelse((y > eps) & (y < (1 - eps)), 
                        1 - exp(-exp(log(-log(1 - y)) - z * se)), 
                        NA)
          ucl <- ifelse((y > eps) & (y < (1 - eps)), 
                        1 - exp(-exp(log(-log(1 - y)) + z * se)), 
                        NA)
        }
        else {
          lcl <- ifelse((y > eps) & (y < (1 - eps)), 
                        invlogit(logit(y) - z * se), NA)
          ucl <- ifelse((y > eps) & (y < (1 - eps)), 
                        invlogit(logit(y) + z * se), NA)
        }
      }
      if (limits & !rgr) 
        data.frame(x = xval, y = y, lcl = lcl, ucl = ucl)
      else data.frame(x = xval, y = y)
    }
  }
  z <- abs(qnorm(1 - alpha/2))
  temp <- predict(x, newdata)
  return(gline(temp))
}
