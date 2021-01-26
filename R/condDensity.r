#' @title Conditional value of fiited density smooth
#'
#' @description  Calculates the estimated density over the range of values of one spatial covariate, at
#' fixed values of all other spatial covariates.
#' 
#' @return A list comprising:
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
#' 
#' @seealso \code{\link{condDensityplot}}
#' 
#' @examples 
#' \dontrun{
#' library(scrmlebook)
#' data("leopard-sim") # simulated Kruger Park data
#' # Fit overly complex density model - just for illustration:
#' fit = secr.fit(kruger.capt,model=list(D~x+s(y,k=3)+s(habitat.cov,k=3)),kruger.mask)
#' smest = condDensity(fit,covariate="y") # calculate conditional density
#' Dscale = 100^2 # scale up to number per sq km
#' # Plot conditional density in Northing dimension
#' plot(smest$preds$x,smest$preds$Dhat*Dscale,ylim=range(0,na.omit(smest$preds[,-1]))*Dscale,
#'      xlab="Northing",ylab=expression(hat(D)),type="n")
#' lines(smest$preds$x,smest$preds$Dhat*Dscale)
#' lines(smest$preds$x,smest$preds$lcl*Dscale,lty=2)
#' lines(smest$preds$x,smest$preds$ucl*Dscale,lty=2)
#' # add locations of camera traps
#' points(smest$detcovs,rep(0,length(smest$detcovs)),pch="|",cex=0.75,col="black")
#' # Plot conditional density in dist.to.water dimension
#' smest = condDensity(fit,covariate="habitat.cov") # calculate conditional density
#' plot(smest$preds$x,smest$preds$Dhat*Dscale,ylim=range(0,na.omit(smest$preds[,-1]))*Dscale,
#'      xlab="habitat.cov",ylab=expression(hat(D)),type="n")
#' lines(smest$preds$x,smest$preds$Dhat*Dscale)
#' lines(smest$preds$x,smest$preds$lcl*Dscale,lty=2)
#' lines(smest$preds$x,smest$preds$ucl*Dscale,lty=2)
#' # add locations of camera traps
#' points(smest$detcovs,rep(0,length(smest$detcovs)),pch="|",cex=0.75,col="black")
#' }
#' @export
#' 
condDensity = function(fit,covariate,otherlevels=NULL,mask=NULL,se=TRUE,cl=TRUE,nx=200) {
  # need some class and value error traps
  if(!inherits(fit,"secr")) stop("Object 'fit' must be of class 'secr'.")
  if(is.null(mask)) mask = fit$mask
  # add x and y to mask
  covariates(mask)$x = mask$x
  covariates(mask)$y = mask$y
  # get covariates at detectors:
  dets = traps(fit$capthist)
  dets = addCovariates(dets,mask,columns=covariate)
  detcovs = covariates(dets)[,covariate]
  # make mask with all necessary covariates:
  covs = covariates(mask)
  covars = modcovs = all.vars(fit$model[["D"]])
#  covs$x = mask$x # always need x
#  covs$y = mask$y # always need y
#  vars = c(covars,"x","y")
  nvar = length(covars)
  if(!is.element(covariate,covars)) stop("Covariate ",covariate," is not in model",sep="")
  condcov = covs[,covariate]
  # mat is twice as long as specified, in order to have contrast in x and y when use read.mask() below, 
  # if want prediction at average x and y
  mat = matrix(data=c(rep(NA,nvar*nx*2)),ncol=nvar) 
  df = as.data.frame(mat)
  covnum = which(covars==covariate)
  names(df) = covars
  df$x = rep(seq(min(covs$x),max(covs$x),length=nx),2)
  df$y = rep(seq(min(covs$y),max(covs$y),length=nx),2)
  df[,covnum] = rep(seq(min(condcov),max(condcov),length=nx),2)
  for(i in 1:nvar) { # don't want to include last two variables, x and y, in this
    if(i!=covnum & covars[i]!="x" & covars[i]!="y") {
      if(is.element(covars[i],names(otherlevels))) df[,covars[i]] = otherlevels[[which(names(otherlevels)==covars[i])]]
      else df[,covars[i]] = mean(covs[,which(names(covs)==covars[i])]) # insert values for conditioned covariates
    }
  }
  if(is.element("x",covars) & is.element("y",covars)) { # need additional mask points with replicated
    if(covariate == "x") df$y[1:nx] = rep(mean(df$y))
    else if(covariate == "y") df$x[1:nx] = rep(mean(df$x))
    else {
      df$x[1:nx] = rep(mean(df$x))
      df$y[1:nx] = rep(mean(df$y))
    }
  }
  predmask = read.mask(data=df,columns=covars[covars!="x" & covars!="y"])
  predmask = subset(predmask,subset=1:nx)
  # predict using constructed mask:
  condest = predictDsurface(fit,mask=predmask,se.D=se,cl.D=cl)
  xcolumn = which(names(df)==covariate)
  condvals = df[1,-xcolumn]
  ses = lcl = ucl = rep(NA,nx)
  if(se) ses = covariates(condest)$SE.0
  if(cl) {
    lcl = covariates(condest)$lcl.0
    ucl = covariates(condest)$ucl.0
  }
  preds = data.frame(x=df[1:nx,xcolumn],Dhat=covariates(condest)$D.0, se=ses, lcl=lcl, ucl=ucl)
  
  return(list(preds=preds, detcovs=detcovs, condvals=condvals, xvar=covariate)) 
}
