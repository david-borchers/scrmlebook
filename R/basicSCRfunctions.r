#' @title SCR negative log-likelihood
#'
#' @description  Calculates the log-likelihood for an SCR model. Currently only multi-catch traps and 
#' count detectors with a single session implemented, and no detection function or density model covariates implemented.
#' 
#' @return Value of negative log-likelihood
#' 
#' @param pars Vector of paramters, named "g0", "sigma" and "D".
#' @param capthist Capture history object, or library \code{secr} class 'capthist'.
#' @param mesh Integration mesh, of library \code{secr} class 'mask': an Mx2 matrix.
#' @param dist KxM matrix of distances between the K detectors in traps(capthist) and the M mesh points.
#' 
#' @export
scr.negloglik=function(pars, capthist, mesh, dist=NULL) {

  M=dim(mesh)[1] # number of mesh points
  a = rep(attributes(mesh)$area,M) # area of each mesh cell (vector of length M)
  
  D = rep(exp(pars["D"]),M) # Density at each mesh point (vector of length M)
  
  # Calculate p..(s) and log(P_i(s)) at each mesh point
  scr = calc.p.Pis(capthist, mesh, pars, dist=dist)

# E[n(s)] = \int a(s) * D(s) * p..(s) ds:
  En.s = a * D * scr$p..s # E[n(s)] = a(s) * D(s) * p..(s)
  En = sum(En.s)
  
  # calculate a(s) * D(s) * P_i(s) at all mesh points
  aDPi.s=t(a * D * t(scr$Pi.s)) # transposing to multiply by row, not colum
  Pch.i=apply(aDPi.s,1,sum) # Prob of observed capture histories (vector of length n)
  
  # log-likelihood
  n = length(Pch.i)
  l = sum(log(Pch.i)) - En - lfactorial(n)
  
  return(loglik = -l)
}


#' @title Calculates p..(s) and P_i(s) for SCR likelihood calculation.
#'
#' @description  Calculates p..(s) (a vector of length n) and P_i(s) an nxM matrix of probabilities for 
#' SCR likelihood calculation.
#'
#' @param pars Vector of paramters, named "g0", "sigma" and "D".
#' @param capthist Capture history object, or library \code{secr} class 'capthist.
#' @param mesh Integration mesh, of library \code{secr} class 'mask': an Mx2 matrix.
#' @param dist KxM matrix of distances between the K detectors in traps(capthist) and the M mesh points.
#' 
#' @return A list with these elements: 
#' \itemize{
#'  \item{p..s}{Probability of detection by any detector, for activity centres at each mesh point 
#'  (vecotr of length M)}
#'  \item{Pi.s}{For each of n detected individuals, the rows of this matrix contain the probability of 
#'  obtaining the observed capture history, for all M possible activity centre locations (mesh points) 
#'  (nxM matrix).}
#'  }
#' 
#' @export
calc.p.Pis=function(capthist, mesh, pars, dist) {
  
  # calculate distances from each detector to each mesh point
  if(is.null(dist)) dist=distances(traps(capthist),mesh)
  
  # Proximity detector with count data
  if (detector(traps(capthist))=="count") {
    # Probabilities for each detector and each mesh point (assumed same at all times)
    g0=exp(pars["g0"])
    sigma=exp(pars["sigma"])
    er <- g0 * exp(-dist^2 / 2 / sigma^2) # Poisson encounter rate at all M distances from each detector
    
    # probability of being caught at least once if at mesh vertex s of M
    p.ts <- 1 - exp(- apply(er,2,sum)) # over a single occasion
    p..s <- 1 - exp(- apply(er * dim(capthist)[2],2,sum)) # over all occasions
    
    # calculate log(P_i(s)): nxM matrix
    #comb.capthist=apply(capthist,c(1,3),sum)
    #log.Pi.s = counts1.log.Pi.s(comb.capthist,er)
    log.Pit.s=counts.log.Pi.s(capthist,er) # log(P_{it}(s))
    log.Pi.s = apply(log.Pit.s,c(1,3),sum) # P_i(s) = \sum_t log(P_{it}(s)) over occasions
  }
  # Multi-catch trap
  if (detector(traps(capthist))=="multi") {
    # Probabilities for each detector and each mesh point (assumed same at all times)
    g0=plogis(pars["g0"])
    sigma=exp(pars["sigma"])
    gtk <- g0 * exp(-dist^2 / 2 / sigma^2)
    log.gtk <- log(gtk) # for use below
    log.gtk1 <- log(1-gtk) # for use below
    
    # probability of being caught at least once if at mesh vertex s of M
    p..s <- 1 - apply(1-gtk, 2, prod) ^ dim(capthist)[2] 
    
    # calculate log(P_i(s)): nxM matrix
    log.Pi.s=multi.log.Pi.s(capthist,log.gtk,log.gtk1)
  }
  
  # put all this crap in a list to pass
  return(list(p..s=p..s, Pi.s=exp(log.Pi.s)))
}



#' @title Calculates distances between two sets of points in the plane
#' 
#' @details Calculates distances between points in X and Y, whic are 
#' 2-column matrices of coordinates (x- and y-coordinates in the two columns).
#' Returns a matrix with dimensions A=dim(X)[1] and B=dim(Y)[1], containing the distances.
#' 
#' @param X is an A by 2 matrix.
#' @param Y is a B by 2 matrix.
#'
#' @references This function written by Murray Efford
#' 
distances <- function (X, Y) {
  onerow <- function (xy) {
    d <- function(xy2) sqrt(sum((xy2 - xy)^2))
    apply(Y, 1, d)
  }
  t(apply(X, 1, onerow))
}


# MULTI-CATCH detectors (detector(traps(capthist))=="multi")
# ---------------------
# Calculates log(P_i(s)) for a single individual for all s in mesh
# wi is ith row of multi-catch trap capture history (ith individual's capture history)
# of length nt, and with detector number in non-zero elements
multi.log.Pi.si = function(wi,log.gtk,log.gtk1) {
  delta=rep(0,dim(log.gtk)[1])
  delta[wi[wi>0]]=1
  log.Pi.si <- delta %*% log.gtk  + (1-delta) %*% log.gtk1 # log(\prod_k g_itk(s))
  return(log.Pi.si)
}

# Calculates log(P_i(s)) for all individuals for all s in mesh
# capthist is n by nt multi-catch trap capture history matrix
# log.gtk and log.gtk1 are K by M matrices
# returns an n by M matrix
multi.log.Pi.s = function(capthist1,log.gtk,log.gtk1) {
  log.Pi.s=apply(capthist1, 1, multi.log.Pi.si,log.gtk=log.gtk,log.gtk1=log.gtk1)
  return(t(log.Pi.s))
}


# COUNT detectors (detector(traps(capthist))=="count")
# ---------------------
# Calculates log(P_i(s)) for a single individual for all s in mesh
# wi is ith row of count capture history (ith individual's capture history)
# of length K, containing capture frequencies at each detector
# er is K by M matrix containing encounter rates for each trap at each mesh point
count.log.Pi.si = function(wi,er) {
  one=rep(1,length(wi))
  log.Pi.si = wi %*% log(er)  - one %*% er - sum(lfactorial(wi)) # log(\prod_k Poisson(n_{ks}))
  return(log.Pi.si)
}

# Calculates log(P_i(s)) for one occasion for all individuals for all s in mesh
# ch1occ is n by K matrix of counts
# returns an n by M matrix
count.log.Pi.s = function(capthist1,er) {
  log.Pi.s=apply(capthist1, 1, count.log.Pi.si,er=er)
  return(t(log.Pi.s))
}

# Calculates log(P_i(s)) for all occasions for all individuals for all s in mesh
# capthist is n by nt by K matrix of counts
# returns an n by M matrix
counts.log.Pi.s = function(capthist,er) {
  n=dim(capthist)[1]
  nt=dim(capthist)[2]
  M=dim(er)[2]
  log.Pit.s=array(dim=c(n,nt,M))
  for(i in 1:nt) { # calculate log(P_i(s)) for each occasion:
    log.Pit.s[,i,]=count.log.Pi.s(capthist[,i,],er)
  }
  return(log.Pit.s)
}


