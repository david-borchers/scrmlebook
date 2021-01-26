#' @name example-data
#' @title An example data set
#' @docType data
#' @description Contains some example data.
#'
#' @usage data("example-data")
NULL

#' @name leopard-sim
#' @title Simulated Kruger Park leopard data
#' @docType data
#' @description Contains these objects created using package
#'     \code{secr} (see that package for data formats): \describe{
#'     \item{\code{kruger.capt}:}{A capture history object (class
#'     \code{capthist}) containing the numbers of captures of each
#'     individual at each of 62 detectors.}
#'     \item{\code{kruger.capt.bin}:}{A capture history object (class
#'     \code{capthist}) containing binary indicators of capture/not of
#'     each individual at each of 62 detectors.}
#'     \item{\code{kruger.cams}:}{A traps object (class \code{traps})
#'     containing the numbers of captures for each individual at each
#'     of 62 detectors of type 'count'
#'     (\code{detector(traps(capthist))=="count"}). This object is
#'     contained in \code{traps(kruger.capt)}.}
#'     \item{\code{kruger.cams.bin}:}{A traps object (class
#'     \code{traps}) containing binary indicators of capture/not for
#'     each individual at each of 62 detectors of type 'proximity'
#'     (\code{detector(traps(capthist))=="proximity"}). This object is
#'     contained in \code{traps(kruger.capt.bin)}.}
#'     \item{\code{kruger.mask}:}{A mask object (class \code{mask})
#'     containing a mask suitable for use with the above objects.}
#'     \item{\code{kruger.popn}:}{A population object, including
#'     information about the entire simulated population (i.e., not
#'     just the detected individuals).}
#'
#'  }
#'  All of the above have these covariates attached:
#'  \describe{
#'    \item{\code{sname}:}{Name of the array. A factor with levels 'Malelane' or 'Pretoriuskop'.}
#'    \item{\code{habitat.cov}:}{Covariate indexing habitat suitability (a continuous 
#'    variable).}
#'    \item{\code{landscape}:}{Habitat class. A factor with levels 'Pretoriuskop Sourveld',
#'    'Mixed Bushwillow Woodlands', 'Malelane Mountain Bushveld', 'Thorn Veld'.}
#'    \item{\code{water}:}{A measure of annual water. A continuous variable.}
#'    \item{\code{high.water}:}{A binary variable which is 1 where \code{water} is above a 
#'    threshold value.}
#'    \item{\code{dist.to.water}:}{Distance to closest \code{high.water} value of 1.}
#'  }
#'  
#' @usage data('leopard-sim')
#' @references 
#' Maputla, N.W. 2014. Drivers of leopard population dynamics in the Kruger National Park, 
#' South Africa. PhD Thesis, University of Pretoria, Pretoria, RSA.
#' @examples
#'  data('leopard-sim')
#'  plotcovariate(kruger.mask,covariate="habitat.cov",asp=1,bty="n")
#'  plot(kruger.cams,add=TRUE)
NULL


#' @name wolverines
#' @title Wolverine baited camera trap dataset
#' @docType data
#' @description Data are from Royle et al. (2011), from a camera trap study of wolverines in shouteastern Alaska in 
#' 2008. Data are in a list that contains the following objects:
#' \describe{
#'     \item{\code{w.ch}:}{A capture history object (class
#'     \code{capthist}) containing the binary capture data for each
#'     of 21 individuals at each of 37 detectors on each of 165 occasions.
#'     The object \code{traps(w.ch)} is an object of class \code{traps}, that has
#'     covariates \code{elevation} and \code{stdelev}, containing 
#'     the elevation and the standardised elevation at each trap (standardised
#'     such that the elevations in \code{w.mesh} have mean zero and variance 1).}
#'     \item{\code{w.ch.count}:}{A capture history object (class
#'     \code{capthist}) containing the counts of captures of each
#'     of 21 individuals at each of 37 detectors across all 165 occasions.
#'     The object \code{traps(w.ch.count)} is an object of class \code{traps}, that has
#'     covariates \code{elevation} and \code{stdelev}, containing 
#'     the elevation and the standardised elevation at each trap (standardised
#'     such that the elevations in \code{w.mesh} have mean zero and variance 1).}
#'     \item{\code{w.mesh}:}{A capture history object (class
#'     \code{mask}) containing a mesh that extends 40km beyond the traps. 
#'     The mesh has covariates \code{elevation} and \code{stdelev}, containing 
#'     the elevation and the standardised elevation at each mesh point (standardised
#'     such that the elevations in \code{w.mesh} have mean zero and variance 1).}
#' }
#'
#' @usage data("wolverines")
#' @references 
#' Royle, J. A. and Kagoun, A. J. and Gardner, B. and Valkenburg, P. and Lowell, R. E. (2011) 
#' Density estimation in a wolverine population using spatial capture-recapture models. 
#' \textit{Journal of Wildlife Management} \textbf{75}: 604-611.
NULL


#' @name wolverinespatial
#' @title Spatial objects to go with \code{wolerines} data
#' @docType data
#' @description Data are from Royle et al. (2011), from a camera trap study of wolverines in shouteastern Alaska in 
#' 2008. Data are in a list that contains the following objects:
#' \describe{
#'     \item{\code{w.elev.spdf}:}{A \code{SpatialPixelsDataFrame}) containing the
#'     elevation data in the vicininty of the woleverine traps.}
#'     \item{\code{w.elev.spgf}:}{A \code{SpatialGridDataFrame}) containing the
#'     elevation data in the vicininty of the woleverine traps.}
#'     \item{\code{w.cl}:}{A \code{SpatialLinesDataFrame}) containing the
#'     Alaska coastline and land border in the vicininty of the woleverine traps.}
#'     \item{\code{w.traps.spdf}:}{A \code{SpatialPixelsDataFrame}) containing the 
#'     trap locations.}
#' }
#'
#' @usage data("wolverinespatial")
NULL
