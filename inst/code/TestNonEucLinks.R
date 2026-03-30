# Get the data
#dat36 = readRDS("sundarbans36.Rds")
dat = readRDS("./analysis/data/sundarbans.Rds")
ch = dat$capthist
cams = dat$cams
mask=dat$mask

# Geometric mean (identity link formulation)
# ------------------------------------------
# need to use link = list(noneuc="identity")
geommean = function(x) exp(mean(x)) 

geomLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = geommean, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

# Geometric mean (log link formulation)
# ------------------------------------------
# need to use (default) link = list(noneuc="log")
geommean.log = function(x) exp(mean(log(x))) 

geomLCdist.log <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = geommean.log, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

# SIGNED geometric mean (identity link)
# -------------------------------------
# need to use link = list(noneuc="identity")
# positive parameter means moving from higher to lower covariate is easier (more conductive)
# and effect is stronger at higher covariate values
signedgeommean = function(x) exp(sign(diff(x))*mean(x)) 

signedgeomLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = signedgeommean, directions = 16,symm=FALSE)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

# SIGN alone
# ----------
# need to use link = list(noneuc="identity")
# positive parameter means moving from higher to lower covariate is easier (more conductive)
# and effect is same at all covariate values (i.e. depends only on difference)
signed = function(x) exp(sign(diff(x))) 

signedLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = signed, directions = 16,symm=FALSE)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}



fit.Geomid <-secr.fit(ch, detectfn="HHN", mask=mask,
                  model=list(D~1, lambda0~1, sigma~1, noneuc~water-1), 
                  details = list(userdist = geomLCdist),
                  start=list(D=exp(-8.119033),lambda0=exp(-3.226995),sigma=exp(8.112891),noneuc=1),
                  link = list(noneuc="identity"))

fit.signedGeomid <-secr.fit(ch, detectfn="HHN", mask=mask,
                        model=list(D~1, lambda0~1, sigma~1, noneuc~water-1), 
                        details = list(userdist = signedgeomLCdist),
                        start=list(D=exp(-9.8),lambda0=exp(-2.1),sigma=exp(8.8),noneuc=0),
                        link = list(noneuc="identity"))

coefficients(fit.Geomid)
coefficients(fit.signedGeomid)
AIC(fit.Geomid,fit.signedGeomid)

fit.signedid <-secr.fit(ch, detectfn="HHN", mask=mask,
                        model=list(D~1, lambda0~1, sigma~1, noneuc~water-1), 
                        details = list(userdist = signedLCdist),
                        start=list(D=exp(-8.119033),lambda0=exp(-3.226995),sigma=exp(8.112891),noneuc=0),
                        link = list(noneuc="identity"))
coefficients(fit.Geomid)
coefficients(fit.signedGeomid)
coefficients(fit.signedid)
AIC(fit.Geomid,fit.signedGeomid,fit.signedid)

fit.log <-secr.fit(ch, detectfn="HHN", mask=mask,
                  model=list(D~1, lambda0~1, sigma~1, noneuc~water-1), 
                  details = list(userdist = geomLCdist.log),
                  start=list(D=exp(-8.119033),lambda0=exp(-3.226995),sigma=exp(8.112891)),
                  link = list(noneuc="log"))

coefficients(fit.id)
coefficients(fit.log)
region.N(fit.id)
region.N(fit.log)

covariates(mask)$water[1:10]
noneuc.id = predictDsurface(fit.id,parameter="noneuc")
covariates(noneuc.id)$noneuc.0[1:10]
noneuc.log = predictDsurface(fit.log,parameter="noneuc")
covariates(noneuc.log)$noneuc.0[1:10]

fit.id.neneg <-secr.fit(ch, detectfn="HHN", mask=mask,
                        model=list(D~1, lambda0~1, sigma~1, noneuc~water-1), 
                        details = list(userdist = geomLCdist),
                        start=list(D=exp(-8.119033),lambda0=exp(-3.226995),sigma=exp(8.112891),noneuc=-4.5),
                        link = list(noneuc="identity"))

fit.log.neneg <-secr.fit(ch, detectfn="HHN", mask=mask,
                         model=list(D~1, lambda0~1, sigma~1, noneuc~water-1), 
                         details = list(userdist = geomLCdist.log),
                         start=list(D=exp(-8.119033),lambda0=exp(-3.226995),sigma=exp(8.112891)),
                         link = list(noneuc="log"))

