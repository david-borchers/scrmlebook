addRaster2Mask = function(ras,mask,maskCRS,covname=NULL) {
  # To transform mask coordinates onto sgdf CRS, need to make sp object with mask CRS
  maskspdf = SpatialPoints(data.frame(x=mask$x,y=mask$y),proj4string=maskCRS)
  # Then transform onto sgdf CRS
  tfmaskspdf = spTransform(maskspdf,CRS(proj4string(ras)))
  cov = extract(ras,tfmaskspdf)
  covariates(mask)$placeholdercovariate = cov
  names(covariates(mask))[length(names(covariates(mask)))] = covname[1]
  return(mask)
}


addSGDF2Mask = function(sgdf,mask,maskCRS,covname=NULL) {
  # To transform mask coordinates onto sgdf CRS, need to make sp object with mask CRS
  maskspdf = SpatialPoints(data.frame(x=mask$x,y=mask$y),proj4string=maskCRS)
  # Then transform onto sgdf CRS
  tfmaskspdf = spTransform(maskspdf,CRS(proj4string(sgdf)))
  cov = over(sgdf,tfmaskspdf)
  covariates(mask)$placeholdercovariate = cov
  names(covariates(mask))[length(names(covariates(mask)))] = covname[1]
  return(mask)
}
