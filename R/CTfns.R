######################################################################
#Misc fns for CT SCR
######################################################################

#' @title Create gam object using mgcv.
#' @description  Creates a GAM object from supplied x and y values that specify a detection hazard.
#' Fits a cyclic cubic regression spline with the specified d.o.f. 
#'    
#' @return Returns the gam object, the spline coefficients, the fitted glm object, and the basis 
#' functions corresponding to both the vector of prediction times and the vector of cycle times.  
#' 
#' @param time.vec_ Vector of x values for the detection hazard cycle. Usually this vector is fairly short
#' since it requires a matching number of haz.vec_ values in order to specify the desired shape.
#' @param haz.vec_ Vector of values for the hazard corresponding to the different hazard cycle times.
#' The values supplied will determine the shape of the detection hazard.
#' @param GAM.k_ The specified degrees of freedom for the regression spline. This number corresponds to what gets passed to mgcv.
#'  Note that one degree of freedom is lost for a cyclic cubic spline and hence the actual number of parameters that get estimated is one less than the supplied value.
#' @param pred.mesh_ Specifies the time vector used for prediction / plotting. This vector is different 
#' to the time.vec_ argument in order to produce smooth plots. 
#' @export
create.GAM <- function(time.vec_, haz.vec_, GAM.k_ = 4, pred.mesh_ = NULL){
  require(mgcv)
  data = data.frame(y = haz.vec_, x = time.vec_)
  y <- haz.vec_
  formula = y ~ s(x, k = GAM.k_, bs='cc')
  smoothsetup = gam(formula, data = data, fit = FALSE)
  X = smoothsetup$X
  var.names = gsub("\\(", "", smoothsetup$term.names)
  var.names = gsub("\\)", "", var.names)
  colnames(X) = var.names
  newformula = formula(paste("y ~", paste(c("-1", colnames(X)), collapse = " + ")))
  fit.glm <- glm(newformula, data = as.data.frame(X), family='poisson')
  class(smoothsetup) = "gam"
  smoothsetup$coefficients = rep(NA, ncol(X))
  if (is.null(pred.mesh_)){pred.mesh_ <- time.vec_}
  predX = mgcv::predict.gam(smoothsetup, data.frame(x = pred.mesh_), 'lpmatrix')
  return(list('GAM obj' = smoothsetup, 'GAM ests ' = fit.glm$coefficients, 'fit.obj' = fit.glm, 'Pred basis' = predX, 'Cycle basis' = X))
}

#' @title Predict detection hazard
#' @description Uses a gam object to predict encounter rate / hazard of detection that correspond to the x values supplied.
#'    
#' @return Returns the hazard values.
#' 
#' @param time.vec_ Vector of x values for the detection hazard cycle. 
#' @param gam.obj_ A gam object created using package mgcv. Required in order to extract the basis functions.
#' @param fit.obj_ A fitted glm object. The spline coefficients are extracted from the glm object when simulating.
#' @param spline.coeff_ A vector of coefficients to be used in the regression spline.
#' 
#' @export
pred.Y <- function(time.vec_, gam.obj_, fit.obj_ = NULL, spline.coeff_ = NULL){
  newdata = data.frame(x = time.vec_)
  class(gam.obj_) = "gam"
  gam.obj_$coefficients = rep(NA, ncol(gam.obj_$X))
  predX = mgcv::predict.gam(gam.obj_, newdata, 'lpmatrix')
  if (is.null(spline.coeff_)) {spline.coeff_ <- fit.obj_$coefficients}
  preds <-  exp(predX %*% spline.coeff_)
  return(preds)
}

#' @title Calcluate distances.
#' @description Calculates euclidean distances between two sets of coordinates. 
#'    
#' @return Returns the matrix of distances.
#' 
#' @param X, Y Matrices of coordinates.  
#' 
#' @export
distances <- function (X, Y) {
  ## X and Y are 2-column matrices of coordinates
  onerow <- function (xy) {
    d <- function(xy2) {
      sqrt(sum((xy2 - xy)^2))
    }
    apply(Y, 1, d)
  }
  t(apply(X, 1, onerow))
}

#' @title AICc
#' @description  Calculates the AICc for a given model or set of models.  
#'    
#' @return Returns the AICc values.  
#' 
#' @param mod.objs_ A list of the model objects.
#' @param decimals_ The number of decimals to round the value to.
#' @param n.inds_ The number of detected individuals in the data set.
#'  
#' @export
Calc.AICc <- function(mod.objs_, decimals_ = 2, n.inds_){
  if (is.list(mod.objs_)){
    n.mods <- length(mod.objs_)
    output1 <- NULL
    output2 <- NULL
    for (i in 1:n.mods){
      K = length(mod.objs_[[i]]$estimate)
      AIC <- 2 * mod.objs_[[i]]$minimum + 2*K
      AICc <- AIC + (2*K*(K+1))/(n.inds_ - K - 1)
      output1[i] <- round(AICc, decimals_)
      output2[i] <- round(AIC, decimals_)
    }
    return(list("AICc" = output1, "AIC" = output2))
  } else {print("Model objects must be a list")}
}

#' @title Calculate the variance for the time-dependent expected encounter rate.
#' @description  Calculates the variance for the estimated time-dependent ER that is used to plot confidence intervals. The 
#' time component of the ER is a simple linear combination of spline parameters, whereas the full encounter rate function is 
#' a non-linear function of the spline and the detection function parameters. Both the Delta method and a parametric bootstrapping 
#' approach is used to estimate the variance of the full encounter rate function. 
#'    
#' @return Returns a matrix of lower and upper limits (for a given alpha value) for the estimated ER or sigma(t) / lambda0(t) at a 
#' vector of time points. If the bootstrap approach is used, the estimated ERs from each bootstrap replication are also returned whereas the
#' if the Delta method is chosen the estimated standard errors (on the ER scale) are also returned. 
#' 
#' @param par.ests The estimated parameters from the chosen model (excluding Density). The order of the parameters (on the link scale)
#' needs to be: lambda0, sigma, and then the spline parameters.
#' @param num.spline.pars The number of estimated spline parameters.
#' @param vcov.mat The estimated variance-covariance matrix, from inverting the hessian. The matrix should exclude sigma if the Delta 
#' method is chosen for the lambda0(t) function. 
#' @param type Specifies the model parameterisation as either "lambda.t" or "sigma.t".
#' @param method This argument specifies the method that is used to approximate the variance, the Delta method (`delta') or a 
#' parametric bootstrap approach (`boot').
#' @param Breps Sets the number of bootstrap replications.
#' @param alpha Sets the desired level of confidence for the intervals that are returned.
#' @param resolution The temporal resolution of the points.
#' @param dist A single value for distance. If this is NULL the function assumes the time-dependent part of the function.  
#' @param pred.mesh Allows the function to be plotted for any vector of points i.e. not just for a single cycle.
#'  
#' @export
Calc.ER.Variance <- function(par.ests, num.spline.pars, vcov.mat, type = "lambda.t", method = "delta", Breps = 1000, alpha = 0.05, resolution = 100, dist = NULL, pred.mesh = NULL){
  
  if(!is.element(method,c("boot","delta"))) stop("method must be `boot` (for parametric bootstrap) or `delta` (for delta method).")
  if(!is.element(type,c("lambda.t","sigma.t"))) stop("type must be `lambda.t` or `sigma.t`.")
  
  cycle.mesh = seq(0,24,length = resolution)
  cycle.mesh.width <- cycle.mesh[2]-cycle.mesh[1]
  smoothsetup.s <- create.GAM(cycle.mesh, rep(0,length(cycle.mesh)),GAM.k_ = num.spline.pars+2, pred.mesh_ = pred.mesh)
  X <- smoothsetup.s$`Pred basis`
  
  num.time.pts <- dim(X)[1]
  
  if (method =="delta"){
    sds <- numeric(num.time.pts)
    
    for (i in 1:num.time.pts){
      if (is.null(dist)){   #dist = NULL indicates the lambda0.t parameterisation
        sds[i] <- sqrt(t(X[i,]) %*% vcov.mat %*% X[i,])
      } else {
        #give list of partial derivatives, the X matrix will be of suitable dimension
        if (type=="lambda.t"){
          pds <- c(1, dist^2*exp(par.ests[2]) * exp(par.ests[2])^(-3), X[i,-1])
        } else {
          pds <- c(1)
          for (j in 1:(num.spline.pars+1)){
            temp <- dist^2*exp(par.ests[-1] %*% X[i,]) * exp(par.ests[-1] %*% X[i,])^(-3) * X[i,j]
            pds <- c(pds, temp)
          }
        }
        sds[i] <- sqrt(t(pds) %*% vcov.mat %*% pds)
      } 
    }
    
    if (type=="lambda.t"){ 
      spline.ests <- c(par.ests[1]) 
      for (j in 1:num.spline.pars){
        temp <- par.ests[2+j]
        spline.ests <- c(spline.ests, temp)
      }
      t <- exp(X %*% spline.ests)
      ER <- t * exp(-dist^2 / (2*exp(par.ests[2])^2))
    } else { 
      spline.ests <- c(par.ests[-1])
      t <- exp(X %*% spline.ests)
      ER <- exp(par.ests[1]) * exp(-dist^2 / (2*t^2))
    }
    
    if (is.null(dist)){
      cis <- matrix(c(t * exp(qnorm((alpha/2))*sds), t * exp(qnorm(1-(alpha/2))*sds)), ncol = 2)
    } else {
      cis <- matrix(c(ER * exp(qnorm((alpha/2))*sds), ER * exp(qnorm(1-(alpha/2))*sds)), ncol = 2)
    }
    
    return(list("Confidence intervals" = cis, "Approximate sds" = exp(sds)))
  } else {
    
    if(is.null(par.ests)) stop("Parameter estimates needed for a parametric bootstrap.")
    
    t.boot <- matrix(rep(NA, resolution*Breps), nrow = Breps)
    ER.boot <- matrix(rep(NA, resolution*Breps), nrow = Breps)
    
    bpars = mvrnorm(Breps, par.ests, vcov.mat)
    
    for (i in 1:Breps){
      if (type=="lambda.t"){
        spline.pars <- c(bpars[i,1])
      } else {
        spline.pars <- c(bpars[i,2])
      }
      for (j in 1:num.spline.pars){
        temp <- bpars[i,2+j]
        spline.pars <- c(spline.pars, temp)
      }
      
      t.boot[i,] <- exp(X %*% spline.pars)  #either lambda0(t) or sigma(t)
      
      if (!is.null(dist)){
        if (type=="lambda.t"){
          ER.boot[i,] <- exp(-dist^2 / (2*exp(bpars[i,2])^2)) * t.boot[i,]  #note that lambda0 par already in t.boot for l(t) param.
        } else {
          ER.boot[i,] <- exp(bpars[i,1])*exp(-dist^2 / (2*t.boot[i,]^2))
        }
      }
    }
    
    if (!is.null(dist)){
      cis <- t(apply(ER.boot,2,quantile,probs=c(alpha/2,1-alpha/2)))
      return(list("Confidence intervals" = cis, "ER fns" = ER.boot, "(t) fns" = t.boot))
    } else {
      cis <- t(apply(t.boot,2,quantile,probs=c(alpha/2,1-alpha/2)))
      return(list("Confidence intervals" = cis, "(t) fns" = t.boot))
    }
  }
}


#' @title Plot cyclical detection hazard.
#' @description Plots the detection hazard from a set of spline coefficients for a specified time period. The function has different ways 
#' of plotting different types of hazards, including a 3D plot for the dependent hazard.
#'    
#' @return Returns a plot.
#' 
#' @param cycle.hours_ Specifies the vector of times for one cycle of the cyclical hazard. 
#' @param time.vec_ Specifies the vector of times for the plotted hazard.
#' @param dist.vec_ A vector of distance values that is used to generate the 3D plot. 
#' @param mesh.size_ Determines the coarseness of the plotted lines.
#' @param det.pars_ The relevant parameters from the detection function or encounter rate function.
#' . If the 'IndU' hazard is specified these will include lambda0 and sigma, if 'IndL' 
#' then g0 replaces lambda0, and if a dependent hazard is specified ('Dep') then only
#' lambda0 is relevant. 
#' @param spline.pars_ The set of spline coefficients that determine the shape of the hazard.
#' @param vcov_ A variance-covariance matrix for the parameters that is used to get std errors for confidence intervals.
#' It must be of the appropriate dimension and exclude sigma if the Delta method is used for the lambda0(t) model. Note
#' the plotting of ci's only currently available for lambda0(t) model. 
#' @param ci.method_ Specifies the method used to calculate the variance and produce intervals, choices include "delta" or "boot".
#' Note that a plot that uses the bootstrap approach with the sigma(t) hazard will produce slightly wonky intervals. Calc.ER.Variance
#' needs to be extended for multiple distances to correct this.
#' @param haz.k_ The appropriate df. This number should correspond to the value passed to mgcv
#' and not to the number of spline coefficients.
#' @param haz.type_ The chosen parameterisation (from 'IndU', 'IndL', and 'Dep').  
#' @param dim_ The dimension argument is used when plotting the dependent hazard and determines what type of plot
#' is drawn. If 'time' is specified, the plot will depict the hazard over the time dimension
#' for a given distance (specified with the dist_ argument), if 'dist' is specified the plot
#' will depict the hazard over distance for a given time (specified with the time_ argument),
#' and if nothing is specified (dim_ is null) then a plot is drawn that includes both dimensions.  
#' @param full.ER_ Specifies whether only the time component of the function is displayed, which is appropriate when 
#' interest is focused on the shape of the hazard through time. If TRUE the full encounter rate function is plotted 
#' for a given distance.
#' @param add_ Specifies whether a new plot should be called or rather if lines should
#' be added to an existing plot.
#' @param dist_ A specified value for distance. More than one value can be supplied in which case a line will be drawn for each.
#' @param time_ A specified value for time. More than one value can be supplied in which case a line will be drawn for each.
#' @param main_ The title for the plot.
#' @param xlabel_ Label for the X axis. 
#' @param ylabel_ Label for the y axias.
#' @param ylim_ Option to specfiy the scale of the Y axis.
#' @param cex.legend_ A scaling factor to control the size of the legend.
#' @param legend.show_ Specifies if a legend should appear in the plot. It is not relevant 
#' for plots showing the scaled hazard.
#' @param legend.pos_ Specifies the position of the legend.
#' @param col_ Specifies the colour to use for the plotted line
#' @param theta_ defines the viewing angle, theta sets the azimuthal direction. 
#' @param phi_ defines the viewing angle, phi sets the colatitude
#' @param rotate_ If this is set to FALSE a static 3 dimensional perspective plot is generated. If it is TRUE
#' the 3D plot can be rotated by the user  (requires the rgl library).
#' @param ydecimals_ Specifies the number of decimals to use for the labelling of the Y axis
#' 
#' @export
hazard.plot <- function(cycle.hours_ = c(0,24), time.vec_ = NULL, dist.vec_ = 0:1000, mesh.size_ = 100, det.pars_, spline.pars_, vcov_ = NULL, ci.method_ = "delta", haz.k_ = c(4), haz.type_ = 'IndU', dim_ = 'time', full.ER_ = T, add_=F, dist_ = NULL, time_ = NULL, main_="", xlabel_ = "Cycle Hour", ylabel_ = "Encounter Rate", ylim_ = NULL, cex.legend_ = 1, legend.show_ = T, legend.pos_ = "topright", col_="black", theta_ = 30, phi_ = 15, rotate_ = F, ydecimals_ = 3){
  cycle.mesh <- seq(cycle.hours_[1], cycle.hours_[2], length = 100)
  cycle.mesh.width = cycle.mesh[2]-cycle.mesh[1]
  if (is.null(time.vec_)) {time.vec_ <- cycle.hours_[1]:cycle.hours_[2]}
  mesh <- seq(min(time.vec_), max(time.vec_), length = mesh.size_)
  mesh.width = mesh[2]-mesh[1]
  at1 <- seq(min(mesh), max(mesh), by = 1)
  xlim = c(min(mesh), max(mesh))  #both at1 and xlim sometimes overwritten below
  
  #det fn pars
  if (haz.type_ == 'Dep'){
    lambda0 = exp(det.pars_[1])
  } else {
    sigma = exp(det.pars_[2])
    if (haz.type_ == 'IndL'){
      g0 = invlogit(det.pars_[1])
    } else {
      lambda0 = exp(det.pars_[1])
    }
  }
  
  #create GAM obj with full haz cycle
  GAM.obj <- create.GAM(time.vec_ = cycle.mesh, haz.vec_ = rep(0,length(cycle.mesh)), GAM.k_ = haz.k_, pred.mesh_ = mesh)
  
  if (haz.type_ == 'Dep'){
    sigma.cycle <- exp(GAM.obj[[4]] %*% c(spline.pars_))  #spline.pars includes B0
    
    if (is.null(dim_)){
      #function to use with outer
      Haz.calc <- function(x_,y_, lambda0_){
        lambda0_*exp(-x_^2 / (2*y_^2))
      }
      haz.z <- outer(dist.vec_,sigma.cycle, Haz.calc, lambda0_=lambda0)
      if (rotate_ == F){
        persp(dist.vec_, mesh, haz.z[,,1], xlab = "Distance (m)", ylab="Time (h)", zlab = "Encounter rate", main = main_, theta = theta_, phi = phi_, expand = 0.65, col = "#1283D4", ticktype = "detailed", nticks = 4, cex.axis = 1, cex.lab = 1.5)
      } else {
        library(rgl)
        plot3d(dist.vec_, sigma.cycle, haz.z[,,1], 
               type="n", xlim=c(min(dist.vec_), max(dist.vec_)), ylim=c(min(time.vec_), max(time.vec_)), zlim=c(0,max(haz.z[,,1])), # specify axis ranges, labels, etc.
               xlab="Distance", ylab="Time", zlab="Encounter Rate")
        
        surface3d(dist.vec_, mesh, haz.z, col="#1283D4", alpha=1)
      }
    } else {
      if (dim_ == 'time'){
        if (is.null(dist_)) {stop("Need a value for distance") }
        lambda.dt <- matrix(NA, nrow = length(dist_), ncol = length(sigma.cycle))
        lambda.dt.ci <- array(dim = c(length(sigma.cycle),2, length(dist_)))
        n.lines <- length(dist_)
        legend.desc <- as.character(n.lines)
        
        for (i in 1:n.lines){
          lambda.dt[i,] <- lambda0*exp(-dist_[i]^2 / (2*sigma.cycle^2))
          legend.desc[i] <- paste(dist_[i],"m")
          if (!is.null(vcov_)){
            lambda.dt.ci[,,i] <- Calc.ER.Variance(par.ests = c(det.pars_,spline.pars_), num.spline.pars = haz.k_-2, vcov.mat = vcov_, type = "sigma.t", method = ci.method_, dist = dist_[i], pred.mesh = mesh)[[1]]
          }
        }
        
        x.label <- "Cycle Hour (h)"
        xlim = c(min(mesh), max(mesh))
        
      } else {
        if (is.null(time_)) {stop("Need a value for time") }
        if (length(dist.vec_)!= length(mesh) ) {dist.vec_ = seq(dist.vec_[1], max(dist.vec_), length = length(mesh)) } #ensures dist.vec of same dimension as sigma.cycle
        lambda.dt <- matrix(NA, nrow = length(time_), ncol = length(sigma.cycle))
        lambda.dt.ci <- array(dim = c(length(dist.vec_),2, length(time_)))
        ci.d <- matrix(NA, nrow = length(sigma.cycle), ncol = 2)
        n.lines <- length(time_)
        legend.desc <- as.character(n.lines)
        
        for (i in 1:n.lines){
          legend.desc[i] <- paste(time_[i],"h")
          sigma.pos <- which.min(abs(mesh - time_[i]))
          lambda.dt[i,] <- lambda0*exp(-dist.vec_^2 / (2*sigma.cycle[sigma.pos]^2)) #ER for a given time across dist.vec
          
          if (!is.null(vcov_)){
            for (j in 1:length(dist.vec_)){ #needed because Calc.ER.Var works with a single distance
              d <- dist.vec_[j]
              ci <- Calc.ER.Variance(par.ests = c(det.pars_,spline.pars_), num.spline.pars = haz.k_-2, vcov.mat = vcov_, type = "sigma.t", method = ci.method_, dist = d, pred.mesh = mesh)[[1]]
              ci.d[j,] <- ci[sigma.pos,]
            }
            lambda.dt.ci[,,i] <- ci.d
          }
        }
        
        x.label <- "Distance (m)"
        mesh <- dist.vec_  #overwrites mesh from earlier
        xlim = c(min(mesh), max(mesh))
        at1 <- seq(min(mesh), max(mesh), length = 8)
      }
      
      if (add_==F){
        if (is.null(ylim_)) {ylim=c(0, 1.25*max(lambda.dt))} else {ylim = ylim_}
        if (!is.null(vcov_)){
          at2 <- seq(from = ylim[1], to = round(max(lambda.dt.ci[[1]]),3), length = 5)
        } else {
          at2 <- seq(from = ylim[1], to = ylim[2], length = 5)
        }
        plot(x=xlim, y=ylim, type="n", xaxt='n', yaxt='n', xlab = x.label, ylab="Expected Encounter Rate", main = main_, las=2, cex.lab = 1.5, cex.axis = 1.5) 
        axis(side = 1, at = at1, labels = T) 
        axis(side = 2, at = round(at2, ydecimals_), labels = T)
        for (j in 1:n.lines){
          lines(lambda.dt[j,] ~ mesh, lwd=2, col=col_[j]) 
          if (!is.null(vcov_)){
            lines(lambda.dt.ci[,1,j] ~ mesh, lty = 3, lwd = 2, col = col_[j])
            lines(lambda.dt.ci[,2,j] ~ mesh, lty = 3, lwd = 2, col = col_[j])
          }
        }
        legend.base <- legend(legend.pos_, legend=legend.desc, col=col_, lty = rep(1,n.lines),  cex = cex.legend_, plot = legend.show_, bty = "n")
      } else {
        for (j in 1:n.lines){
          lines(lambda.dt[j,] ~ mesh, lwd=2, col=col_[j]) 
          lines(lambda.dt.ci[,1,j] ~ mesh, lty = 3, lwd = 2, col = col_[j])
          lines(lambda.dt.ci[,2,j] ~ mesh, lty = 3, lwd = 2, col = col_[j])
        }
        legend.base <- legend(legend.pos_, legend=legend.desc, col=col_, lty = rep(1,n.lines), cex = cex.legend_, plot=F)
        legend(x = legend.base$rect$left, y = legend.base$rect$top - 
                 legend.base$rect$h, legend=legend.desc, col = col_,lty = rep(1,n.lines), cex = cex.legend_, plot = legend.show_, bty = "n")
      }
    }
  } else {
    lambda.t.cycle <- exp(GAM.obj[[5]] %*% c(det.pars_[1],spline.pars_))
    lambda.t.full <- exp(GAM.obj[[4]] %*% c(det.pars_[1],spline.pars_))
    cycle.int <- sum(lambda.t.cycle)*cycle.mesh.width
    
    if (full.ER_ != T){ #this specifies a plot for lambdat vs for the full ER fn
      
      if (is.null(ylim_)) {ylim=c(0, round(1.25*max(lambda.t.cycle),3))} else {ylim = ylim_}
      at2 <- seq(from = 0, to = round(max(lambda.t.cycle), 3), length = 5)
      if (!is.null(vcov_)){
        lambda.t.ci <- Calc.ER.Variance(par.ests = c(det.pars_, spline.pars_), num.spline.pars = haz.k_-2, vcov.mat = vcov_, method = ci.method_, dist = NULL, pred.mesh = mesh)
        at2 <- seq(from = 0, to = round(max(lambda.t.ci[[1]]),3), length = 5)
      }
      if (add_==F){
        plot(x=xlim, y=ylim, type="n", xaxt='n', yaxt='n', xlab = xlabel_, ylab = ylabel_, main = main_, las=2, cex.lab = 1.5, cex.axis = 1.5) 
        axis(side = 1, at = at1, labels = T) 
        axis(side = 2, at = round(at2,3), labels = T)
      }
      lines(lambda.t.full ~ mesh, lwd=2, col=col_[1])
      if (!is.null(vcov_)){
        lines(lambda.t.ci[[1]][,1] ~ mesh, lty = 3, lwd = 2, col = col_[1])
        lines(lambda.t.ci[[1]][,2] ~ mesh, lty = 3, lwd = 2, col = col_[1])
      }
    } else {  #this loop is for the full ER function
      if (is.null(dist_)){stop("Need a value for distance") }
      n.lines <- length(dist_)
      lambda.d <- as.numeric(n.lines)
      legend.desc <- as.character(n.lines)
      
      if (add_==F){ #new plot requested
        if (is.null(ylim_)) {ylim=c(0, round(1.25*max(lambda.t.cycle),3))} else {ylim=ylim_}
        at2 <- seq(from = ylim[1], to = ylim[2], length = 5)
        plot(x=xlim, y=ylim, type="n", xaxt='n', yaxt='n', xlab = xlabel_, ylab = ylabel_, main = main_, las=2, cex.lab = 1.5, cex.axis = 1.5) 
        axis(side = 1, at = at1, labels = T) 
        axis(side = 2, at = round(at2, ydecimals_), labels = T)
        
        for (j in 1:n.lines){
          if (haz.type_ == 'IndL'){
            lambda.d[j] <- log(1/(1-g0*exp(-dist_[j]^2 / (2*sigma^2)))) / cycle.int
          } else {
            lambda.d[j] <- exp(-dist_[j]^2 / (2*sigma^2)) #lambda0 already included in lambda.t
          }
          
          ER <- lambda.t.full * lambda.d[j]
          
          lines(ER ~ mesh, lwd=2, col=col_[j])
          if (!is.null(vcov_)){
            ER.ci <- Calc.ER.Variance(par.ests = c(det.pars_, spline.pars_), num.spline.pars = haz.k_-2, vcov.mat = vcov_, dist = dist_[j], method = ci.method_, resolution = mesh.size_, pred.mesh = mesh)
            lines(ER.ci[[1]][,1] ~ mesh, lty = 3, lwd = 2, col = col_[j])
            lines(ER.ci[[1]][,2] ~ mesh, lty = 3, lwd = 2, col = col_[j])
          }
          H.period <- sum(ER)*mesh.width
          cap.prob <- round(1-exp(-H.period),3)
          legend.desc[j] <- paste("Dist of ", dist_[j], "\n Prob = ", cap.prob)
        }
        legend(legend.pos_, legend=legend.desc, col=col_, lty = rep(1,n.lines), y.intersp = 1.5, cex = cex.legend_, plot = legend.show_, bty = "n")
      } else {
        for (j in 1:n.lines){
          if (haz.type_ == 'IndL'){
            lambda.d[j] <- log(1/(1-g0*exp(-dist_[j]^2 / (2*sigma^2)))) / cycle.int
          } else {
            lambda.d[j] <- exp(-dist_[j]^2 / (2*sigma^2)) #lambda0 already included in lambda.t
          }
          
          ER <- lambda.t.full * lambda.d[j]
          
          lines(ER ~ mesh, lwd=2, col=col_[j])
          if (!is.null(vcov_)){
            ER.ci <- Calc.ER.Variance(par.ests = c(det.pars_, spline.pars_), num.spline.pars = haz.k_-2, vcov.mat = vcov_, dist = dist_[j], method = "delta", resolution = mesh.size_, pred.mesh = mesh)
            lines(ER.ci[[1]][,1] ~ mesh, lty = 3, lwd = 2, col = col_[j])
            lines(ER.ci[[1]][,2] ~ mesh, lty = 3, lwd = 2, col = col_[j])
          }
          H.period <- sum(ER)*mesh.width
          cap.prob <- round(1-exp(-H.period),3)
          legend.desc[j] <- paste("Dist of ", dist_[j], "\n Prob = ", cap.prob)
        }
        legend.base <- legend(legend.pos_, legend=legend.desc, col=col_, lty = rep(1,n.lines), cex = cex.legend_, plot=F)
        legend(x = legend.base$rect$left, y = legend.base$rect$top - 
                 legend.base$rect$h, legend=legend.desc, col = col_,lty = rep(1,n.lines), cex = cex.legend_, plot = legend.show_, bty = "n")
      }
    }
  }
}

#' @title Calculates the probability of detection from a chosen activity centre, for a given time point and for a given duration. 
#' @description The function calculates the detection probability from a chosen pixel for a given time point and for a 
#' given duration. The probability depends on the mask of the study area, the covariate attached to the mask, and a given 
#' conductance parameter.
#'    
#' @return Returns a vector of capture probabilities associated with each pixel in the mask.
#' 
#' @param pixel_ The two dimensional coordinate that corresponds to the pixel of interest.
#' @param mask.obj_ A mask object that needs to have one numeric covariate attached to it.
#' @param beta1_ The conductance parameter.
#' @param duration_ The length of time that the detection probability is caluclated for.
#' @param time_ A supplied time in the 24 hour cycle.
#' @param lambda0_ A value for the base encounter rate at a distance of zero/
#' @param spline.coeff_ A vector of spline coefficients that determine the shape of the temporal encounter rate function.
#' @param sigma_ A value for the sigma parameter in the encounter rate function. This needs to be supplied for the lambda(t) 
#' function. If it is NULL the sigma(t) form will be assumed.
#' 
#' @export
prob.detect <- function(pixel_, mask.obj_, beta1_, duration_, time_, lambda0_, spline.coeff_, sigma_ = NULL){
  #calculate LC dists for a given point
  LC <- noneuc.distances(pixel_, mask.obj_, conduct.par_ = beta1_)
  
  #create GAM object, needed to extract basis fns
  cycle.hours_=c(0,24)
  cycle.mesh <- seq(0,24, length = 100)
  fine.mesh <- seq(0,24, length = 1000)
  if (is.null(sigma_)){
    GAM.obj <- create.GAM(time.vec_ = cycle.mesh, haz.vec_ = rep(0,length(cycle.mesh)), GAM.k_ = length(spline.coeff_)+1, pred.mesh_ = fine.mesh)
    sig.t <- pred.Y(time_, GAM.obj[[1]], spline.coeff_ = spline.coeff_)
    lambda.dt <- lambda0_*exp(-as.numeric(LC)^2 / (2*sig.t^2))
  } else {
    GAM.obj <- create.GAM(time.vec_ = cycle.mesh, haz.vec_ = rep(0,length(cycle.mesh)), GAM.k_ = length(spline.coeff_)+2, pred.mesh_ = fine.mesh)
    haz.t <- pred.Y(time_, GAM.obj[[1]], spline.coeff_ = c(2,spline.coeff_))
    lambda.dt <- lambda0_*exp(-as.numeric(LC)^2 / (2*sigma_^2))*haz.t
  }
  return(list("probs" = 1 - exp(-duration_*lambda.dt), "LC dists" = as.numeric(LC)))
}

#' @title Construct the hazard vector for the full study.
#' @description Constructs the hazard vector for the full study duration from the hazard for a single cycle. Used in both the LL
#' function (CT.SCR) and in the simulation function (sim.ctsCH).
#'    
#' @return Returns a vector of detection hazards of the same dimension as the supplied study time vector.
#' 
#' @param cycle.vec_ A vector of times that correspond to a single hazard cycle.
#' @param y.vec.full_ A vector of detection hazards that correspond to a single hazard cycle.
#' @param spline.coeff_ A vector of spline coefficients.
#' @param release.hour_ Primarily applicable to SC traps when the release event defines the occasions. Can also be used to 
#' adjust the first cycle if it starts after the cycle's start point.
#' @param endpoint_ Specifies the duration of the simulation survey.
#' @param n.haz.cycles_ The number of full hazard cycles in the full study period.
#' @param haz.interval_ The duration of a single hazard cycle.
#' @param gam.obj_ GAM obj used by function pred.Y to predict detection hazards.
#' @param fit.obj_ A fitted glm object used by function pred.Y to predict detection hazards.
#' 
#' @export
study.haz.vec <- function(cycle.vec_, y.vec.full_, spline.coeff_ = NULL, release.hour_, endpoint_, n.haz.cycles_, haz.interval_, gam.obj_, fit.obj_ = NULL){
  
  #need 2 partial vecs for bits of haz cycle at start and end of release occs
  if (release.hour_ > 0){
    if (is.null(spline.coeff_)){
      y.vec.part.1 <- pred.Y(time.vec_ = cycle.vec_[cycle.vec_ > release.hour_], gam.obj_ = gam.obj_, fit.obj_ = fit.obj_)
    } else {
      y.vec.part.1 <- pred.Y(time.vec_ = cycle.vec_[cycle.vec_ > release.hour_], gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)
    }
  } else {
    y.vec.part.1 <- NULL
  }
  
  part2 <- endpoint_ %% haz.interval_
  
  if (part2 != 0) {
    if (is.null(spline.coeff_)){
      y.vec.part.2 <- pred.Y(time.vec_ = cycle.vec_[cycle.vec_ < part2], gam.obj_ = gam.obj_, fit.obj_ = fit.obj_)
    } else {
      y.vec.part.2 <- pred.Y(time.vec_ = cycle.vec_[cycle.vec_ < part2], gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)
    }
  } else {y.vec.part.2 <- NULL}
  
  y.vec.tot <- c(y.vec.part.1, rep(y.vec.full_, times = n.haz.cycles_), y.vec.part.2)
  
  excl.time.pts <- which(y.vec.part.1==tail(y.vec.full_,1))
  
  if (n.haz.cycles_ >= 1){
    for (k in 1:n.haz.cycles_){
      excl.time.pt <- length(y.vec.part.1) + length(cycle.vec_)*k
      excl.time.pts <- c(excl.time.pts, excl.time.pt)
    }
  }
  
  y.vec.tot <- y.vec.tot[-c(excl.time.pts)]
  return(y.vec.tot)
}

######################################################################
#Simulation functions
######################################################################

#' @title Simulate CT data.
#' @description  Simulates continuous-time datasets.  
#'    
#' @return Returns the capture histories with continuous times as a data frame, and with an occasion structure as a secr capthist object.  
#' 
#' @param traps.obj_ An secr traps object.
#' @param mask.obj_ An secr mask object. If not supplied make.mask will be used to create the mask. 
#' @param buffer_ Buffer radius added to traps.obj_ to create the mask (if not supplied).
#' @param haz.type_ The type of detection hazard to use in the simulation. Choices c("IndL", "IndU", "Dep") include the two 
#' independent hazard parameterisations and the dependent hazard parameterisation.
#' @param detector_type_ Specifies the type of detector. Defaults to using a "proximity" or passive detector but can also simulate
#' from an array of single-catch traps.
#' @param gam.obj_ GAM obj to simulate capture times from according to a time varying spline hazard.
#' @param fit.obj_ A fitted glm object. The spline coefficients are extracted from the object when simulating.
#' @param cycle.hrs_ The start and end of a single detection hazard cycle, defaults to c(0,24) to represent a daily cycle.
#' @param mesh.size_ Determines the resolution of the time vector. This argument determines the vector length from the seq() command.
#' @param dens.surf_ The density to use in the simulation. If a single number then the data will be simulated with a constant density.
#' The argument can also be a matrix of values of the same dimension as the mask that is supplied.  
#' @param g0_ The intercept parameter of the detection function. Supplied when appropriate (for an independent linked hazard).
#' @param lambda0_ The intercept parameter of the encounter rate function. Supplied when appropriate (for an independent unlinked 
#' or dependent hazard).
#' @param sigma_ The range parameter of the detection / encounter rate function. Will be NULL when a dependent hazard is used.
#' @param endpoint_ Specifies the duration of the simulation survey.
#' @param setup.hour_ Used when the survey starts at some point within the first cycle.
#' @param occ.interval_ The duration of an occasion. Used to construct the discrete-time secr capthist object.
#' 
#' @export
sim.ctsCH <- function (traps.obj_, mask.obj_ = NULL, buffer_ = NULL, haz.type_ = "IndL", detector.type_ = "Prox", gam.obj_ = NULL, fit.obj_ = NULL, cycle.hrs_ = c(0,24), mesh.size_ = 1000, dens.surf_, g0_ = NULL, lambda0_ = NULL, sigma_ = NULL, endpoint_, setup.hour_ = 0, occ.interval_ = 24){
  
  if (is.null(mask.obj_)) {mask.obj_ <- make.mask (traps.obj_, buffer_)}
  
  if (length(dens.surf_) > 1){
    Pop <- sim.popn (dens.surf_, core = mask.obj_, model2D ='IHP')
  } else {
    Pop <- sim.popn (D = dens.surf_, core = traps.obj_, buffer = buffer_)
  }
  
  n.inds <- dim(Pop)[1]
  n.traps <- dim(traps.obj_)[1]
  trap.names <- rownames(traps.obj_)
  n.occs <- ceiling((endpoint_ - setup.hour_) / occ.interval_)
  if (is.null(gam.obj_)){Const.Haz = T} else {Const.Haz = F}
  
  Detector <- NULL 
  Caught <- NULL  
  id <- NULL
  
  haz.interval <- cycle.hrs_[2] - cycle.hrs_[1]
  
  if (setup.hour_ %% haz.interval != 0){
    n.haz.cycles <- ((endpoint_ - setup.hour_)-(haz.interval - setup.hour_ %% haz.interval)-(endpoint_ %% haz.interval)) / haz.interval
  } else {
    n.haz.cycles <- ((endpoint_ - setup.hour_)-(endpoint_ %% haz.interval)) / haz.interval
  } 
  
  if (n.inds>0){
    act.dists <- distances(traps.obj_, Pop) #traps * n.inds
    if (dim(act.dists)[1]==1) {act.dists=t(act.dists) ; colnames(act.dists) <- 1 }    
    
    if (Const.Haz == F){
      #integrate the gam hazard, uses fitted gam obj
      #hazt is the time-specific component, equivalent to sigma vec for Dep haz
      cycle.mesh <- seq(min(cycle.hrs_), max(cycle.hrs_), length=mesh.size_)
      cycle.mesh.width = cycle.mesh[2] - cycle.mesh[1]
      
      hazt.vec.cycle <- pred.Y(time.vec_ = cycle.mesh, gam.obj_ = gam.obj_, fit.obj_ = fit.obj_)
      hazt.vec.tot <- study.haz.vec(cycle.vec_ = cycle.mesh, y.vec.full_ = hazt.vec.cycle, release.hour_ = setup.hour_, endpoint_ = endpoint_, n.haz.cycles_ = n.haz.cycles, haz.interval_ = haz.interval, gam.obj_ = gam.obj_, fit.obj_ = fit.obj_)
      study.time <- seq(setup.hour_,endpoint_, length=length(hazt.vec.tot))
      total.mesh.width <- study.time[2]-study.time[1]
      
      hazt.cycle.occ1.int <- sum(hazt.vec.tot[study.time >= setup.hour_ & study.time <= (setup.hour_ + occ.interval_)])*total.mesh.width
      hazt.total.int <- sum(hazt.vec.tot)*total.mesh.width
    } else { #else loop is for a constant hazard
      study.time <- seq(setup.hour_,endpoint_, length=mesh.size_*n.haz.cycles)
      hazt.cycle.occ1.int <- occ.interval_
      hazt.total.int <- endpoint_ - setup.hour_
    }
    
    if (haz.type_=='Dep') {
      #calculate full hazard matrix, depends on both distance and time
      #uses actual distances
      full.haz = array(NA, dim =c(n.traps,length(hazt.vec.tot),n.inds))
      dimnames(full.haz) <- list(trap.names,NULL,seq(1,n.inds))
      for (i in 1:n.inds){
        for (j in 1:n.traps){
          full.haz[j,,i] <- lambda0_*exp(-act.dists[j,i]^2 / (2*hazt.vec.tot^2))
        }
      }
    } else {
      if (haz.type_=='IndL'){ 
        hd <- log(1/(1-g0_*exp(-act.dists^2 / (2*sigma_^2)))) / hazt.cycle.occ1.int
      } else {
        hd <-  lambda0_*exp(-act.dists^2 / (2*sigma_^2))
      }
    }
    
    if (detector.type_=="Prox"){ #independence means can sim per ind and trap
      cap.times <- matrix(list(), nrow=n.inds, ncol=n.traps)
      colnames(cap.times)<- rownames(traps.obj_)
      captures <- matrix(NA, nrow=n.inds, ncol=n.traps)  
      
      #still a NH PP and so ft can be factorised into P(omega)*f(t | Omega)
      #independence mneans can simulate per ind and trap
      #CDF needs to end at 1 and always produce a time
      for (ind in 1:n.inds){ 
        for (det in 1:n.traps) {
          if (haz.type_=='Dep') {   #lambda depends on time and so must be integrated
            lambda <- sum(full.haz[det,,ind])*total.mesh.width
            CDF <- CDF.gam.cyc(ind.haz_ = full.haz[det,,ind], study.vec_ = study.time, hazt.vec.tot_ = NULL, detector.type_ = "Prox", haz.type_ = "Dep")
          } else {
            if (Const.Haz == T){
              lambda <- hd[det, ind]*(endpoint_ - setup.hour_)
            } else {
              lambda <- hd[det, ind] * hazt.total.int
              CDF <- CDF.gam.cyc(ind.haz_ = hd[det,ind], study.vec_ = study.time, hazt.vec.tot_ = hazt.vec.tot, detector.type_ = "Prox", haz.type_ = "Ind")
            }
          }
          
          captures[ind,det] <- rpois(1, lambda) 
          time.vec <- NULL
          
          if (captures[ind,det]>0){
            if (Const.Haz == T){
              time.vec <- c(time.vec, round(runif(captures[ind,det], min=0, max=endpoint_),4))  
            } else {
              for (det.tim in 1:captures[ind,det]) {
                times.all.list <- sim.time.gam(CDF.tot_ = CDF, study.vec_ = study.time, occ.interval_ = occ.interval_,  haz.interval_ = haz.interval, release.hour_ = setup.hour_)
                cap.time <- times.all.list[["Time"]]
                hour.time <- times.all.list[["Haz cycle time"]]
                
                if (is.null(time.vec)){
                  time.vec <- round(cap.time,2)
                } else {
                  time.vec <- c(time.vec,round(cap.time,2))
                }
              }
            }
            
            cap.times[ind,det] <- list(sort(unique(time.vec))) 
            
            times.ind <- cap.times[ind,det][[1]]  
            detect.ind <- rep(rownames(traps.obj_)[det],length(cap.times[ind,det][[1]]))    
            ind.vec <- rep(ind,length(cap.times[ind,det][[1]]))  
            
            Caught <- c(Caught, times.ind)    
            Detector <- c(Detector, detect.ind)  
            id  <- c(id,ind.vec) 
          }
        }
      }
    } else {  #this loop for single catch traps, need to move forward with time
      
      detect.num <- NULL
      times.release <- NULL
      release.times.vec <- NULL
      trap.sat.vec <- rep(0, n.occs)
      
      occ.cycle.num <- 1
      
      if (haz.type_=='Dep') { #1st arg is a vector of hazt per ind
        Tot.haz.t <- apply(full.haz, c(2,3), sum) #summed across traps, times * inds matrix
        Tot.haz.d <- apply(full.haz, c(1,3), sum)*total.mesh.width #integrated over time, traps * inds matrix (like hd)
        CDF.all <- apply(Tot.haz.t,2, CDF.gam.cyc, study.vec_ = study.time, hazt.vec.tot_ = hazt.vec.tot_, detector.type_ = "SC", haz.type_="Dep")
        tot.haz <- colSums(Tot.haz.d) #total hazard exp per individual, at start of loop all traps operating
        prob.mat <- t(t(Tot.haz.d)/tot.haz) #gives traps * n  matrix of trap specific probs
      } else { #1st arg is ind specific vector of hd values
        if (Const.Haz == T){
          CDF.all <- apply(hd,2, function(x) 1 - exp(-sum(x)*(study.time - setup.hour_)))
        } else {
          CDF.all <- apply(hd,2, CDF.gam.cyc, study.vec_ = study.time, hazt.vec.tot_ = hazt.vec.tot, detector.type_ = "SC", haz.type_="Ind")
        }
        tot.hd <- colSums(hd) #total hd exp per individual, at start of loop all traps operating
        prob.mat <- t(t(hd)/tot.hd) #gives traps * n  matrix of trap specific probs
      } 
      
      repeat{
        if (is.null(Caught)){ #this for 1st capture only
          times.all.list <- apply(CDF.all,2, sim.time.gam, study.vec_ = study.time, occ.interval_ = occ.interval_, haz.interval_ = haz.interval, release.hour_ = setup.hour_)
        } else {  #this at start of each new occasion
          occ.cycle.num = occ.cycle.num + 1 
          if (occ.cycle.num > n.occs) {break}  
          
          if (haz.type_=="Dep"){
            full.haz[,study.time <= release.time,] <- 0
            Tot.haz.t.red <- apply(full.haz, c(2,3), sum) #summed across traps, times * inds matrix
            CDF.all <- apply(Tot.haz.t.red,2, CDF.gam.cyc, study.vec_ = study.time, hazt.vec.tot_ = NULL, detector.type_ = "SC", haz.type_="Dep")
            tot.haz <- colSums(Tot.haz.d) #total hazard exp per individual, at start of loop all traps operating
            prob.mat <- t(t(Tot.haz.d)/tot.haz) #gives traps * n  matrix of trap specific probs
          } else {
            if (Const.Haz == T){
              time.left <- study.time - release.time
              time.left[time.left < 0] <- 0
              CDF.all <- apply(hd,2, function(x) 1 - exp(-sum(x)*time.left))
            } else {
              hazt.vec.tot[study.time <= release.time] <- 0
              CDF.all <- apply(hd,2, CDF.gam.cyc, study.vec_ = study.time, hazt.vec.tot_ = hazt.vec.tot, detector.type_ = "SC", haz.type_="Ind")
            }
          }
          times.all.list <- apply(CDF.all,2, sim.time.gam, study.vec_ = study.time, occ.interval_ = occ.interval_, haz.interval_ = haz.interval, release.hour_ = setup.hour_)
        }
        
        times.sim <- sapply(times.all.list,function(x){x[["Time"]]})
        if (all(is.na(times.sim))) {break} 
        
        #this checks if more than 1 min time simulated
        #if so randomly keeps one and sets others to NA
        min.times <- which(times.sim==min(times.sim, na.rm = T))
        num.times <- length(min.times)
        if (num.times > 1) {
          time.incl <- sample(min.times,1)
          times.sim[min.times[which(min.times != time.incl)]] <- NA
        }
        
        #cycles.sim and .num refers to occ cycles
        cycles.sim <- sapply(times.all.list,function(x){x[["Num occ"]]})
        hours.sim <- sapply(times.all.list,function(x){x[["Haz cycle time"]]})
        cap.time <- min(times.sim, na.rm=T)
        
        occ.cycle.num <- cycles.sim[[which(times.sim==cap.time)]]
        hour.time <- hours.sim[[which(times.sim==cap.time)]]
        
        if (haz.type_ == "Dep"){
          ind.caught <- as.numeric(dimnames(full.haz)[[3]][which(times.sim==cap.time)])
        } else {
          ind.caught <- as.numeric(colnames(hd)[which(times.sim==cap.time)])
        } 
        
        
        detect.k <- rmultinom(1,1,prob.mat[,which(times.sim==cap.time)])
        detect <- rownames(detect.k)[which(detect.k==1)]
        Detector <- c(Detector,detect)
        
        release.time <- setup.hour_ + (occ.interval_*occ.cycle.num)
        times.release <- c(times.release, release.time)
        Caught <- c(Caught, cap.time)
        id <- c(id, ind.caught)
        trap.sat <- 1/n.traps
        trap.sat.vec[occ.cycle.num] <- trap.sat
        
        #now need a loop to check for other captures in that "occasion"
        n.detect.occ <- 1
        occ.detect.vec <- tail(Detector, n=1)
        occ.release <- tail(times.release, n=1)
        inds.caught.occ <- ind.caught
        detect.vec.num <- NULL
        #hour.time.occ <- hour.time
        
        repeat{
          if (tail(Caught,1)==occ.release){break}
          if (n.detect.occ==n.inds | n.detect.occ==n.traps) {
            break
          } else {
            if (haz.type_=="Dep"){
              detect.num <- which(dimnames(full.haz)[[1]]==tail(Detector, n=1))
              detect.vec.num <- c(detect.vec.num,detect.num)
              full.haz[,study.time <= tail(Caught,1),] <- 0
              full.haz.red <- full.haz[-detect.vec.num,,-inds.caught.occ]
              
              #deals with one trap or individual l left
              if (is.na(dim(full.haz.red)[3])){
                if (n.inds-length(inds.caught.occ)==1){  #one ind left, full haz is a traps * times
                  CDF.all <- CDF.gam.cyc(ind.haz_ = colSums(full.haz.red), study.vec_ = study.time, hazt.vec.tot_ = NULL, detector.type_ = "SC", haz.type_="Dep")
                  haz.k.ind <- apply(full.haz.red,1,sum)*total.mesh.width
                  tot.haz.red <- sum(haz.k.ind)
                  prob.mat.red <- t(t(haz.k.ind)/tot.haz.red)
                } else {       #one trap left, full haz is a times * inds matrix
                  CDF.all <- apply(full.haz.red,2, CDF.gam.cyc, study.vec_ = study.time, hazt.vec.tot_ = NULL, detector.type_ = "SC", haz.type_="Dep")
                  tot.haz.red <- apply(full.haz.red,2,sum)*total.mesh.width 
                  prob.mat.red <- tot.haz.red/tot.haz.red
                } 
              } else {
                Tot.haz.t.red <- apply(full.haz.red, c(2,3), sum)
                CDF.all <- apply(Tot.haz.t.red,2, CDF.gam.cyc, study.vec_ = study.time, hazt.vec.tot_ = NULL, detector.type_ = "SC", haz.type_="Dep")
                Tot.haz.d.red <- apply(full.haz.red, c(1,3), sum)*total.mesh.width
                tot.haz.red <- colSums(Tot.haz.d.red) 
                prob.mat.red <- t(t(Tot.haz.d.red)/tot.haz.red) 
              }
            } else {
              detect.num <- which(rownames(hd)==tail(Detector, n=1))
              detect.vec.num <- c(detect.vec.num,detect.num)
              hd.red <- as.matrix(hd[-detect.vec.num,-inds.caught.occ])
              
              if (n.inds-length(inds.caught.occ)>1){
                if (dim(hd.red)[2] == 1 ) {hd.red <- t(hd.red)} #transposes h mat if only 1 trap left
              }
              
              if (Const.Haz == T){
                time.left <- study.time - tail(Caught,1)
                time.left[time.left < 0] <- 0 
                CDF.all <- apply(hd.red,2, function(x) 1 - exp(-sum(x)*time.left))
              } else {
                hazt.vec.tot[study.time <= tail(Caught,1)] <- 0
                CDF.all <- apply(hd.red,2, CDF.gam.cyc, study.vec_ = study.time, hazt.vec.tot_ = hazt.vec.tot, detector.type_ = "SC", haz.type_="Ind")
              }
              
              tot.hd.red <- colSums(hd.red)
              prob.mat.red <- t(t(hd.red)/tot.hd.red)
            }
            
            if (n.inds-length(inds.caught.occ)==1){
              times.all.occ <- sim.time.gam(CDF.tot_ = CDF.all, study.vec_ = study.time, occ.interval_ = occ.interval_, haz.interval_ = haz.interval, release.hour_ = setup.hour_)
              times.occ <- times.all.occ[["Time"]]
            } else {
              times.all.occ <- apply(CDF.all, 2, sim.time.gam, study.vec_ = study.time, occ.interval_ = occ.interval_, haz.interval_ = haz.interval, release.hour_ = setup.hour_)
              times.occ <- sapply(times.all.occ,function(x){x[["Time"]]})
            }
            
            if (all(is.na(times.occ))) {break} 
            
            min.times.occ <- which(times.occ==min(times.occ, na.rm = T))
            num.times.occ <- length(min.times.occ)
            if (num.times.occ > 1) {
              time.incl.occ <- sample(min.times.occ,1)
              times.occ[min.times.occ[which(min.times.occ != time.incl.occ)]] <- NA
            }
            
            next.time <- min(times.occ, na.rm=T)
            if (next.time > occ.release){break} 
            
            Caught <- c(Caught, next.time)
            
            #the if statement handles the case where only 1 ind or trap left
            if (n.inds-length(inds.caught.occ)==1){
              if (haz.type_=="Dep"){
                new.ind.caught <- as.numeric(setdiff(dimnames(full.haz)[[3]], inds.caught.occ))
              } else {
                new.ind.caught <- as.numeric(setdiff(colnames(hd), inds.caught.occ))
              }
              detect.k <- rmultinom(1,1,prob.mat.red)
              detect <- rownames(prob.mat.red)[which(detect.k==1)]
            } else {
              if (haz.type_=="Dep"){
                if (n.traps - n.detect.occ == 1) {
                  new.ind.caught <- as.numeric(names(prob.mat.red)[which(times.occ==next.time)])
                  detect <- setdiff(dimnames(full.haz)[[1]], occ.detect.vec) 
                } else {
                  new.ind.caught <- as.numeric(dimnames(full.haz.red)[[3]][which(times.occ==next.time)])
                  detect.k <- rmultinom(1,1,prob.mat.red[,colnames(prob.mat.red)==new.ind.caught])
                  detect <- rownames(prob.mat.red)[which(detect.k==1)]
                }
              } else {
                new.ind.caught <- as.numeric(names(tot.hd.red[which(times.occ==next.time)]))
                if (dim(hd.red)[1]==1){
                  detect <- setdiff(rownames(hd), occ.detect.vec) 
                } else {
                  detect.k <- rmultinom(1,1,prob.mat.red[,colnames(prob.mat.red)==new.ind.caught])
                  detect <- rownames(prob.mat.red)[which(detect.k==1)]
                }
              }
            }
            
            inds.caught.occ <- c(inds.caught.occ, new.ind.caught)
            n.detect.occ <- n.detect.occ+1
            occ.detect.vec <- c(occ.detect.vec, detect)
            trap.sat <- length(occ.detect.vec)/n.traps
            trap.sat.vec[occ.cycle.num] <- trap.sat
            Detector <- c(Detector,detect)
            id <- c(id,new.ind.caught)
            times.release <- c(times.release, occ.release)
          } #closes else statement for any inds left
        }#closes inner repeat loop (checks for other caps in 'occasion')
      }#closes outer loop
    }#closes loop for SC
    
    if (is.null(Caught)){return('No captures')
    } else {
      session <- rep(1,length(Caught))
      
      if (detector.type_=="Prox"){
        data <- data.frame(session, id, Caught, Detector) 
      } else { 
        data <- data.frame(session, id, Caught,'Release'=times.release, Detector) 
      }
      
      data <- data[order(data$id,data$Caught),]
      data$id <- reorder(data$id)
      data <- remove.duplicates(data.frame_ = data, n.inds_ = length(unique(data$id)))
      rownames(data) <- seq(1,dim(data)[1])
      
      #create capthist object from cts CH  
      breaks <- seq(setup.hour_, endpoint_, occ.interval_)
      times.cat <- findInterval(data$Caught, breaks, rightmost.closed=T)
      capthist <- data.frame(data$session, data$id, times.cat, data$Detector)
      CH.prox <- make.capthist(capthist, traps.obj_, noccasions = n.occs)
      CH <- secr::reduce(CH.prox, dropunused=F)    #reduce counts > 1 to equal 1
      return(list('cts' = data, 'CH prox' = CH.prox, 'CH binary' = CH))
    } 
  } else {return("No animals")}
}

#' @title Remove duplicate times.
#' @description Gets rid of any duplicate times that may have been simulated, for the same individual at the same detector. 
#' This function is called by sim.ctsCH. Note that it requires variables called "Caught" and "id" in the data frame.
#'    
#' @return Returns the capture histories without the duplicate times..  
#' 
#' @param data.frame_ The capture history data.
#' @param n.inds_ The number of unique individuals in the capture history data.
#' 
#' @export
remove.duplicates <- function(data.frame_, n.inds_){
  new.data <- NULL
  id.col <- which(names(data.frame_)=="id")
  id.names <- unique(data.frame_[, id.col])
  
  for (i in 1:n.inds_){
    data.ind <- data.frame_[data.frame_$id == id.names[i],]
    if (dim(data.ind)[1]>1){
      new.data.ind <- data.ind[match(unique(data.ind$Caught), data.ind$Caught),]
      new.data <- rbind(new.data, new.data.ind)
    } else {new.data <- rbind(new.data, data.ind)}
  }
  return(new.data)
}

#' @title Reorders a vector of id's.
#' @description This function is called by sim.ctsCH and reorders the id numbers in the simulated data so that the numbers run
#' from 1 to n.
#'    
#' @return Returns a vector of reordered ids.  
#' 
#' @param id.vec A vector of individual ids.
#' 
#' @export
reorder <- function(id.vec){
  j=1
  newid=NULL
  for (k in 1:length(unique(id.vec))){
    match <- unique(id.vec)[k]
    reorder=rep(j,length(id.vec[id.vec==match]))
    newid=c(newid,reorder)
    j=j+1
  }
  newid	
}

#' @title Construct the CDF for capture times.
#' @description Constructs the Cumulative Distribution Function (CDF) for capture times conditional on a frequency or 
#' count of captures. The conditional nature opf the CDF means it must end at 1 and hence a time will always be generated.
#' Called by the simctsCH function.
#'    
#' @return Returns a matrix of dimension number traps * time steps that represents F(t) or the prob of being caught before time t.
#' 
#' @param ind.haz_ The individual and detector-specific detection hazard(s). For a dependent hazard, the argument will be a vector 
#' of the actual encounter rates whereas for an independent hazard it will be a scalar representing lambda_d. For SC traps the
#' argument will be a vector of values across trpas that gets summed for an overall hazard.
#' @param study.vec_ A vector of times that correspond to the time period one wants to simulate captures for.
#' @param hazt.vec.tot_ The h(t) component of the encounter rate function i.e. the detection hazard at different cycle times.
#' Only relevant for independent hazards. 
#' @param detector.type_ Specifies the type of detector. Defaults to using a "proximity" or passive detector but can also simulate
#' from an array of single-catch traps. 
#' @param haz.type_ The type of detection hazard to use in the simulation. The choice is between an independent hazard ("Ind) 
#' or a dependent hazard ("Dep).
#' 
#' @export
CDF.gam.cyc <- function(ind.haz_, study.vec_, hazt.vec.tot_, detector.type_ = "Prox", haz.type_ = "Ind"){
  
  ft <- vector(mode='numeric', length = length(study.vec_))
  Ft <- vector(mode='numeric', length = length(study.vec_))
  
  mesh.width <- study.vec_[2]-study.vec_[1]
  
  if (detector.type_ == "SC"){
    sum.tot <- 0
    if (haz.type_ == "Dep"){
      for (i in 1:length(ind.haz_)){
        sum.tot <- sum.tot + ind.haz_[i]
        Ft[i] <-  1 - exp(-(sum.tot*mesh.width))
      }
    } else {
      tot.hd <- sum(ind.haz_)
      for (i in 1:length(study.vec_)){
        if (is.null(hazt.vec.tot_)){  #indicates constant hazard
          sum.tot <- sum.tot + tot.hd
        } else {
          sum.tot <- sum.tot + tot.hd*hazt.vec.tot_[i]
        }
        Ft[i] <-  1 - exp(-(sum.tot*mesh.width))
      }
    }
  } else {  #proximity detectors
    if (haz.type_ == "Dep"){
      haz.tot.int <- sum(ind.haz_)*(mesh.width)
      ft <- ind.haz_ / haz.tot.int
    } else {
      haz.tot.int <- ind.haz_*sum(hazt.vec.tot_)*(mesh.width)
      ft <- (ind.haz_*hazt.vec.tot_) / haz.tot.int
    }
    Ft <- cumsum(ft)*mesh.width
  }
  return("CDF" = Ft)
}

#' @title Simulate capture time.
#' @description Called by the simctsCH function to simulate individual capture times using the CDF of capture times.
#'    
#' @return Returns both the actual capture time in total time units, and in terms of the hazard cycle. Also returns the 
#' capture occasion.
#' 
#' @param CDF.tot_ The CDF that is used to simulate capture times via the probability integral transform approach.
#' @param study.vec_ A vector of times that correspond to the time period one wants to simulate captures for.
#' @param occ.interval_ The duration of an occasion.
#' @param haz.interval_ The duration of a single hazard cycle.
#' @param release.hour_ Relevant for SC traps when individuals are released at a particular time. Usually in such cases the 
#' release time will define the occasions.
#' 
#' @export
sim.time.gam <- function(CDF.tot_, study.vec_, occ.interval_ = 24, haz.interval_=24, release.hour_ = 0){
  
  time.pos <- NULL
  time <- NULL
  occ.num <- NULL
  haz.cycle.time <- NULL
  last.zero <- 0
  
  if (min(CDF.tot_)==0) {last.zero <- tail(which(CDF.tot_==0),1)}
  
  p <- runif(1)
  
  if (max(CDF.tot_) > p){
    if (p < CDF.tot_[last.zero+1] / 2) {time.pos <- last.zero + 1 } else {
      time.pos <- which.min(abs(CDF.tot_ - p))
    }
    
    time <- study.vec_[time.pos] 
    
    occ.num <- ceiling((time - release.hour_)/occ.interval_)
    
    if(occ.num==0){occ.num <- 1}
    
    haz.cycle.time <- time %% haz.interval_
    
  } else {time <- NA} 
  return(list("Num occ" = occ.num, "Haz cycle time" = haz.cycle.time, "Time" = time))
}

######################################################################
#fns related to CT LL
######################################################################

#' @title Constructs the Continuous-Time (CT) SCR likelihood.
#' @description Constructs the CT likelihood that is required to fit these models using maximum likelihood and a half normal form 
#' for the encounter rate / detection function.
#'    
#' @return Returns the negative log likelihood value. The chosen optimiser used to find parameter estimates is actually a minimiser
#' and minimising the negative log likelihood is equivalent to maximising the log likelihood.  
#' 
#' @param beta_ A vector of parameter values. These contain the starting values for the optimisation when passed to nlm.
#' @param data_ The data frame. The capture times variable needs to be called "Caught", the detector column "Detector", and 
#' the identifier column "id".
#' @param trap.objs_ The secr trap object. This argument needs to be a list when a model is fitted to more than one session. 
#' @param mask.objs_ The secr mask object. This argument needs to be a list when a model is fitted to more than one session. 
#' @param dists_ A matrix of distances between the detectors and the mask points. If this is NULL it will be created by the function.
#' @param D.type_ A string vector that specifies the density model to be fitted for each session. Options include constant, exponential 
#' and quadratic density ( c("Cons","Exp", "Quad"). The length of the vector must equal the number of sessions in the data. 
#' @param D.mod_ A numeric vector that specifies whether the same density model and parameter(s) is used for the different sessions.
#' A different number specifies that different parameters are estimated. The length of the vector must equal the number of sessions in 
#' the data. 
#' @param lambda0.mod_ A numeric vector that specifies whether the same intercept parameter for the encounter rate / detection function
#'  is used for the different sessions. A different number specifies that different parameters are estimated. The length of the vector 
#'  must equal the number of sessions in the data.  
#' @param sigma.mod_ A numeric vector that specifies whether the same range parameter for the encounter rate / detection function
#'  is used for the different sessions. A different number specifies that different parameters are estimated. The length of the vector 
#'  must equal the number of sessions in the data. 
#' @param haz.mod_ A numeric vector that specifies whether the same spline parameters for the encounter rate / detection function
#'  is used for the different sessions. A different number specifies that different parameters are estimated. The length of the vector 
#'  must equal the number of sessions in the data.  
#' @param lambda0.b_ Specifies if an overall behavioural effect on the intercept parameter for the encounter rate / detection function
#' is required. Not developed yet for trap-specific effects or for the range parameter.
#' @param haz.type_ The type of detection hazard to use in the model. Choices c("IndL", "IndU", "Dep") include the two 
#' independent hazard parameterisations and the dependent hazard parameterisation.
#' @param haz.k_ A numeric vector that specifies the df for the regression spline for the different sessions. A different number 
#' specifies different df's for different sessions. A constant haz is specified by giving a value of 1,and can be used for either
#' the linked or unlinked independent hazard. The length of the vector must equal the number of sessions in the data.
#' @param conduct.mod_ A numeric vector that specifies whether the same conductance parameter is used for the different sessions. 
#' A different number specifies that different parameters are estimated. The length of the vector must equal the number of sessions in the data.
#' @param fixed.b0_ A scalar for the b0 spline parameter. Although it is relevant for the independent hazards, it is redundant
#' due to the fact that lambda_d scales the estimated encounter rate. 
#' @param usage_ A list that supplies information about trap usage i.e. between what time periods each detector was NOT operating.
#' @param endpoint_ The endpoint of the study that determines the study duration. The length of the vector must equal the number of sessions in the data. 
#' @param setup.hour_ The starting point of the study. Values of zero are assumed unless alternative values are supplied, in which case
#' the length of the vector must equal the number of sessions in the data.
#' Can be used if the survey started after time 0 i.e. if it does not start at the start of a hazard cycle.
#' @param occ.interval_ The duration of an occasion. Primarily used when integrating the hazard function over one occasion.
#' @param cycle.vec_ A vector of times that correspond to a single hazard cycle.
#' @param mesh.size_ Determines the resolution of the time vector. This argument determines the vector length from the seq() command.
#' 
#' @export
CT.SCR <- function (beta_, data_, trap.objs_, mask.objs_, dists_ = NULL, D.type_=c('Cons'), D.mod_ = c(1), lambda0.mod_ = c(1), sigma.mod_ = NULL, lambda0.b_ = c('n'), haz.mod_=c(1), haz.type_ = 'IndL', haz.k_ = c(6), conduct.mod_ = NULL, fixed.b0_ = 0, usage_ = NULL, endpoint_, setup.hour_ = 0, occ.interval_= 24, cycle.vec_ = seq(0,24,by=1), mesh.size_=100, trace_ = FALSE) {
  
  num.sessions <- length(haz.k_)
  
  #calc position of density and spline parameters
  #for dens it is the final position
  
  pos.vecs <- Calc.pos.Haz(D.mod_ = D.mod_, D.type_ = D.type_, lambda0.mod_ = lambda0.mod_, b.effect_ = lambda0.b_, sigma.mod_ = sigma.mod_, haz.type_ = haz.type_, haz.mod_ = haz.mod_, haz.k_ = haz.k_, conduct.mod_ = conduct.mod_)

  dens.pos.vec <- pos.vecs[["Density"]]
  spline.pos.vec <- pos.vecs[["Spline"]]
  
  if (haz.type_=='Dep') {
    lambda0.pos.vec <- pos.vecs[["Lambda0"]]
  } else {
    if (haz.type_=='IndL') {
      sigma.pos.vec <- pos.vecs[["Sigma"]] 
      g0.pos.vec <- pos.vecs[["G0"]]
      beta.1 <- fixed.b0_
    } else {
      sigma.pos.vec <- pos.vecs[["Sigma"]] 
      lambda0.pos.vec <- pos.vecs[["Lambda0"]]
      beta.1 <- fixed.b0_
    }
  }
  
  LL <- 0
  
  for (i in 1:num.sessions){
    Const.Haz = NULL
    
    if (haz.k_[i] > 1 & haz.k_[i] < 4){print("Dof for haz must be a min of 4") ; break}
    
    if (haz.k_[i] == 1){Const.Haz = T} else {Const.Haz = F}
    
    if (num.sessions==1){
      trap.obj <- trap.objs_
      mask.obj <- mask.objs_
      usage.info <- usage_
      endpoint <- endpoint_
      setup.hour <- setup.hour_
    } else {
      trap.obj <- trap.objs_[[i]]
      mask.obj <- mask.objs_[[i]]
      usage.info <- usage_[[i]]
      endpoint <- endpoint_[[i]]
      if (length(setup.hour_)==1){
        setup.hour <- 0
      } else {
        setup.hour <- setup.hour_[[i]]
      }
    }
    
    trap.names <- rownames(trap.obj)
    
    #get counts of ids
    n.s <- length(unique(data_[data_$session==i,which(names(data_)=="id")]))
    n.traps <- dim(trap.obj)[1]
    
    #get appropriate dist matrix
    
    if (is.null(conduct.mod_)){
      if (is.null(dists_)) {
        dist.s <- scrmlebook:::distances(trap.obj, mask.obj)
      } else {  #if distance matrix is provided in a list
        dist.s <- dists_[[i]]
      }
    } else {
      conduct.pos.vec <- pos.vecs[["Conductance"]]
      conduct <- beta_[conduct.pos.vec[i]]
      dist.s <- noneuc.distances(traps.obj_ = trap.obj , mask.obj_ = mask.obj, transitionFn_ = 1, directions_ = 16, conduct.par_ = conduct)
    }
    
    ##standardise X coords for D models
    Zx.s <- (mask.obj$x-mean(mask.obj$x)) / sd(mask.obj$x)
    Zx2.s <- Zx.s^2
    A <- attr(mask.obj, 'area')
    
    ##Set spline hazard and detection fn parameters
    
    if (haz.type_=='Dep') {
      lambda0 <- exp(beta_[lambda0.pos.vec[i]])
      if (lambda0.b_[i]=='y'){
        lambda0.b <- exp(beta_[lambda0.pos.vec[i]+1])
      } 
      hazt.beta.vec.s <- c(beta_[spline.pos.vec[i]:(spline.pos.vec[i] + haz.k_[i] - 2)]) #one less if hazard != Dep
    } else {
      if (haz.type_=='IndL') {
        sigma <- exp(beta_[sigma.pos.vec[i]])
        g0 <- invlogit(beta_[g0.pos.vec[i]])
        if (lambda0.b_[i]=='y'){
          g0.b <- invlogit(beta_[g0.pos.vec[i]+1])
        }
        beta.1 <- fixed.b0_
        hazt.beta.vec.s <- c(beta.1, beta_[spline.pos.vec[i]:(spline.pos.vec[i] + haz.k_[i] - 3)])
      } else {
        sigma <- exp(beta_[sigma.pos.vec[i]])
        lambda0 <- exp(beta_[lambda0.pos.vec[i]])
        if (lambda0.b_[i]=='y'){
          lambda0.b <- exp(beta_[lambda0.pos.vec[i]+1])
        } 
        beta.1 <- fixed.b0_
        hazt.beta.vec.s <- c(beta.1, beta_[spline.pos.vec[i]:(spline.pos.vec[i] + haz.k_[i] - 3)])
      }
    }
    
    #Now density parameters
    if (D.type_[i]=='Exp'){
      D1 <- beta_[dens.pos.vec[i]-1]
      D2 <- beta_[dens.pos.vec[i]]
      Dvec <- exp(D1 + D2*Zx.s) #density vector
      Dvec[Dvec==0]=2.225074e-308 
      Dvec[is.infinite(Dvec) ] = 1.797693e+308
      log.Dvec <- log(Dvec)
      log.Dvec[is.infinite(log.Dvec)]=(-1.797693e+308)/50
    } else {
      if (D.type_[i]=='Quad'){
        D1 <- beta_[dens.pos.vec[i]-2]
        D2 <- beta_[dens.pos.vec[i]-1]
        D3 <- beta_[dens.pos.vec[i]]
        Dvec <- exp(D1 + D2*Zx.s + D3*Zx2.s) #density vector
        Dvec[Dvec==0]=2.225074e-308 
        Dvec[is.infinite(Dvec) ] = 1.797693e+308
        log.Dvec <- log(Dvec)
        log.Dvec[is.infinite(log.Dvec)]=(-1.797693e+308)/50
      } else {
        D <- exp(beta_[dens.pos.vec[i]])
      }
    }
    
    cycle.interval = max(cycle.vec_) - min(cycle.vec_)
    cycle.mesh <- seq(min(cycle.vec_), max(cycle.vec_), length=mesh.size_)
    cycle.mesh.width = cycle.mesh[2] - cycle.mesh[1]
    
    if (setup.hour %% cycle.interval != 0){
      n.haz.cycles <- ((endpoint - setup.hour)-(cycle.interval-setup.hour %% cycle.interval)-(endpoint %% cycle.interval)) / cycle.interval
    } else {
      n.haz.cycles <- ((endpoint - setup.hour)-(endpoint %% cycle.interval)) / cycle.interval
    } #note that n.haz.cycles is the number of complete haz cycles
    
    if (Const.Haz==T){
      study.time <- seq(setup.hour,endpoint, length=length((cycle.mesh-1)*n.haz.cycles))
      cycle.int.occ1 <- occ.interval_
      full.cycle.int.s <- endpoint
      smoothsetup.s = NULL
    } else {
      #create time dependent spline hazard as gam object
      smoothsetup.s <- create.GAM(cycle.vec_, rep(0,length(cycle.vec_)),haz.k_[i])[[1]]
      hazt.cycle.s <- pred.Y(time.vec_ = cycle.mesh, gam.obj_ = smoothsetup.s, spline.coeff_ = hazt.beta.vec.s)
      if (any(is.infinite(hazt.cycle.s))) {LL <- -1e10 ; break}
      
      #now generate the time-dependent sigma vector over full study
      hazt.study.s <- study.haz.vec(cycle.vec_ = cycle.mesh, y.vec.full_ = hazt.cycle.s, spline.coeff_ = hazt.beta.vec.s, release.hour_ = setup.hour, endpoint_ = endpoint, n.haz.cycles_ = n.haz.cycles, haz.interval_ = cycle.interval, gam.obj_ = smoothsetup.s)
      study.time <- seq(setup.hour,endpoint, length=length(hazt.study.s))
      tot.mesh.width <- (study.time[2] - study.time[1])
      
      cycle.int.occ1 <- sum(hazt.study.s[study.time >= setup.hour & study.time <= (setup.hour + occ.interval_)])*tot.mesh.width
      full.cycle.int.s <- sum(hazt.study.s)*tot.mesh.width
    }
    
    #Calc appropriate haz, depends on haz.type
    #all scenarios end up with the cumulative haz vector H.s (1 * K)
    if (haz.type_=='Dep') {
      #calculate full hazard matrix, depends on both distance and time
      lambda.dt.s = array(NA, dim =c(n.traps,length(study.time),dim(mask.obj)[1]))
      dimnames(lambda.dt.s) <- list(trap.names,round(study.time,2),NULL)
      for (k in 1:dim(mask.obj)[1]){
        for (j in 1:n.traps){
          lambda.dt.s[j,,k] <- lambda0*exp(-dist.s[j,k]^2 / (2*hazt.study.s^2))
        }
      }
      lambda.dt.s[is.infinite(lambda.dt.s) ] = 1.797693e+308
      H.k.s <- apply(lambda.dt.s,c(1,3),function(x) sum(x)*tot.mesh.width)
      rm(lambda.dt.s)  #removes array as not used again, for memory purposes
      gc()
      H.s <- colSums(H.k.s)
    } else {
      if (haz.type_=='IndL') {  #h0 and time seperated
        lambda.d.s <- (log(1/(1-g0*exp(-dist.s^2 / (2*sigma^2))))) / cycle.int.occ1
        lambda.d.s[is.infinite(lambda.d.s) ] = 1.797693e+308
        log.lambda.d.s <- log(lambda.d.s)
        log.lambda.d.s[is.infinite(log.lambda.d.s)]=(-1.797693e+308)/50
        
        if (lambda0.b_[i]=='y'){
          lambdab.d.s <- (log(1/(1-g0.b*exp(-dist.s^2 / (2*sigma^2))))) / cycle.int.occ1
          lambdab.d.s[is.infinite(lambdab.d.s) ] = 1.797693e+308
          log.lambdab.d.s <- log(lambdab.d.s)
          log.lambdab.d.s[is.infinite(log.lambdab.d.s)]=(-1.797693e+308)/50
        }
      } else {
        lambda.d.s <- lambda0*exp(-dist.s^2 / (2*sigma^2))
        lambda.d.s[is.infinite(lambda.d.s) ] = 1.797693e+308
        log.lambda.d.s <- log(lambda.d.s)
        log.lambda.d.s[is.infinite(log.lambda.d.s)]=(-1.797693e+308)/50
        
        if (lambda0.b_[i]=='y'){
          lambdab.d.s <- lambda0.b*exp(-dist.s^2 / (2*sigma^2))
          lambdab.d.s[is.infinite(lambdab.d.s) ] = 1.797693e+308
          log.lambdab.d.s <- log(lambdab.d.s)
          log.lambdab.d.s[is.infinite(log.lambdab.d.s)]=(-1.797693e+308)/50
        }
      }
      H.k.s <- full.cycle.int.s * lambda.d.s 
      H.k.s[is.infinite(H.k.s) ] = 1.797693e+308
      H.s <- colSums(H.k.s)  #gives cumulative/integrated hazard for all traps over study, for uncaught 
    }    
    
    #create usage intervals for b effect
    if (lambda0.b_[i]=='y'){
      if (is.null(usage.info)){
        outage.intervals=NULL
      } else {
        outage.intervals <- t(sapply(trap.names, intervals.outage, usage.info_=usage.info))
        names(outage.intervals) <- trap.names  
      }
    }
    
    if (is.null(usage.info)) {
      Surv.total.s <- exp(-H.s)     
      #so this gives the prob of no capture, at any trap, given X
      #uses 1st lambda0 if Mb, as does calculation of exposure below
      p.m.s <- 1 - Surv.total.s
    } else {
      Trap.exposures.s <- matrix(NA, nrow = dim(trap.obj)[1], ncol=dim(mask.obj)[1])
      if (haz.type_=='Dep'){
        Trap.exposures.s <- t(sapply(trap.names, LL.Exposure, haz.type_ = haz.type_, lambda.d_ = NULL, full.exp_= H.k.s, dists_ = dist.s, traps.fail.file_ = usage.info, lambda0_ = lambda0, const.haz_= Const.Haz, study.time_ = study.time, gam.obj_ = smoothsetup.s, spline.coeff_ = hazt.beta.vec.s))
      } else {
        Trap.exposures.s <- t(sapply(trap.names, LL.Exposure, haz.type_ = haz.type_, lambda.d_ = lambda.d.s, full.exp_= H.k.s, dists_ = NULL, traps.fail.file_ = usage.info, lambda0_ = NULL, const.haz_= Const.Haz, study.time_ = study.time, gam.obj_ = smoothsetup.s, spline.coeff_ = hazt.beta.vec.s))
      }
      Surv.total.s <- colSums(Trap.exposures.s)
      p.m.s <- 1 - exp(-Surv.total.s)
    }
    
    if (D.type_[i] == 'Cons'){
      a.s <- sum(p.m.s)*A      #effective survey area
      if (a.s==0) a.s=2.225074e-308  #replaces a with a tiny number if =0 
      Da.s <- D*a.s
    } else {
      Da.s <- sum(Dvec*p.m.s)*A   #filtered density (after integrating out distance)
    }
    
    #compute likelihood term for each capture time and for each trap
    #depends on LL.cap function
    LL.ind.s <- matrix(NA, nrow = n.s, ncol=dim(mask.obj)[1])
    sum.t.s <- NULL
    int.x.s <- NULL
    
    #now log likelihood for capture times
    data.s <- data_[data_$session==i,]
    id.col <- which(names(data.s)=="id")
    id.names <- unique(data.s[, id.col])
    
    for (j in 1:n.s){
      ind.data <- data.s[data.s$id==id.names[j],]
      if (lambda0.b_[i]=='y'){
        b <- rep(0, dim(ind.data)[1])
        b[1] <- 1  
        ind.data <- cbind(ind.data,b)
      }
      
      LL.haz.times<-matrix(NA,nrow=dim(ind.data)[1],ncol=dim(mask.obj)[1])
      
      if (haz.type_=='Dep'){
        if (lambda0.b_[i] == 'y'){
          LL.haz.times <- t(apply(ind.data, 1, LL.captimes, dists_ = dist.s, log.lambda.d_ = NULL, lambda0_ = lambda0, lambda0.b_ = lambda0.b, gam.obj_=smoothsetup.s, spline.coeff_ = hazt.beta.vec.s))
          Hk.s.i <- LL.survival.ind(haz.type_ = haz.type_, trap.obj_ = trap.obj, first.cap_ = ind.data[1,cap.col], usage.intervals_ = outage.intervals, lambda.d_ = NULL, lambda.b.d_ = NULL, lambda0_ = lambda0, lambda0.b_ = lambda0.b,  const.haz_= Const.Haz, dists_ = dist.s, gam.obj_ = smoothsetup.s, spline.coeff_ = hazt.beta.vec.s, study.time_ = study.time, setup.hour_ = setup.hour, endpoint_ = endpoint, mask.pts_ = dim(mask.obj)[1])
        } else {
          LL.haz.times <- t(apply(ind.data, 1, LL.captimes, dists_ = dist.s, log.lambda.d_ = NULL, lambda0_ = lambda0, gam.obj_=smoothsetup.s, spline.coeff_ = hazt.beta.vec.s))
        }
      } else {
        if (lambda0.b_[i] == 'y'){
          LL.haz.times <- t(apply(ind.data, 1, LL.captimes, dists_ = NULL, log.lambda.d_ = log.lambda.d.s, log.lambda.b.d_ = log.lambdab.d.s, lambda0_ = NULL, const.haz_ = Const.Haz, gam.obj_=smoothsetup.s, spline.coeff_ = hazt.beta.vec.s))
          Hk.s.i <- LL.survival.ind(haz.type_ = haz.type_, trap.obj_ = trap.obj, first.cap_ = ind.data[1,cap.col], usage.intervals_ = outage.intervals, lambda.d_ = lambda.d.s, lambda.b.d_ = lambdab.d.s, lambda0_ = NULL, lambda0.b_ = NULL,  const.haz_= Const.Haz, dists_ = dist.s, gam.obj_ = smoothsetup.s, spline.coeff_ = hazt.beta.vec.s, study.time_ = study.time, setup.hour_ = setup.hour, endpoint_ = endpoint, mask.pts_ = dim(mask.obj)[1])
        } else {
          LL.haz.times <- t(apply(ind.data, 1, LL.captimes, dists_ = NULL, log.lambda.d_ = log.lambda.d.s, log.lambda.b.d_ = NULL, lambda0_ = NULL, const.haz_ = Const.Haz, gam.obj_=smoothsetup.s, spline.coeff_ = hazt.beta.vec.s))
        }
      }
      
      if (lambda0.b_[i] == 'y'){
        if (D.type_[i] == 'Cons'){
          sum.t.s <- colSums(LL.haz.times) - colSums(Hk.s.i)
        } else {
          sum.t.s <- colSums(LL.haz.times) - colSums(Hk.s.i) + log.Dvec
        }
      } else {
        if (D.type_[i] == 'Cons'){
          sum.t.s <- colSums(LL.haz.times) - H.s
        } else {
          sum.t.s <- colSums(LL.haz.times) - H.s + log.Dvec
        }
      }
      LL.ind.s[j,] <- sum.t.s
    }
    
    int.x.s <- rowSums(exp(LL.ind.s))*A #unlogs and integrates across x, L1 relogs and sums across individuals
    int.x.s[which(int.x.s==0)]=2.225074e-308
    
    L1.s <- sum(log(int.x.s))
    
    if (D.type_[i]=='Exp'|D.type_[i]=='Quad'){
      L2.s <- -Da.s
      LL <- LL + L1.s + L2.s
    } else {
      L2.s <- -n.s*log(a.s)
      L3.s <- dpois (n.s, Da.s, log = TRUE)
      LL <- LL + L1.s + L2.s + L3.s
    }
  }
  
  if (!is.finite(LL))
    1e10
  else
    -LL   ## return negative log likelihood
}

#' @title Calculate the extent of exposure to each detector. 
#' @description Called by the CT.SCR function and used to calculate the total exposure for each trap. Usage information will
#' obviously reduce this exposure for traps not operating for the full survey duration. 
#'    
#' @return Returns a vector of total exposure
#' 
#' @param x_ The function is used with apply() and the trap names and so the first argument x_ is the 
#' trap name that gets passed by apply().
#' @param haz.type_ The type of detection hazard to use in the model. Choices c("IndL", "IndU", "Dep") include the two 
#' independent hazard parameterisations and the dependent hazard parameterisation.
#' @param lambda.d_ If this is Null then it implies that a dependent hazard is used
#' @param full.exp_ The total exposure (integrated hazard) corresponding to a full hazard cycle for all detectors.
#' @param const.haz_ Indicates if a constant hazard is used.
#' @param dists_ A matrix of distances between detectors and mask points.
#' @param traps.fail.file_ A list of trap usage information.
#' @param lambda0_ Value for the intercept of the encounter rate function. Only needed when a dependent hazard is specified.
#' @param study.time_ A vector corresponding to the study duration.
#' @param gam.obj_ GAM obj used to predict the detection hazard for the outage periods.
#' @param spline.coeff_ A vector of coefficients to be used in the regression spline.
#' 
#' @export
LL.Exposure <- function(x_, haz.type_ = "IndL", lambda.d_ = NULL, full.exp_, const.haz_ = F, dists_ = NULL, traps.fail.file_, lambda0_ = NULL, study.time_, gam.obj_, spline.coeff_){
  
  if (haz.type_ !="Dep") {lambda.d.k <- lambda.d_[rownames(lambda.d_)==x_,]}
  
  H.k <- full.exp_[rownames(full.exp_)==x_,]
  
  haz.mesh.width <- study.time_[2]-study.time_[1]
  Survk <- NULL
  excl <- NULL
  
  if (is.null(traps.fail.file_[[x_]])){
    excl <- 0
  } else {
    
    if (haz.type_=="Dep"){dist.k <- dists_[rownames(dists_)==x_,]}
    
    for (i in 1:length(traps.fail.file_[[x_]])){
      mesh.start <- traps.fail.file_[[x_]][[i]][[1]]
      mesh.end <- traps.fail.file_[[x_]][[i]][[2]]
      period.mesh <- study.time_[study.time_ > mesh.start & study.time_ <= mesh.end]
      
      if (const.haz_ == T){
        period = (mesh.end - mesh.start) * lambda.d.k
      } else {
        if (haz.type_ !="Dep"){
          period <- sum(pred.Y(time.vec_ = period.mesh, gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)) * haz.mesh.width * lambda.d.k
        } else {
          #find vector of sigma values that corrrespond to outage time
          sig.vec <- pred.Y(time.vec_ = period.mesh, gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)
          
          #calc total hazard values for outage
          #for a specific trap but has distance * time
          lambda.outage <- matrix(NA, nrow = length(H.k), ncol = length(period.mesh))
          
          for (d in 1:length(dist.k)){
            lambda.outage[d,] <- lambda0_*exp(-dist.k[d]^2 / (2*sig.vec^2))
          }
          
          #integrate across time
          period <- apply(lambda.outage,1,sum)*haz.mesh.width
        }
      }
      
      if (is.null(excl)){
        excl <- period            
      } else {
        excl <- excl+period
      }
    }
  }
  Survk <- H.k-excl
  return(Survk)
}

#' @title Calculate the log likelihood value for captures. 
#' @description calculates the contribution to the log likelihood term for each capture history. Called by the 
#' CT.SCR function and used with the apply() function. 
#'    
#' @return Returns a vector for each capture time of dimension equal to the number of mask points.
#' 
#' @param x_ The function is used with apply() on each individual's capture history. The first argument x_ is a  
#' row of data from the capture history data frame.
#' @param dists_ A matrix of distances between detectors and mask points. Used when a dependent hazard is specified.
#' @param log.lambda.d_ A matrix with the log encounter rates for each detector corresponding to each mask point, i.e. the rates 
#' exclude the time component. Only supplied for independent hazards, if it is NULL a dependent hazard is implied.
#' @param log.lambda.b.d_ Another log lambda matrix corresponding to distances but using a different parameter estimate after 
#' first capture. Only relevant for models that have a behavioural effect specified.
#' @param lambda0_ The intercept parameter for the dependent hazard encounter rate function.
#' @param lambda0.b_ A second intercept parameter for the dependent hazard encounter rate function that is applicable after first
#' capture. Only relevant for models that have a behavioural effect specified.
#' @param const.haz_ Indicates a hazard that does not cvary through time.
#' @param gam.obj_ GAM obj used to find the detection hazard for the capture times.
#' @param spline.coeff_ A vector of coefficients to be used in the regression spline.
#' 
#' @export
LL.captimes <- function(x_, dists_ = NULL, log.lambda.d_ = NULL, log.lambda.b.d_ = NULL, lambda0_, lambda0.b_,  const.haz_= F, gam.obj_, spline.coeff_){
  
  trap.col <- which(names(x_)=="Detector")
  cap.col <- which(names(x_)=="Caught")
  b.col=which(names(x_)=="b")		
  
  cap.trap <- as.character(x_[trap.col])
  cap.time <- as.numeric(x_[cap.col])
  
  if (is.null(log.lambda.d_)){
    dist.k <- dists_[rownames(dists_)==cap.trap,]
    sig.t <- pred.Y(time.vec_ = cap.time, gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)
    if (length(b.col)>0){
      if (as.numeric(x_[b.col])==1){
        haz.ll <- log(lambda0_*exp(-dist.k^2 / (2*rep(sig.t^2,length(dist.k)))))
      } else {
        haz.ll <- log(lambda0.b_*exp(-dist.k^2 / (2*rep(sig.t^2,length(dist.k)))))
      }
    } else {
      haz.ll <- log(lambda0_*exp(-dist.k^2 / (2*rep(sig.t^2,length(dist.k)))))
    }
    haz.ll[is.infinite(haz.ll)] <- (-1.797693e+308)/50
    LL.det <- haz.ll 
  } else {
    if (const.haz_==F){
      hazt <- pred.Y(time.vec_ = cap.time, gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)
      if (hazt==0){hazt <- 1e-100}
      if (length(b.col)>0){
        if (as.numeric(x_[b.col])==1){
          log.lambda.d.k <- log.lambda.d_[rownames(log.lambda.d_)==cap.trap,]
        } else {
          log.lambda.d.k <- log.lambda.b.d_[rownames(log.lambda.b.d_)==cap.trap,]
        }
      } else {
        log.lambda.d.k <- log.lambda.d_[rownames(log.lambda.d_)==cap.trap,]
      }
      LL.det <- log.lambda.d.k + c(rep(log(hazt),length(log.lambda.d.k)))
    } else {
      if (length(b.col)>0){
        if (as.numeric(x_[b.col])==1){
          LL.det <- log.lambda.d_[rownames(log.lambda.d_)==cap.trap,]
        } else {
          LL.det <- log.lambda.b.d_[rownames(log.lambda.b.d_)==cap.trap,]
        }
      } else {
        LL.det <- log.lambda.d_[rownames(log.lambda.d_)==cap.trap,]
      }
    }
  }
  return(LL.det)
}

#' @title Create outage intervals. 
#' @description creates intervals for the outage times that are supplied with the usage information using the intervals R package.
#' The function is used with sapply and the trap names. 
#'    
#' @return Returns the intervals that are of class "Intervals" created by the intervals package.
#' 
#' @param trap.name The trap name that is passed to this function by sapply. The trap names are extracted as rownames of the trap object.
#' @param usage_ A list that supplies information about trap usage i.e. between what time periods each detector was NOT operating. 
#' 
#' @export
intervals.outage <- function(trap.name, usage_){
  
  vec.outage <- NULL
  Trap.usage <- usage_[[trap.name_]]
  
  if (is.null(Trap.usage)){return(NULL)} else {
    
    for (i in 1:length(Trap.usage)){
      vec.outage <- c(vec.outage,Trap.usage[[i]])
    }
    
    outage.interval <- Intervals(matrix(vec.outage, ncol=2, byrow=T),closed=F)
    return(outage.interval)
  }
} 

#' @title Calcluates the exposure to capture for individuals. 
#' @description Calcluates the exposure to capture for individuals and is required for models with a behavioural effect whereby
#' each individual's exposure to detection changes depending on when they were first caught. The function also uses the usage 
#' to reduce exposure related to trap outages.
#'    
#' @return Returns a vector of total exposure.
#' 
#' @param haz.type_ The type of detection hazard to use in the model. Choices c("IndL", "IndU", "Dep") include the two 
#' independent hazard parameterisations and the dependent hazard parameterisation.
#' @param trap.obj_ An secr trap object.
#' @param first.cap_ The time of first capture for a specific individual.
#' @param usage.intervals_ Intervals of outages that are of class "Intervals". Required to subtract out periods of exposure that
#' correspond to detectors being out of action.
#' @param lambda.d_ A matrix with the encounter rates for each detector corresponding to each mask point, i.e. the rates 
#' exclude the time component. Only supplied for independent hazards, if it is NULL a dependent hazard is implied.
#' @param lambda.b.d_ Another lambda matrix corresponding to distances but using a different parameter estimate after 
#' first capture. Only relevant for models that have a behavioural effect specified.
#' @param lambda0_ The intercept parameter for the dependent hazard encounter rate function.
#' @param lambda0.b_ A second intercept parameter for the dependent hazard encounter rate function that is applicable after first
#' capture. Only relevant for models that have a behavioural effect specified.
#' @param const.haz_ Indicates if a constant hazard is used.
#' @param dists_ A matrix of distances between detectors and mask points. Used when a dependent hazard is specified.
#' @param gam.obj_ GAM obj used to find the detection hazard for the capture times.
#' @param spline.coeff_ A vector of coefficients to be used in the regression spline.
#' @param study.time_ A vector corresponding to the study duration.
#' @param setup.hour_ Can be used if the survey started after time 0 i.e. if it does not start at the start of a hazard cycle.
#' @param endpoint_ The endpoint of the study, this vzlue determines the study duration.
#' @param mask.pts_ Gives the number of rows or cells in the mask object.
#' 
#' @export
LL.survival.ind  <- function(haz.type_, trap.obj_, first.cap_, usage.intervals_ = NULL, lambda.d_ = NULL, lambda.b.d_ = NULL, lambda0_ = NULL, lambda0.b_ = NULL, const.haz_= F, dists_ = NULL, gam.obj_ = NULL, spline.coeff_ = NULL, study.time_, setup.hour_=0, endpoint_, mask.pts_){
  num.traps <- dim(trap.obj_)[1]
  traps.list <- vector("list", num.traps)
  trap.names <- rownames(trap.obj_)
  
  haz.mesh.width <- study.time_[2]-study.time_[1]
  
  #subtract matrix to be num traps * M with periods of time that need to be removed
  subtract <- matrix(0,nrow=num.traps, ncol=mask.pts_)
  rownames(subtract) <- trap.names
  
  #Create ind specific intervals around 1st capture
  int.1 <- Intervals(matrix(c(setup.hour_,first.cap_), ncol=2, byrow=T),closed=F)
  int.2 <- Intervals(matrix(c(first.cap_, endpoint_), ncol=2, byrow=T),closed=F)
  
  #calc meshes for total exposure excl outages
  total.mesh1 <- study.time_[study.time_ > setup.hour_ & study.time_ <= first.cap_]
  total.mesh2 <- study.time_[study.time_ > first.cap_ & study.time_ <= endpoint_]
  
  #Used to calculate total exposure based on b effect
  Surv.1 <- 0
  Surv.2 <- 0
  
  #1) calculate total exposure ignoring outages, based on haz type, for all traps
  if (haz.type_ == "Dep"){  #Dep hazard
    
    #find vector of sigma values that corrrespond to the two periods
    sig.vec.1 <- pred.Y(time.vec_ = total.mesh1, gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)
    sig.vec.2 <- pred.Y(time.vec_ = total.mesh2, gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)
    
    #calculate full hazard matrix, depends on both distance and time
    lambda.dt.1 = array(NA, dim =c(num.traps,length(total.mesh1),mask.pts_))
    dimnames(lambda.dt.1) <- list(trap.names,round(total.mesh1,2),NULL)
    for (k in 1:mask.pts_){
      for (j in 1:num.traps){
        lambda.dt.1[j,,k] <- lambda0_*exp(-dists_[j,k]^2 / (2*sig.vec.1^2))
      }
    }
    lambda.dt.1[is.infinite(lambda.dt.1)] = 1.797693e+308
    
    lambda.dt.2 = array(NA, dim =c(num.traps,length(total.mesh2),mask.pts_))
    dimnames(lambda.dt.2) <- list(trap.names,round(total.mesh2,2),NULL)
    for (k in 1:mask.pts_){
      for (j in 1:num.traps){
        lambda.dt.2[j,,k] <- lambda0.b_*exp(-dists_[j,k]^2 / (2*sig.vec.2^2))
      }
    }
    lambda.dt.2[is.infinite(lambda.dt.2)] = 1.797693e+308
    
    #integrate across time
    Surv.1 <- apply(lambda.dt.1,c(1,3),sum)*haz.mesh.width
    Surv.2 <- apply(lambda.dt.2,c(1,3),sum)*haz.mesh.width
    Tot.Surv <- Surv.1 + Surv.2
  } else {
    if (const.haz_==F){ #Ind hazard
      Surv.1 <- sum(pred.Y(time.vec_ = total.mesh1, gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)) * haz.mesh.width * lambda.d_
      Surv.2 <- sum(pred.Y(time.vec_ = total.mesh2, gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)) * haz.mesh.width * lambda.b.d_
      Tot.Surv <- Surv.1 + Surv.2
    } else {  #constant hazard
      Surv.1 <- (first.cap_ - setup.hour_) * lambda.d_
      Surv.2 <- (endpoint_ - first.cap_) * lambda.b.d_
      Tot.Surv <- Surv.1 + Surv.2
    }
  }
  
  #only do second part if there is some usage info
  if (length(usage.intervals_)>0){
    
    #2) calc overlaps to subtract for each trap depending on usage 
    for (trap.num in 1:num.traps){
      if (is.null(usage.intervals_[trap.names[trap.num]][[1]])) {
        subtract.1 <- 0
        subtract.2 <- 0
      } else {
        #extract times of overlap with trap k
        period.1 <- interval_intersection(int.1,usage.intervals_[[trap.names[trap.num]]])
        period.2 <- interval_intersection(int.2,usage.intervals_[[trap.names[trap.num]]])     
        
        if (haz.type_ == "Dep"){  #Dep haz
          if (length(period.1) > 0){
            subtract.1 <- 0
            for (i in 1:dim(period.1@.Data)[1]){
              mesh.start <- period.1@.Data[i,1]
              mesh.end <- period.1@.Data[i,2]
              
              #extract haz for trap over the period
              lambda.dt.1.part <- lambda.dt.1[trap.num,c(as.character(round(total.mesh1[total.mesh1 > mesh.start & total.mesh1 <= mesh.end],2))),]
              
              #integrate across time
              #if statement for cases where the overlap is at a single time point
              if (class(lambda.dt.1.part)=="numeric"){
                subtract.1 <- subtract.1 + (lambda.dt.1.part*haz.mesh.width)
              } else {
                subtract.1 <- subtract.1 + apply(lambda.dt.1.part,2,sum)*haz.mesh.width
              }
              rm(lambda.dt.1.part)  #removes array as not used again, for memory purposes
              gc()
            }
          } else {subtract.1 <- 0}
          
          if (length(period.2) > 0){
            subtract.2 <- 0
            for (i in 1:dim(period.2@.Data)[1]){
              mesh.start <- period.2@.Data[i,1]
              mesh.end <- period.2@.Data[i,2]
              
              #extract haz for trap over the period
              lambda.dt.2.part <- lambda.dt.2[trap.num,c(as.character(round(total.mesh2[total.mesh2 > mesh.start & total.mesh2 <= mesh.end],2))),]
              
              #integrate across time
              #if statement for cases where the overlap is at a single time point
              if (class(lambda.dt.2.part)=="numeric"){
                subtract.2 <- subtract.2 + (lambda.dt.2.part*haz.mesh.width)
              } else {
                subtract.2 <- subtract.2 + apply(lambda.dt.2.part,2,sum)*haz.mesh.width
              }
              
              rm(lambda.dt.2.part)  #removes array as not used again, for memory purposes
              gc()
            }
          } else {subtract.2 <- 0}
        } else {
          if (const.haz_==F){ #Ind hazard
            if (length(period.1) > 0){
              subtract.1 <- 0
              for (i in 1:dim(period.1@.Data)[1]){
                mesh.start <- period.1@.Data[i,1]
                mesh.end <- period.1@.Data[i,2]
                period1.mesh <- study.time_[study.time_ > mesh.start & study.time_ <= mesh.end]
                subtract.1 <- subtract.1 + sum(pred.Y(time.vec_ = period1.mesh, gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)) * haz.mesh.width * lambda.d_[rownames(lambda.d_)==trap.names[trap.num]]
              }
            } else {subtract.1 <- 0}
            
            if (length(period.2) > 0){
              subtract.2 <- 0
              for (i in 1:dim(period.2@.Data)[1]){
                mesh.start <- period.2@.Data[i,1]
                mesh.end <- period.2@.Data[i,2]
                period2.mesh <- study.time_[study.time_ > mesh.start & study.time_ <= mesh.end]
                subtract.2 <- subtract.2 + sum(pred.Y(time.vec_ = period2.mesh, gam.obj_ = gam.obj_, spline.coeff_ = spline.coeff_)) * haz.mesh.width * lambda.b.d_[rownames(lambda.b.d_)==trap.names[trap.num]]
              }
            } else {subtract.2 <- 0}
          } else {  #Constant hazard
            
            if (length(period.1) > 0){
              subtract.1 <- 0
              for (i in 1:dim(period.1@.Data)[1]){
                subtract.1 <- subtract.1 + ((period.1@.Data[i,2] - period.1@.Data[i,1]) * lambda.d_[rownames(lambda.d_)==trap.names[trap.num]])
              }
            } else {subtract.1 <- 0}
            
            if (length(period.2) > 0){
              subtract.2 <- 0
              for (i in 1:dim(period.2@.Data)[1]){
                subtract.2 <- subtract.2 + ((period.2@.Data[i,2] - period.2@.Data[i,1]) * lambda.b.d_[rownames(lambda.b.d_)==trap.names[trap.num]])
              }
            } else {subtract.2 <- 0}
          }
        }
      }
      tot.subtract <- subtract.1 + subtract.2
      if (all(tot.subtract == 0)){
        subtract[trap.num,] <- 0
      } else {
        subtract[trap.num,] <- tot.subtract
      }
    }
    return('Ind Exposure'=Tot.Surv - subtract) 
  } else {
    return('Ind Exposure'=Tot.Surv) 
  }
}

#' @title Calcluates the position in the beta vector for all parameters. 
#' @description The function calcluates the position in the beta vector for density, encounter rate / detection function, 
#' and spline parameters that depends on the type of model that is specified. The function also accommodates the sharing 
#' of parameters across different sessions if desired.
#'    
#' @return Returns a list with the positions in the beta vector for density and encounter rate / detection function parameters.
#' 
#' @param D.mod_ A numeric vector that specifies what density parameters are required for the different sessions. A different 
#' number specifies different parameters e.g. if there are 2 sessions: c(1,1) specifies the same density parameters for the two 
#' sessions whereas c(1,2) specifies that the 2nd session needs its own density parameters. The length of the vector must equal 
#' the number of sessions.
#' @param D.type_ A string vector that specifies the density model for each session. Options include constant, exponential 
#' and quadratic density ( c("Cons","Exp", "Quad"). The length of the vector must equal the number of sessions in the data.
#' @param lambda0.mod_ A numeric vector that specifies what intercept parameters are required for the different sessions. A different 
#' number specifies a different parameter e.g. if there are 2 sessions: c(1,1) specifies the same lambda0 parameter for the two 
#' sessions whereas c(1,2) specifies that the 2nd session needs its own lambda0 parameter. The length of the vector must equal 
#' the number of sessions.
#' @param b.effect_ A string vector  that specifies if there is a behavioural effect on the intercept of the encounter rate / 
#' detection function. Options include c("n","b") where "b" indicates the effect. The length of the vector must equal the number 
#' of sessions in the data.
#' @param sigma.mod_ A numeric vector that specifies what range parameters are required for the different sessions. A different 
#' number specifies a different parameter e.g. if there are 2 sessions: c(1,1) specifies the same sigma parameter for the two 
#' sessions whereas c(1,2) specifies that the 2nd session needs its own sigma parameter. The length of the vector must equal 
#' the number of sessions.
#' @param haz.type_ The type of detection hazard to use in the model. Choices c("IndL", "IndU", "Dep") include the two 
#' independent hazard parameterisations and the dependent hazard parameterisation.
#' @param haz.mod_ A numeric vector that specifies what spline parameters are required for the different sessions. A different 
#' number specifies different parameters e.g. if there are 2 sessions: c(1,1) specifies the same spline parameters for the two 
#' sessions whereas c(1,2) specifies that the 2nd session needs its own spline parameters. The length of the vector must equal 
#' the number of sessions.
#' @param haz.k_ The values here match what is passed to mgcv and so do not correspond with the number of spline parameters that
#' are estimated. The length of the vector must equal the number of sessions in the data. A value of 1 specifies a constant hazard 
#' and then no spline parameters are involved.
#' @param conduct.mod_ A numeric vector that specifies whether the same conductance parameter is used for the different sessions. 
#' A different number specifies that different parameters are estimated. The length of the vector must equal the number of sessions.
#' 
#' @export
Calc.pos.Haz <- function(D.mod_=c(1), D.type_=c('Cons'), lambda0.mod_= c(1), b.effect_ = c('n'), sigma.mod_ = c(1), haz.type_ = 'Dep', haz.mod_=c(1), haz.k_=c(4), conduct.mod_=NULL){
  dens.pos.vec <- NULL
  spline.pos.vec <- NULL
  lambda0.pos.vec <- NULL
  sigma.pos.vec <- NULL
  conduct.pos.vec <- NULL
  
  dens.pos <- NULL
  lambda0.pos <- NULL
  sigma.pos <- NULL
  conduct.pos <- NULL
  
  num.sessions <- length(D.mod_)
  session.ids <- 1:num.sessions
  
  for (z in 1:num.sessions){
    if (z==1){
      dens.pos <- 1
      if (D.type_[z]=="Cons"){  #constant density (1 par)
        lambda0.pos <- 2
        if (b.effect_[z] == 'n'){ #no behaviour effect specified
          if (haz.type_ != 'Dep'){
            sigma.pos <- 3
            spline.pos <- 4
          } else {  #for the dependent hazard the spline pars start earlier
            spline.pos <- 3
          }
        } else {  #with a behaviour effect
          if (haz.type_ != 'Dep'){
            sigma.pos <- 4
            spline.pos <- 5
          } else {  #for the dependent hazard the spline pars start earlier
            spline.pos <- 4
          }
        }
      } else {
        if (D.type_[z]=="Exp"){ #exp density (2 pars)
          lambda0.pos <- 3
          if (b.effect_[z] == 'n'){ #no behaviour effect specified
            if (haz.type_ != 'Dep'){
              sigma.pos <- 4
              spline.pos <- 5
            } else {  #for the dependent hazard the spline pars start earlier
              spline.pos <- 4
            }
          } else {  #with a behaviour effect
            if (haz.type_ != 'Dep'){
              sigma.pos <- 5
              spline.pos <- 6
            } else {  #for the dependent hazard the spline pars start earlier
              spline.pos <- 5
            }
          }
        } else {  #quadratic density (3 pars)
          lambda0.pos <- 4
          if (b.effect_[z] == 'n'){ #no behaviour effect specified
            if (haz.type_ != 'Dep'){
              sigma.pos <- 5
              spline.pos <- 6
            } else {  #for the dependent hazard the spline pars start earlier
              spline.pos <- 5
            }
          } else {  #with a behaviour effect
            if (haz.type_ != 'Dep'){
              sigma.pos <- 6
              spline.pos <- 7
            } else {  #for the dependent hazard the spline pars start earlier
              spline.pos <- 6
            }
          }
        }
      } #closes density type loop
      
      if (haz.type_ == 'Dep'){  #correct num spline pars, 1 df lost with dep haz and 2 otherwise
        if (is.null(conduct.mod_) ){
          counter <- spline.pos + haz.k_[z] - 1 #counter calculates the position at the end of the parameter vector, where the next parameter goes
        } else {
          conduct.pos <- spline.pos + haz.k_[z] - 1
          counter <- conduct.pos + 1
        }
      } else {  
        if (is.null(conduct.mod_)){
          if (haz.k_[z]==1){
            counter <- spline.pos
          } else {
            counter <- spline.pos + haz.k_[z] - 2            
          }
        } else {
          if (haz.k_[z]==1){
            conduct.pos <- spline.pos
            counter <- conduct.pos + 1
          } else {
            conduct.pos <- spline.pos + haz.k_[z] - 2
            counter <- conduct.pos + 1
          }
        }
      }
    } else {  #below is code for sessions other than the first
      if (any(session.ids[1:(z-1)] == D.mod_[z])){  #if density.s est = density est from previous session 
        session.D.par <- which(session.ids == D.mod_[z])
        dens.pos <- dens.pos.vec[session.D.par]
      } else {
        dens.pos <- counter
        if (D.type_[z]=='Cons'){
          counter <- counter + 1
        } else {
          if (D.type_[z]=='Exp'){
            counter <- counter + 2
          } else {
            counter <- counter + 3
          }
        }
      }
      
      if (any(session.ids[1:(z-1)] == lambda0.mod_[z])){  #if lambda0.s est = est from previous session
        session.lambda.par <- which(session.ids == lambda0.mod_[z])
        lambda0.pos <- lambda0.pos.vec[session.lambda.par]
      } else {
        lambda0.pos <- counter
        if (b.effect_[z] == 'n'){
          counter <- counter + 1
        } else {
          counter <- counter + 2
        }
      }
      
      if (haz.type_ != 'Dep'){  #for dep haz sigma is a fn of the spline pars
        if (any(session.ids[1:(z-1)] == sigma.mod_[z])){  #if sigma.s est = est from previous session
          session.sigma.par <- which(session.ids == sigma.mod_[z])
          sigma.pos <- sigma.pos.vec[session.sigma.par]
        } else {
          sigma.pos <- counter
          counter <- counter + 1
        }
      }
      
      if (any(session.ids[1:(z-1)] == haz.mod_[z])){  #if spline pars = pars from previous session
        session.spline.par <- which(session.ids == haz.mod_[z])
        spline.pos <- spline.pos.vec[session.spline.par]
      } else {
        spline.pos <- counter
        if (haz.type_ == 'Dep'){
          counter <- spline.pos + haz.k_[z] - 1
        } else {
          if (haz.k_[z]==1){
            counter <- spline.pos
          } else {
            counter <- spline.pos + haz.k_[z] - 2            
          }          
        }
      }
      
      if (any(session.ids[1:(z-1)] == conduct.mod_[z])){  #if conduct par = par from previous session
        session.conduct.par <- which(session.ids == conduct.mod_[z])
        conduct.pos <- conduct.pos.vec[session.conduct.par]
      } else {
        if (!is.null(conduct.mod_)){
          conduct.pos <- counter
          counter <- conduct.pos + 1
        }     
      }
      
    } #closes loop for sessions other than first    
    
    #constructs the vectors that record what position the various par appear in, will execute after s = 1 for first time
    dens.pos.vec <- c(dens.pos.vec, dens.pos)
    lambda0.pos.vec <- c(lambda0.pos.vec, lambda0.pos)
    sigma.pos.vec <- c(sigma.pos.vec, sigma.pos)
    spline.pos.vec <- c(spline.pos.vec, spline.pos)
    conduct.pos.vec <- c(conduct.pos.vec, conduct.pos)
  } #closes the loop over the z sessions
  
  if (haz.type_ == 'Dep'){
    if (!is.null(conduct.mod_)){
      return(list("Density" = dens.pos.vec, "Lambda0" = lambda0.pos.vec, "Spline" = spline.pos.vec, "Conductance" = conduct.pos.vec)) 
    } else {
      return(list("Density" = dens.pos.vec, "Lambda0" = lambda0.pos.vec, "Spline" = spline.pos.vec)) 
    }
    
  } else {
    if (haz.type_ == 'IndU'){
      if (!is.null(conduct.mod_)){
        return(list("Density" = dens.pos.vec, "Lambda0" = lambda0.pos.vec, "Sigma" = sigma.pos.vec, "Spline" = spline.pos.vec, "Conductance" = conduct.pos.vec)) 
      } else {
        return(list("Density" = dens.pos.vec, "Lambda0" = lambda0.pos.vec, "Sigma" = sigma.pos.vec, "Spline" = spline.pos.vec))  
      }
    } else {
      if (!is.null(conduct.mod_)){
        return(list("Density" = dens.pos.vec, "G0" = lambda0.pos.vec, "Sigma" = sigma.pos.vec, "Spline" = spline.pos.vec, "Conductance" = conduct.pos.vec)) 
      } else {
        return(list("Density" = dens.pos.vec, "G0" = lambda0.pos.vec, "Sigma" = sigma.pos.vec, "Spline" = spline.pos.vec))
      }
    }
  }
}

#' @title Calcluates the least cost distances. 
#' @description The function calcluates the least cost distance matrix for a given conductance parameter and spatial covariate.
#' Assumes one covariate attached to the mask object. exp(Beta0) is sigma and beta0 has been excluded. Requires the gdistance
#' package.
#'    
#' @return Returns a matrix containing least cost distances between all cells.
#' 
#' @param traps.obj_ A traps object that is used to get least cost distances from every trap to every mask cell.
#' @param mask.obj_ A mask object that needs to have one numeric covariate attached to it.
#' @param transitionFn_ Specifies how the conductance between two cells is calculated. Defaults to 1 for the arithmetic mean.
#' 2 specifies the geometric mean. Used in gdistance's transition function.
#' @param directions_ Specifies how many different transitions are possible from each cell. Used in gdistance's transition function.
#' @param conduct.par_ Supplies the conductance parameter that is used with the spatial copvariate to calculate the noneuc 
#' linear predictor.
#' 
#' @export
noneuc.distances <- function(traps.obj_, mask.obj_, transitionFn_ = 1, directions_ = 16, conduct.par_ = 0){
  #Calc the linear predictor called noneuc
  covariates(mask.obj_)$noneuc <- conduct.par_ * covariates(mask.obj_)[,1]
  Sraster <- raster(mask.obj_, 'noneuc')
  require(gdistance)
  if (transitionFn_==1){
    trans <- transition(Sraster, transitionFunction = function(x) mean(exp(x)),
                        directions = directions_)
  } else {
    trans <- transition(Sraster, transitionFunction =  function(x) exp(mean(x)),
                        directions = directions_)
  }
  
  trans.geo <- geoCorrection(trans)
  LC <- costDistance(trans.geo, as.matrix(traps.obj_), as.matrix(mask.obj_))
  return(LC)
}
