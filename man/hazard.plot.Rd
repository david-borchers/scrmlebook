% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CTfns.R
\name{hazard.plot}
\alias{hazard.plot}
\title{Plot cyclical detection hazard.}
\usage{
hazard.plot(cycle.hours_ = c(0, 24), time.vec_ = NULL,
  dist.vec_ = 0:1000, mesh.size_ = 100, det.pars_, spline.pars_,
  vcov_ = NULL, ci.method_ = "delta", haz.k_ = c(4),
  haz.type_ = "IndU", dim_ = "time", full.ER_ = T, add_ = F,
  dist_ = NULL, time_ = NULL, main_ = "", xlabel_ = "Cycle Hour",
  ylabel_ = "Encounter Rate", ylim_ = NULL, cex.legend_ = 1,
  legend.show_ = T, legend.pos_ = "topright", col_ = "black",
  theta_ = 30, phi_ = 15, rotate_ = F, ydecimals_ = 3)
}
\arguments{
\item{cycle.hours_}{Specifies the vector of times for one cycle of the cyclical hazard.}

\item{time.vec_}{Specifies the vector of times for the plotted hazard.}

\item{dist.vec_}{A vector of distance values that is used to generate the 3D plot.}

\item{mesh.size_}{Determines the coarseness of the plotted lines.}

\item{det.pars_}{The relevant parameters from the detection function or encounter rate function.
. If the 'IndU' hazard is specified these will include lambda0 and sigma, if 'IndL' 
then g0 replaces lambda0, and if a dependent hazard is specified ('Dep') then only
lambda0 is relevant.}

\item{spline.pars_}{The set of spline coefficients that determine the shape of the hazard.}

\item{vcov_}{A variance-covariance matrix for the parameters that is used to get std errors for confidence intervals.
It must be of the appropriate dimension and exclude sigma if the Delta method is used for the lambda0(t) model. Note
the plotting of ci's only currently available for lambda0(t) model.}

\item{ci.method_}{Specifies the method used to calculate the variance and produce intervals, choices include "delta" or "boot".
Note that a plot that uses the bootstrap approach with the sigma(t) hazard will produce slightly wonky intervals. Calc.ER.Variance
needs to be extended for multiple distances to correct this.}

\item{haz.k_}{The appropriate df. This number should correspond to the value passed to mgcv
and not to the number of spline coefficients.}

\item{haz.type_}{The chosen parameterisation (from 'IndU', 'IndL', and 'Dep').}

\item{dim_}{The dimension argument is used when plotting the dependent hazard and determines what type of plot
is drawn. If 'time' is specified, the plot will depict the hazard over the time dimension
for a given distance (specified with the dist_ argument), if 'dist' is specified the plot
will depict the hazard over distance for a given time (specified with the time_ argument),
and if nothing is specified (dim_ is null) then a plot is drawn that includes both dimensions.}

\item{full.ER_}{Specifies whether only the time component of the function is displayed, which is appropriate when 
interest is focused on the shape of the hazard through time. If TRUE the full encounter rate function is plotted 
for a given distance.}

\item{add_}{Specifies whether a new plot should be called or rather if lines should
be added to an existing plot.}

\item{dist_}{A specified value for distance. More than one value can be supplied in which case a line will be drawn for each.}

\item{time_}{A specified value for time. More than one value can be supplied in which case a line will be drawn for each.}

\item{main_}{The title for the plot.}

\item{xlabel_}{Label for the X axis.}

\item{ylabel_}{Label for the y axias.}

\item{ylim_}{Option to specfiy the scale of the Y axis.}

\item{cex.legend_}{A scaling factor to control the size of the legend.}

\item{legend.show_}{Specifies if a legend should appear in the plot. It is not relevant 
for plots showing the scaled hazard.}

\item{legend.pos_}{Specifies the position of the legend.}

\item{col_}{Specifies the colour to use for the plotted line}

\item{theta_}{defines the viewing angle, theta sets the azimuthal direction.}

\item{phi_}{defines the viewing angle, phi sets the colatitude}

\item{rotate_}{If this is set to FALSE a static 3 dimensional perspective plot is generated. If it is TRUE
the 3D plot can be rotated by the user  (requires the rgl library).}

\item{ydecimals_}{Specifies the number of decimals to use for the labelling of the Y axis}
}
\value{
Returns a plot.
}
\description{
Plots the detection hazard from a set of spline coefficients for a specified time period. The function has different ways 
of plotting different types of hazards, including a 3D plot for the dependent hazard.
}
