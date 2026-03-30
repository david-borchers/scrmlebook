# log-likelihood for binary proximity detector data
# =================================================

hn=function(d,sigma2=1,l0=1) l0*exp(-0.5*d^2/sigma2)
p4hn=function(d,sigma2=1,l0=1,T=1) 1-exp(-l0*exp(-0.5*d^2/sigma2)*T)


p.hn = function(x,y,sigma2,l0,traps){
  dists = scrmlebook:::distances(matrix(c(x,y),nrow=1,ncol=2),as.matrix(traps[,c("x","y")]))
  p = p4hn(dists,sigma2,l0)
  return(p)
}

negllik = function(theta,count,nocc,traps) {
  # unpack the parameters
  nind = nrow(count)
  x = theta[1:nind]
  y = theta[(nind+1):(2*nind)]
  sigma2 = exp(theta[2*nind+1])
  l0 = exp(theta[2*nind+2])
  
  # loop through individuals
  nll = 0
  for(i in 1:nind) {
    pcount = p.hn(x[i],y[i],sigma2,l0,traps)
    binp = dbinom(count[i,],nocc,pcount)
    logbinp = log(binp)
    if(is.nan(sum(logbinp))) {
      cat("Par: ",theta,"\n")
      cat("logbinp:  ",logbinp,"\n")
      cat("Individual: ",i,"\n")
      stop("Stopped because of NaN")
    }
    nll = nll - sum(binp)
  }
  return(nll)
}


getfitx = function(fit) {
  nx = length(fit$par) - 2
  estx = matrix(fit$par[1:nx],ncol=2)
  row.names(estx) = 1:(nx/2)
  colnames(estx) = c("x","y")
  return(estx)
}

#================== end of functions ==================

# testing stuff
xk = x[1,1]
yk = x[1,2]
sigma2 = 0.3
l0 = 1.5
  
p.hn(xk,yk,sigma2,l0,simtraps)

dbinom(count[1,],nocc,p.hn(xk,yk,sigma2,l0,simtraps),log=TRUE)

x1=c(2.4,2.6) # individual activity centre
x2=c(4.1,4.2) # individual activity centre
x3=c(0.71,0.76) # individual activity centre
x4=c(1.52,4.9) # individual activity centre
x = rbind(x1,x2,x3,x4)
colnames(x) = c("x","y")

theta = c(x[,1],x[,2],log(sigma2),log(l0))
fit = optim(theta,negllik,count=count,nocc=nocc,traps=simtraps,hessian=TRUE)
vcv = solve(fit$hessian)

estx = getfitx(fit)

plot(simtraps,border=6)
points(x[,1],x[,2])
points(estx[,1],estx[,2],pch="*")




# try just the first individual
theta1 = c(x[1,1],x[1,2],log(sigma2),log(l0))
fit1 = optim(theta1,negllik,count=count[1,,drop=FALSE],nocc=nocc,traps=simtraps,hessian=TRUE)
vcv = solve(fit$hessian)

estx = getfitx(fit)

plot(simtraps,border=6)
points(x[1,1],x[1,2])
points(estx[1,1],estx[1,2],pch="*")


# Try with eyeballed starting values for AC locations:
theta = c(c(2.4,4.1,0.71,1.52),c(2.6,4.2,0.76,4.9),log(sigma2),log(l0))
#theta = c(rep(2.5,4),rep(2.5,4),log(sigma2),log(l0)) # does not - need good starting values!
fit = optim(theta,negllik,count=count,nocc=nocc,traps=simtraps,hessian=TRUE,control=list(trace=5))
#nlmfit = nlm(negllik,theta,count=count,nocc=nocc,traps=simtraps,hessian=TRUE)

vcv = solve(fit$hessian)
#nlmvcv = solve(nlmfit$hessian)

estx = getfitx(fit)

plot(simtraps,border=6)
points(x[,1],x[,2])
points(estx[,1],estx[,2],pch="*")

# Try plotting ACs and det prob with estimated parameters:

sigma2hat = exp(fit$par[9])
l0hat = exp(fit$par[10])

plotppn = FALSE # true if want to plot heatmap of proportion of time at location, rather than of detection probs
# plot the detection probabilities
trimmaskEn = trimmaskprobs =trimmaskdists = matrix(rep(0,4*ntrimmask),nrow=4)
if(dopdf) pdf("./keepfigure/spatialdetprobest.pdf",h=10,w=10)
par(mfrow=c(2,2))
for(i in 1:4) {
  # distances from ac to mask
  trimmaskdists[i,] = scrmlebook:::distances(matrix(estx[i,],nrow=1,ncol=2),trimmask) # distances from xi
  # detection probabilities on mask
  trimmaskprobs[i,] = p4hn(trimmaskdists[i,],sigma2=sigma2hat,l0=l0hat,T=sT)
#  # nocc-occasion probs:
#  trimmaskprobs[i,] = 1 - (1 - trimmaskprobs[i,])^nocc
  # expected encounters
  trimmaskEn[i,] = hn(trimmaskdists[i,],sigma2=sigma2hat,l0=l0hat)*sT
  # put coordinates and probs in a data frame:
  if(plotppn) {
    totEn = sum(trimmaskEn[1,]) # animal i=1 is wholly within mask, so we get integral of En from it
    hframe=data.frame(x=trimmask$x,y=trimmask$y,z=trimmaskEn[i,]/totEn)
    hframe$z = hframe$z * 100 # *100 to make it a percentage
    if(i==1) zlim = max(hframe$z) # max for animal wholly within mask
    prob4image = prep4image(hframe,breaks=seq(0,zlim,length=26),col=viridis(25))
  }else {
    hframe=data.frame(x=trimmask$x,y=trimmask$y,z=trimmaskprobs[i,])
    prob4image = prep4image(hframe,breaks=seq(0,1,length=26),col=viridis(25))
  }
  points(estx[i,1],estx[i,2],pch=19,col="white")
  points(estx[i,1],estx[i,2])
  plot(simtraps,border=0,add=TRUE)
  title(paste("Individual ",i,sep=""))
}
if(dopdf) dev.off()
