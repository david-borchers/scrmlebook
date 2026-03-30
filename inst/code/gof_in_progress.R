load("./analysis/objects/Kruger-popn.RData")

library(secr)
library(scrmlebook)

data("krugerdata")

par(mfrow=c(1,2))
splotcovariate(kruger.mask,covariate="habitat.cov",asp=1,what="image")
plot(kruger.ch.count, tracks=TRUE, add=TRUE) # count data
plot(kruger.cams,add=TRUE)
splotcovariate(kruger.mask,covariate="habitat.cov",asp=1,what="image")
plot(kruger.ch.bin, tracks=TRUE, add=TRUE) # binary data
plot(kruger.cams,add=TRUE)


fit.hn.l0 <- secr.fit(kruger.ch.count, mask=kruger.mask,
                      model=list(lambda0~1,sigma~1), detectfn="HHN", trace=FALSE)
fit.hr.l0 <- secr.fit(kruger.ch.count, mask=kruger.mask,
                 model=list(lambda0~1,sigma~1), detectfn="HHR", trace=FALSE)
fit.hn.g0 <- secr.fit(kruger.ch.count, mask=kruger.mask,
                 model=list(g0~1,sigma~1), detectfn="HN", trace=FALSE)
fit.hr.g0 <- secr.fit(kruger.ch.count, mask=kruger.mask,
                 model=list(g0~1,sigma~1), detectfn="HR", trace=FALSE)
AIC(fit.hn.l0,fit.hr.l0,fit.hn.g0,fit.hr.g0)
nsim = 99
gof.FT.hn.l0 = gof_simtest(fit.hn.l0,statfn=Freeman_Tukey,nbins=20,plot=TRUE,nsim=nsim,seed=1)
gof.Dd.hn.l0 = gof_simtest(fit.hn.l0,statfn=Devdf,nbins=20,plot=TRUE,nsim=nsim,seed=1)

# results from secr.test are implausible. Also checked with null model with mean
# kruger density and results there implausible too. Something not right with
# secr.test, so don't use it.
#fit.hn.l0.gofdbn = secr.test(fit.hn.l0,fit=TRUE,nsim=nsim)
#1-fit.hn.g0.gofdbn$output$p
#plot(fit.hn.g0.gofdbn)


bigfit  <- secr.fit(kruger.ch.count,mask=kruger.mask,model=list(D~habitat.cov, lambda0~h2, sigma~h2),
                    hcov="sex", detectfn="HHN", trace=FALSE)
nsim = 99
# This does not work - simulations don't cope with this model
#bigfit.gofdbn = secr.test(bigfit,fit=TRUE,nsim=nsim)
#1-bigfit.gofdbn$output$p
#plot(bigfit.gofdbn)
biggof.Dd = gof_simtest(bigfit,statfn=Devdf,nbins=20,plot=TRUE,nsim=99,seed=1)


# check gof test for uniform density population
fit = df.hn.fit
meanD = derived(fit)["D","estimate"]
ch = fit$capthist
dets = traps(ch)
mask = kruger.mask
set.seed(1)
simpop = sim.popn(D=meanD, core=mask, buffer=0, Ndist="poisson")
simch = sim.capthist(dets, popn=simpop, 
                     detectfn=fit$detectfn, detectpar(fit),
                     noccasions=1)
simfit = secr.fit(simch,model=fit$model,mask=mask,detectfn=fit$detectfn,trace=0)

nsim = 500
simfit.FT = FT_test(simfit,mask,plot=TRUE,nsim=nsim,seed=1)
# Monte Carlo deviance distributions
simfit.gofdbn = secr.test(df.hn.fit,fit=TRUE,nsim=nsim)
plot(simfit.gofdbn)
1-simfit.gofdbn$output$p
simfit.FTgofdbn = secr.test(df.hn.fit,statfn=Freeman_Tukey, fit=TRUE,nsim=nsim)
plot(simfit.FTgofdbn)


# testing gof stuff
####################
fit=simfit


gof.FT = gof_simtest(fit,statfn=Freeman_Tukey,nbins=20,plot=TRUE,nsim=99,seed=1)
gof.Dd = gof_simtest(fit,statfn=Devdf,nbins=20,plot=TRUE,nsim=99,seed=1)
biggof.Dd = gof_simtest(bigfit,statfn=Devdf,nbins=20,plot=TRUE,nsim=99,seed=1)

