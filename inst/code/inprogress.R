library(secr)
library(scrmlebook)

# Get some data to illustrate plotting
load("./analysis/objects/Grizzlybear.fit.h2.Rdata")

AIC(bear.fit.null.lambh2,bear.fit.null.sigh2,bear.fit.null.lamh2sigh2,
    bear.fit.best.lambh2,bear.fit.best.sigh2,bear.fit.best.lambh2.sigh2,
    bear.fit.best.lambnoh2.signoh2)

pdf("Grizzlybear_eg_Dline.pdf",h=5,w=10)
plotDline(bear.fit.best.lambh2.sigh2,"BEI")
dev.off()

# Try with Rishi data:
# save(Model1,Model2, file="C:/Users/Rishi/Rishi-fits.RData") 
load("/Users/dlb/Dropbox/Analyses/snowleopards/Spiti_Rishi/Rishi-fits.RData")

pdf("Rishi_Models_1_2.pdf",h=5,w=15)
layout(matrix(c(1,2),nrow=1))
plotDline(Model1,"Ruggedness")
plotDline(Model2,"WildPrey")
dev.off()


