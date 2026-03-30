library(secr)
library(scrmlebook)

# (1) 
# ===
data("leopard-sim")
summary(kruger.mask)
splotcovariate(kruger.mask,covariate="landscape") # does not work - not apparently able to cope with Factors

# (2) 
# ===
# Get rid of plotcovariate function (superceded by splotcovariate and implotcovariate)

