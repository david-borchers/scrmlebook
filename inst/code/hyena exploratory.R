# Read in and format Hyena SCR data
# =================================

library(secr)
library(sf)
library(sp)
library(scrmlebook)

CH <- read.capthist(captfile="./analysis/data/Hyeana data/capt.txt", trapfile="trap.txt", fmt="trapID", detector="single")
# summarise capthist object\par
summary (CH)
# build habitat mask\par
mask <- make.mask (traps(CH), type="trapbuffer", buffer=100, nx=64)


test1 <- st_read("./analysis/data/Hyeana data/1km.shp")
plot(test1)
buff10k <- st_read("./analysis/data/Hyeana data/buff 10k.shp")
plot(buff10k)
test3 <- st_read("./analysis/data/Hyeana data/buff 10km.shp")
plot(test3)
buffvill <- st_read("./analysis/data/Hyeana data/buffvill.shp")
plot(buffvill)
test5 <- st_read("./analysis/data/Hyeana data/fishnet.shp")
plot(test5)
test6 <- st_read("./analysis/data/Hyeana data/Hyena_distri.shp")
plot(test6)
test7 <- st_read("./analysis/data/Hyeana data/hyena-cap-prob.shp")
plot(test7)
test8 <- st_read("./analysis/data/Hyeana data/mcp of traps.shp")
plot(test8)
# slope file seems to be corrupt
#test9 <- st_read("./analysis/data/Hyeana data/slope.shp")
#plot(test9)
# slope12 file seems to be corrupt
#test10 <- st_read("./analysis/data/Hyeana data/slope12.shp")
#plot(test10)
test11 <- st_read("./analysis/data/Hyeana data/trap_location.shp")
plot(test11)
test12 <- st_read("./analysis/data/Hyeana data/trap-id.shp")
# plot geneerates error message:
plot(test12)
test13 <- st_read("./analysis/data/Hyeana data/trap-loc.shp")
plot(test13)
d2vill <- st_read("./analysis/data/Hyeana data/vill-d.shp")
plot(d2vill)
test15 <- st_read("./analysis/data/Hyeana data/village.shp")
plot(test15)

trapdf = read.csv(file="./analysis/data/Hyeana data/traps.csv")
names(trapdf) = c("Trap.ID","x","y")
cams = read.traps(data=trapdf,detector="proximity",trapID="Trap.ID")

trapdf = read.csv(file="./analysis/data/Hyeana data/Trap location.csv")
ncam = nrow(trapdf)
ndays = apply(trapdf[,-c(1:3)],1,sum)
ndaytot = sum(ndays)
ndaytot
usage = rep("",length=ncam)
for(i in 1:ncam) usage[i] = paste(trapdf[i,-(1:3)],collapse="")
trapdf = cbind(trapdf[,1:3],usage) # 61 days of camera trapping
names(trapdf)[1:3] = c("#detectorID","x","y")
write.csv(trapdf,file="./analysis/data/Hyeana data/trapfile.csv",row.names=FALSE)
# Now manually edit file to get rid of all quotation marks and replace commas with spaces. Then:
cams = read.traps(file="./analysis/data/Hyeana data/trapfile.csv",detector="proximity")



plot(cams)
mask <- make.mask(cams, type="trapbuffer", buffer=10000, nx=120)
plot(mask)
plot(cams,add=TRUE)

buff10k.spdf = as_Spatial(buff10k)
splot(buff10k.spdf)
mask = addCovariates(mask,buff10k.spdf)
pdf(file="./analysis/data/Hyeana data/buff10k.pdf")
splotcovariate(mask,covariate="BUFFERDIS",what="image",col=parula(40))
points(cams,pch="+",col="black")
dev.off()
plot(cams,add=TRUE)


# Take session info from trap data and add to capthist data
# =========================================================
trapdf = read.csv(file="./analysis/data/usable Hyena data/hyenaTraps.csv")
capthistdf = read.csv(file="./analysis/data/usable Hyena data/hyenaCH.csv")
ncam = nrow(trapdf)
for(i in 1:ncam) {
  session = trapdf$BLOCK[i]
  blockch = which(capthistdf$Detector==trapdf$TRAP_ID[i])
  if(!is.null(blockch)) capthistdf$X.Session[blockch] = session
}
# capture history data below are for traps not in trapdf
capthistdf[which(capthistdf$X.Session == 2010),]



# old code:
ndays = apply(trapdf[,-c(1:3)],1,sum)
ndaytot = sum(ndays)
ndaytot
usage = rep("",length=ncam)
for(i in 1:ncam) usage[i] = paste(trapdf[i,-(1:3)],collapse="")
trapdf = cbind(trapdf[,1:3],usage) # 61 days of camera trapping
names(trapdf)[1:3] = c("#detectorID","x","y")

usagetraps = read.csv("./analysis/data/usable Hyeana data/trapfile_with_usage.csv")
usedtraps = read.csv("./analysis/data/usable Hyeana data/trap-locations.csv")


# New shapefiles:
# --------------
mesh10k.settle <- st_read("./analysis/data/usable Hyena data/shapefiles/buffer_around-settlements-10k.shp")
plot(mesh10k.settle)

mesh10k.traps <- st_read("./analysis/data/usable Hyena data/shapefiles/buffer_arround_traps-10k.shp")
plot(mesh10k.traps)

traplocs <- st_read("./analysis/data/usable Hyena data/shapefiles/trap-locations.shp")
plot(traplocs)

villages <- st_read("./analysis/data/usable Hyena data/shapefiles/Villages.shp")
plot(villages)

