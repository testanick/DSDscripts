### read in dependencies ----------------------------------
library(geomorph)
library(car)
library(Morpho)


### set working directory at images folder  ------------------------
setwd("Documents/DSD Shape Data/")  #work computer
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/DSD Shape Data/")   #home computer

### read in source files ----------------------------
source("scripts/defineSliders.R")   #contains semi-landmark and link configurations
source("scripts/WRP_FUNCTIONS.R")   #contains some custom scripts from Will that I need
source("scripts/NDT_functions.R")   #contains custom scripts for analysis
source("scripts/dataRead_blocking.R")        #reads in raw data, performs GPA, and other useful data processing






### Plot specimens to check that everything's working properly ------------------
plotAllSpecimens(gpa.2$coords, mean = T, links = cell.link.2, plot.param = list(pt.bg = as.factor(paste(day, method = "points", sep="-")), pt.cex = 0.8, mean.bg = "gray", mean.cex = 1.2, link.col = "gray", link.lwd = 2, link.lty = 1))
plotAllSpecimens(gpa.4$coords, links = cell.link.4, mean = TRUE, plot.param = list(pt.bg = as.factor(paste(ID, method = "points", sep="-")), pt.cex = 0.8, mean.bg = "gray", mean.cex = 1.2, link.col = "gray", link.lwd = 2, link.lty = 1))

### as a reminder... this shows the landmark configuration
plotAllSpecimens(gpa.2$coords, mean = T, links = cell.link.2, label = T, plot.param = list(pt.cex = 0, mean.bg = "gray", mean.cex = 1.2, link.col = "gray", link.lwd = 2, link.lty = 1))
plotAllSpecimens(gpa.4$coords, mean = T, links = cell.link.4, label = T, plot.param = list(pt.cex = 0, mean.bg = "gray", mean.cex = 1.2, link.col = "gray", link.lwd = 2, link.lty = 1))




###plot PCA space by species ----------------------
ID2 <- as.factor(paste(day, method = "points", sep="-"))
col.gp <- rainbow(length(levels(ID2)))
names(col.gp) <- levels(ID2)
col.gp <- col.gp[match(ID2, names(col.gp))] # col.gp must NOT be a factor

plotTangentSpace(gpa.2$coords, groups = col.gp, axis1 = 1, axis2 = 2, legend = T)
plotTangentSpace(gpa.4$coords, groups = col.gp, axis1 = 1, axis2 = 2, legend = T)

###plot shape change from mean to each species ---------------------

par(mfrow=c(2, 2), 
    mar=c(0.5, 0.5, 0.5, 0.2), 
    oma = c(0.2, 3, 3, 0.2))

GP3 <- gridPar(pt.bg = "gray", pt.size = 1) 

##for 'elegans', N2
plotRefToTarget(ref.2, lsmeans.2[,,4], links = cell.link.2, mag = 3, method = "vector", gridPars = GP3)  
mtext( dimnames(lsmeans.2)[[3]][2], side = 2, font = 2, line = 0)
mtext("Two-cell", side = 3, line = 0, font = 2)
plotRefToTarget(ref.4, lsmeans.4[,,2], links = cell.link.4, mag = 3, method = "vector", gridPars = GP3)  
mtext("Four-cell", side = 3, line = 0, font = 2)


##for 'briggsae', AF16
plotRefToTarget(ref.2, lsmeans.2[,,1], links = cell.link.2, mag = 3, method = "vector", gridPars = GP3)
mtext(dimnames(lsmeans.2)[[3]][1], side=2, font=2, line = 0)
plotRefToTarget(ref.4, lsmeans.4[,,1], links = cell.link.4, mag = 3, method = "vector", gridPars = GP3)

par(mfrow=c(1,1), 
    mar=c(4,4,4,1),
    oma=c(0.5,0.5,0.5,0.5))

# plot differences in centroid for each cell by species -------
ref <- mshape( gpa.4$coords )
lsmeans <- arrayspecs(( rowsum( two.d.array( gpa.4$coords ), ID ) / as.vector( table( ID ))), 18, 2) #gp means

ABa <- apply( lsmeans[ c(1,5,8,18,4,11),,], c(3,2), mean)
ABp <- apply( lsmeans[ c(1,12,2,13,3,6,5),,], c(3,2), mean)
EMS <- apply( lsmeans[ c(5,6,10,16,9,17,8),,], c(3,2), mean)
P2  <- apply( lsmeans[ c(3,14,7,6,15,10),,], c(3,2), mean)

GP3 <- gridPar(pt.bg = "gray", pt.size = 1) 
plotRefToTarget(lsmeans.4[,,1], lsmeans.4[,,3], links = cell.link.4, mag = 1, method = "points", gridPars = GP3)  

points( ABa, pch = 21, bg = c("black","gray") )
points( ABp, pch = 21, bg = c("black","gray") ) 
points( EMS, pch = 21, bg = c("black","gray") )
points( P2, pch = 21, bg = c("black","gray") )

# procrustes distance linear model to test for differences in shape based on ID (i.e. species) -----
procD.lm(gpa.2$coords ~ day * gpa.4$Csize, iter=999, RRPP=T)
procD.lm(gpa.4$coords ~ day * gpa.4$Csize, iter=999, RRPP=T)


temp1 <- data.frame(two.d.4, gpa.4$Csize)
day1.data <- subset(temp1, day=="day1")
day2.data <- subset(temp1, day=="day2")
day3.data <- subset(temp1, day=="day3")
day4.data <- subset(temp1, day=="day4")

GetDayShape <- function(day1.data)
{    day1.shape <- as.matrix(day1.data[,c(1:(ncol(day1.data)-2))])
    day1.size <- day1.data$gpa.4.Csize
    day1.allometry.model <- lm(day1.shape ~ day1.size)
    day1.Beta <- day1.allometry.model$coefficients[2,] 
    day1.shape <- ShapeScore(t(day1.Beta), day1.shape)
    day1 <- mean(day1.shape)
    return(day1)
}

day1.shape.4 <- t(replicate(1000, getShapeByDay(land = two.d.4, gpa = gpa.4, days = "day1")))
day2.shape.4 <- t(replicate(1000, getShapeByDay(land = two.d.4, gpa = gpa.4, days = "day2")))
day3.shape.4 <- t(replicate(1000, getShapeByDay(land = two.d.4, gpa = gpa.4, days = "day3")))
day4.shape.4 <- t(replicate(1000, getShapeByDay(land = two.d.4, gpa = gpa.4, days = "day4")))

day.shapes <- data.frame(
  "day" = seq(1:4),
  "day.mean" = c(GetDayShape(day1.data), GetDayShape(day2.data), GetDayShape(day3.data), GetDayShape(day4.data)),
  "day.low" = c(quantile(day1.shape.4, 0.025), quantile(day2.shape.4, 0.025), quantile(day3.shape.4, 0.025), quantile(day4.shape.4, 0.025)),
  "day.high" = c(quantile(day1.shape.4, 0.975), quantile(day2.shape.4, 0.975), quantile(day3.shape.4, 0.975), quantile(day4.shape.4, 0.975)))

plot(day.shapes$day, day.shapes$day.mean, xlim = c(0.5, 4.5), ylim = c(-0.041, 0.041), xlab = "Day", ylab = "Shape Score", main = "Effect of Dissecting Block on Shape", xaxt = 'n')
axis(side = 1, at = seq(1:4))
for (i in 1:4)
{
  lines(c(i,i), c(day.shapes$day.low[i], day.shapes$day.high[i]))
}

### ugh... 95% CIs overlap (nearly identical)... yeah... shapescore is not different at all...





#phylogenetic ANOVA/regression for shape data
s <- "worms(N2:4.2,AF16:4.2);"
cat(s, file = "ex.tre", sep = "\n")
tree.worms <- read.tree("ex.tre", keep.multi = T)
str(tree.worms)

gdf <- geomorph.data.frame(gpa.4, phy = tree.worms)
procD.pgls( two.d.array( lsmeans.2) ~ as.character(ID), phy = tree.worms, RRPP = T)



### morphological disparity between across species ------------------------
#this shows pairwise morphological disparity amongst groups
morphol.disparity(gpa.4$coords ~ day, groups= ~day, iter=999)

#just joins together blocks of data to see if things are still significant (just having a looksee)
#change less.days back to day when finished
less.days <- factor(x = rep(NA, length(day)), levels = c("first","second"))
less.days[which(day == "day1" | day == "day2")] <- "first"
less.days[which(day == "day3" | day == "day4")] <- "second"

procD.lm(gpa.4$coords ~ less.days * gpa.4$Csize, iter=999, RRPP=T)

cell.lm <- advanced.procD.lm(gpa.4$coords 
                            ~ less.days,
                            ~ less.days * ID)
summary(cell.lm)


### modularity testing ----------------------------------
#renaming landmark data
Y.2 <- gpa.2$coords
Y.4 <- gpa.4$coords


land.2 <- rep('a', 12); land.2[c(1,4)]<-'b' #set modules based on landmarks (in this case, P2 vs the rest of the embryo)
MT.2 <- modularity.test(Y.2, land.2, CI = F, iter = 999) #run modularity test to see if above-defined modules differ from one another in shape
summary(MT.2) # Test summary

land.4 <- rep('a', 18); land.4[c(5,6)]<-'b'
MT.4 <- modularity.test(Y.4, land.4, CI = F, iter = 999)
summary(MT.4)

# Significant result implies modularity present

### comparing 4-cell ebryo cell-type-sizes ----------------
ID.nums <- data.frame("ID"=day, "num"=1:length(day))
ABa <- gpa.4$coords[c(4,11,1,5,8,18),,]
ABp <- gpa.4$coords[c(1,12,2,13,3,6,5),,]
EMS <- gpa.4$coords[c(5,6,10,16,9,17,8),,]
P2 <- gpa.4$coords[c(3,14,7,15,10,6),,]

#looks like cSize() takes x/y coords from one specimen at a time. must use mshape() to get mean shape of all specimens
cellSize_func <- function(x)
{
  temp.x <- x[,,sample(dim(x)[3], replace=T)]
  output <- cSize(mshape(temp.x))
  return(output)
}

#to change cell type, use (ABa, ABp, EMS, P2) outside of brackets; to change species, use (N2, AF16, etc...) inside quotes of ID == ""
get.cellSize <- function(cell.type, species){
  X <- cell.type[,,subset(ID.nums, ID==species)$num]
  X.mean <- cSize(mshape(X))
  X.boot <- replicate(1000, cellSize_func(X))
  X.CI <- quantile(X.boot, probs = c(0.025, .975))
  return(c(X.mean, X.CI))
}

# type <- c("ABa","ABa","ABp","ABp","EMS","EMS","P2","P2")
type <- rep(c("ABa", "ABp", "EMS", "P2"), 4)
species <- c(rep(c("day1"),4), rep(c("day2"),4), rep(c("day3"),4), rep(c("day4"),4))
# species <- rep(c("N2", "AF16"),4)
cell.points <-   rbind(
      get.cellSize(ABa, "day1"),
      get.cellSize(ABp, "day1"),
      get.cellSize(EMS, "day1"),
      get.cellSize(P2, "day1"),
      get.cellSize(ABa, "day2"),
      get.cellSize(ABp, "day2"),
      get.cellSize(EMS, "day2"),
      get.cellSize(P2, "day2"),
      get.cellSize(ABa, "day3"),
      get.cellSize(ABp, "day3"),
      get.cellSize(EMS, "day3"),
      get.cellSize(P2, "day3"),
      get.cellSize(ABa, "day4"),
      get.cellSize(ABp, "day4"),
      get.cellSize(EMS, "day4"),
      get.cellSize(P2, "day4"))
cells <- data.frame("type" = type, "species" = species, "mean" = cell.points[,1], "low" = cell.points[,2], "high" = cell.points[,3])
cells <- cells[with(cells, order(-species, type)), ]

par(mfrow=c(1,1))
plot(0,0, main = "Cell Sizes by Day", xlim = c(0.5,nrow(cells)+5), ylim = c(min(cells$low), max(cells$high)), xaxt = 'n', xlab = NA, ylab = "Centroid Size")
j <- -0.5
for (i in 1:nrow(cells)){
  if (cells$type[i] == "ABa"){ cell.pch  <- 5}
  if (cells$type[i] == "ABp"){ cell.pch  <- 0}
  if (cells$type[i] == "EMS"){ cell.pch  <- 1}
  if (cells$type[i] == "P2"){ cell.pch  <- 2}
  if (cells$species[i] == "day1") { cell.color <- "limegreen"}  #elegans
  if (cells$species[i] == "day2") { cell.color <- "blue"}     #day2
  if (cells$species[i] == "day3") { cell.color <- "bisque2"} #day3
  if (cells$species[i] == "day4") { cell.color <- "red1"}    #day4  
  
  #this part separates the species' points from one another
  if (i %% length(unique(cells$species)) == 0) {k <- 0 }
  if (i %% length(unique(cells$species)) == 1) {k <- 1 }
  j <- j + 0.5 + k
  
  points(j, cells$mean[i], pch = cell.pch, col = cell.color)
  lines(x = c(j,j), y = c(cells$low[i], cells$high[i]), pch = cell.pch, col = cell.color)
}
legend(x = 0.25, y = 0.3, c("ABa","ABp","EMS","P2"), pch = c(5,0,1,2), bty = 'n')
legend(x = 1, y = 0.3, c("elegans","briggsae", "brenneri","remanei.a","remanei.b"), pch = c(15,15,15,15,15), col = c("limegreen","blue","bisque3","red1","red4"), bty = 'n')
par(new=T, mar=c(3,8,13,8))
GP3 <- gridPar(pt.bg = "blue", link.col = "blue", pt.size = 1, tar.pt.bg = "limegreen", tar.link.col = "limegreen") 
plotRefToTarget(lsmeans.4[,,1], lsmeans.4[,,2], links = cell.link.4, mag = 1, method = "points", gridPars = GP3)  

par(mar=c(4,4,4,1))



### comparing 2-cell ebryo cell-type-sizes ----------------
AB <- gpa.2$coords[c(2,5,6,1,4,11,12),,]
P1 <- gpa.2$coords[c(1,7,8,3,9,10,4),,]


type2 <- c("AB","AB","P1","P1")
species2 <- rep(c("N2", "AF16"),2)
cell.points2 <-   rbind(
  get.cellSize(AB, "N2"),
  get.cellSize(AB, "AF16"),
  get.cellSize(P1, "N2"),
  get.cellSize(P1, "AF16"))
cells2 <- data.frame("type" = type2, "species" = species2, "mean" = cell.points2[,1], "low" = cell.points2[,2], "high" = cell.points2[,3])


plot(0,0, main = "Cell Sizes by Species", xlim = c(0.5,nrow(cells2)), ylim = c(min(cells2$low), max(cells2$high)), xaxt = 'n', xlab = NA, ylab = "Centroid Size")
j <- -0.5
for (i in 1:nrow(cells2)){
  if (cells2$type[i] == "AB"){ cell.pch  <- 5}
  if (cells2$type[i] == "P1"){ cell.pch  <- 2}
  if (cells2$species[i] == "N2") { cell.color <- "limegreen"}  #elegans
  if (cells2$species[i] == "AF16") { cell.color <- "blue"}     #briggsae
  if (cells2$species[i] == "LKC28") { cell.color <- "bisque2"} #brenneri
  if (cells2$species[i] == "EM464") { cell.color <- "red1"}    #remanei.a  
  if (cells2$species[i] == "PP219") { cell.color <- "red4"}    #remanei.b
  if (i %% length(unique(cells2$species)) == 0) {k <- 0 }
  if (i %% length(unique(cells2$species)) == 1) {k <- 1 }
  j <- j + 0.5 + k
  points(j, cells2$mean[i], pch = cell.pch, col = cell.color)
  lines(x = c(j,j), y = c(cells2$low[i], cells2$high[i]), pch = cell.pch, col = cell.color)
}
legend(x = 0.35, y = 0.525, c("AB","P1"), pch = c(5,2), bty = 'n')
legend(x = 0.75, y = 0.53, c("elegans","briggsae", "brenneri","remanei.a","remanei.b"), pch = c(15,15,15,15,15), col = c("limegreen","blue","bisque3","red1","red4"), bty = 'n')
par(new=T, mar=c(3,8,13,8))
GP3 <- gridPar(pt.bg = "blue", link.col = "blue", pt.size = 1, tar.pt.bg = "limegreen", tar.link.col = "limegreen") 
plotRefToTarget(lsmeans.2[,,1], lsmeans.2[,,2], links = cell.link.2, mag = 2, method = "points", gridPars = GP3)  
par(mar=c(4,4,4,1), mfrow=c(1,1))


#may want to look into using phylo.modularity() 
#find out how to get phylogenetic tree


### calculate Allometry-free values for shape dimorphism ---------------
dat_temp <- data.frame(two.d.4, ID)
csize <- gpa.4$Csize
shape_residuals <- lm(as.matrix(dat_temp[,1:24]) ~ csize * ID, data=dat_temp)$resid

N2_shape_mean <- colMeans(shape_residuals[ID == "N2",])
AF16_shape_mean <- colMeans(shape_residuals[ID == "AF16",])
PD((N2_shape_mean - AF16_shape_mean))  

procD.lm(shape_residuals ~ ID)


### calculate vector correlations ---------------------------------------

shape <- as.matrix(two.d.4)
fitboot <- lm(shape ~ day)
estimates <- coef(fitboot)

# treatment groups 
day1_shape <- estimates[1,]
day2_shape <- estimates[1,] + estimates[2,]
day3_shape <- estimates[1,] + estimates[3,]
day4_shape <- estimates[1,] + estimates[4,]

d1.2 <- ang.vec.abs(day1_shape, day2_shape)
d1.3 <- ang.vec.abs(day1_shape, day3_shape)
d1.4 <- ang.vec.abs(day1_shape, day4_shape)
d2.3 <- ang.vec.abs(day2_shape, day3_shape)
d2.4 <- ang.vec.abs(day2_shape, day4_shape)
d3.4 <- ang.vec.abs(day3_shape, day4_shape)

day_vectors <- data.frame(
  "groups" = c("day1", "day2","day3","day4"),
  "day1" = c(1,       d1.2[1],    d1.3[1], d1.4[1]),
  "day2" = c(d1.2[1], 1,          d2.3[1], d2.4[1]),
  "day3" = c(d1.3[1], d2.3[1],    1,       d3.4[1]),
  "day4" = c(d1.4[1], d2.4[1],    d3.4[1], 1))

### trajectory analysis -----------------------------------------------
# this analysis will be necessary for determining how much each species reacts to the change in temperature.
# performs a sort of vector-based analysis for differences in magnitude (path distance) and direction, etc...
# (to be used alongside traditional vector correlations)
TA <- trajectory.analysis(two.d.4 ~ day)
summary(TA)
plot(TA)
