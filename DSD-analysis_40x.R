### read in dependencies ----------------------------------
library(geomorph)
library(car)
library(Morpho)

### set working directory at images folder  ------------------------
setwd("Documents/DSD Shape Data/")  #work computer
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/DSD Shape Data/")   #home computer

### read in source files ----------------------------
source("scripts/defineSliders.R")   #contains semi-landmark and link configurations
source("scripts/NDT_functions.R")   #contains custom scripts for analysis
source("scripts/WRP_FUNCTIONS.R")
source("scripts/dataRead_40x.R")        #reads in raw data, performs GPA, and other useful data processing

par(mfrow=c(1,1))
### Plot specimens to check that everything's working properly ------------------
plotAllSpecimens(gpa.cell$coords, links = cell.link.4, mean = TRUE, plot.param = list(pt.bg = as.factor(paste(ID, devTemp, sep="-")), pt.cex = 0.8, mean.bg = "gray", mean.cex = 1.2, link.col = "gray", link.lwd = 2, link.lty = 1))

### as a reminder... this shows the landmark configuration
plotAllSpecimens(gpa.cell$coords, mean = T, links = cell.link.4, label = T, plot.param = list(pt.cex = 0, mean.bg = "gray", mean.cex = 1.2, link.col = "gray", link.lwd = 2, link.lty = 1))




###plot PCA space by species ----------------------
ID2 <- as.factor(paste(ID, devTemp, sep="-"))
col.gp <- rainbow(length(levels(ID2))); names(col.gp) <- levels(ID2)
col.gp <- col.gp[match(ID2, names(col.gp))] # col.gp must NOT be a factor

plotTangentSpace(gpa.cell$coords, groups = col.gp, axis1 = 1, axis2 = 2, legend = T)


#displays all of the min (gray) vs max (black) shapes associated with each Principal Component
par(mfrow=c(4,3), mar=c(1,1,1,0))
k <- 0
for (i in seq(1,24,2))
{
  j <- i
  j2 <- j+1
  plotRefToTarget(PCAshapes$pc.shapes[j][[1]], PCAshapes$pc.shapes[j2][[1]], method = "points", links = cell.link.4)
  k <- k+1
  mtext(side = 2, text = paste("PC", k , sep=""), line = -2)
}
par(mfrow=c(1,1), mar=c(4,4,4,1))


### Squish/buffer test ------------
ID3 <- ID2[ID2=="N2-20" | ID2=="N2-BM" | ID2=="N2-M9"]
gpa.squish <- gpagen(arrayspecs(subset(two.d.cell, ID == "N2" & devTemp != "30"), p = 18, k = 2), curves = curvepts.4, ProcD = T)

plotAllSpecimens(gpa.squish$coords, links = cell.link.4, mean = TRUE, plot.param = list(pt.bg = as.factor(ID3), pt.cex = 0.8, mean.bg = "gray", mean.cex = 1.2, link.col = "gray", link.lwd = 2, link.lty = 1))

col.gp3 <- rainbow(length(levels(ID3))); names(col.gp3) <- levels(ID3)
col.gp3 <- col.gp3[match(ID3, names(col.gp3))] # col.gp must NOT be a factor

par(mfrow=c(1,2))
PCAshapes <- plotTangentSpace(gpa.cell$coords, groups = col.gp, axis1 = 1, axis2 = 2, legend = T)
squishPCA <- plotTangentSpace(gpa.squish$coords, groups = col.gp3, axis1 = 1, axis2 = 2, legend = T)
par(mfrow=c(1,1))


procD.lm(gpa.squish$coords ~ ID3)
pc.squish <- data.frame(prcomp(subset(two.d.cell, ID=="N2" & devTemp != "30"))$x[,1:26], ID3) #26 is the point at which >99.9% of cumulative variance is explained
k.squish <- kmeans(pc.squish[,1:2], 3, 99, 20)


par(mfrow=c(1,2))
plot(pc.squish[,1], pc.squish[,2], col = as.factor(col.gp3), pch=20, main = "PCA")
plot(pc.squish[,1], pc.squish[,2], col = as.factor(k.squish$cluster), pch=20, main = "k-means")
par(mfrow=c(1,1))



#compare actual clusters with kmeans clustering ------------
pcCluster <- kmeans(pc.cell[,1:2], 5, 99, 20)

par(mfrow=c(1,2))
plot(pc.cell[,1], pc.cell[,2], col = as.factor(col.gp), pch=20, main = "PCA")
plot(pc.cell[,1], pc.cell[,2], col = as.factor(pcCluster$cluster), pch=20, main = "k-means")
par(mfrow=c(1,1))

###plot shape change from mean to each species ---------------------

par(mfrow=c(2, 2), 
    mar=c(0.5, 0.5, 0.5, 0.2), 
    oma = c(0.2, 3, 3, 0.2),
    las = 1)

GP3 <- gridPar(pt.bg = "gray", pt.size = 1) 

##for 'briggsae', AF16
plotRefToTarget(ref, lsmeans[,,1], links = cell.link.4, mag = 3, method = "TPS", gridPars = GP3)  
mtext( dimnames(lsmeans)[[3]][1], side = 3, font = 2, line = 0)
plotRefToTarget(ref, lsmeans[,,2], links = cell.link.4, mag = 3, method = "TPS", gridPars = GP3)  
mtext( dimnames(lsmeans)[[3]][2], side = 3, font = 2, line = 0)
##for 'elegans', N2
plotRefToTarget(ref, lsmeans[,,3], links = cell.link.4, mag = 3, method = "TPS", gridPars = GP3)  
mtext( dimnames(lsmeans)[[3]][3], side = 3, font = 2, line = 0)
plotRefToTarget(ref, lsmeans[,,4], links = cell.link.4, mag = 3, method = "TPS", gridPars = GP3)  
mtext( dimnames(lsmeans)[[3]][4], side = 3, font = 2, line = 0)

par(mfrow=c(1,1), 
    mar=c(4,4,4,1),
    oma=c(0.5,0.5,0.5,0.5))

# plot differences in centroid for each cell by species -------
ref <- mshape( gpa.cell$coords )

ABa <- apply( lsmeans[ c(1,5,8,18,4,11),,], c(3,2), mean)
ABp <- apply( lsmeans[ c(1,12,2,13,3,6,5),,], c(3,2), mean)
EMS <- apply( lsmeans[ c(5,6,10,16,9,17,8),,], c(3,2), mean)
P2  <- apply( lsmeans[ c(3,14,7,6,15,10),,], c(3,2), mean)

GP3 <- gridPar(pt.bg = "gray", pt.size = 1) 
plotRefToTarget(lsmeans.4[,,1], lsmeans.4[,,2], links = cell.link.4, mag = 1, method = "points", gridPars = GP3)  

points( ABa, pch = 21, bg = c("black","gray") )
points( ABp, pch = 21, bg = c("black","gray") ) 
points( EMS, pch = 21, bg = c("black","gray") )
points( P2, pch = 21, bg = c("black","gray") )

# procrustes distance linear model to test for differences in shape based on ID (i.e. species) -----

procD.lm(gpa.cell$coords ~ ID * devTemp * gpa.cell$Csize, iter=999, RRPP=T)
advanced.procD.lm(gpa.cell$coords ~ ID + devTemp,
                  ~ ID * devTemp)

temp1 <- data.frame(two.d.cell, gpa.cell$Csize)
N2.data <- subset(temp1, ID=="N2")
AF16.data <- subset(temp1, ID=="AF16")

N2.shape <- as.matrix(N2.data[,c(1:(ncol(N2.data)-2))])
N2.size <- N2.data$gpa.cell.Csize
N2.allometry.model <- lm(N2.shape ~ N2.size)
N2.Beta <- N2.allometry.model$coefficients[2,] 
N2.shape <- ShapeScore(t(N2.Beta), N2.shape)
N2 <- mean(shape)
N2.shape.4 <- t(replicate(1000, getShape(land = two.d.cell, gpa = gpa.cell, id = "N2")))
AF16.shape.4 <- t(replicate(1000, getShape(land = two.d.cell, gpa = gpa.cell, id = "AF16")))
quantile(N2.shape.4, c(0.025, 0.975))
quantile(AF16.shape.4, c(0.025, 0.975))
### need to update this section to do (N2-20 - N2-30), (AF16-20 - AF16-30)

#phylogenetic ANOVA/regression for shape data
s <- "worms(N2:4.2,AF16:4.2);"
cat(s, file = "ex.tre", sep = "\n")
tree.worms <- read.tree("ex.tre", keep.multi = T)
str(tree.worms)

gdf <- geomorph.data.frame(gpa.cell, phy = tree.worms)
procD.pgls( two.d.array( lsmeans.2) ~ as.character(ID), phy = tree.worms, RRPP = T)



### morphological disparity between across species ------------------------
#this shows pairwise morphological disparity amongst groups
morphol.disparity(gpa.cell$coords ~ ID * devTemp, groups= ~ID * devTemp, iter=999)

cell.lm <- advanced.procD.lm(gpa.cell$coords 
                             ~ ID,
                             ~ ID * devTemp)
summary(cell.lm)


### modularity testing ----------------------------------
#renaming landmark data
Y.2 <- gpa.2$coords
Y.4 <- gpa.cell$coords


land.2 <- rep('a', 12); land.2[c(1,4)]<-'b' #set modules based on landmarks (in this case, P2 vs the rest of the embryo)
MT.2 <- modularity.test(Y.2, land.2, CI = F, iter = 999) #run modularity test to see if above-defined modules differ from one another in shape
summary(MT.2) # Test summary

land.4 <- rep('a', 18); land.4[c(5,6)]<-'b'
MT.4 <- modularity.test(Y.4, land.4, CI = F, iter = 999)
summary(MT.4)

# Significant result implies modularity present

### comparing 4-cell ebryo cell-type-sizes ----------------
ID.nums <- data.frame("ID"=ID, "num"=1:length(ID))
ABa <- gpa.cell$coords[c(4,11,1,5,8,18),,]
ABp <- gpa.cell$coords[c(1,12,2,13,3,6,5),,]
EMS <- gpa.cell$coords[c(5,6,10,16,9,17,8),,]
P2 <- gpa.cell$coords[c(3,14,7,15,10,6),,]

#looks like cSize() takes x/y coords from one specimen at a time. must use mshape() to get mean shape of all specimens
cellSize_func <- function(x)
{
  temp.x <- x[,,sample(dim(x)[3], replace=T)]
  output <- cSize(mshape(temp.x))
  return(output)
}

#to change cell type, use (ABa, ABp, EMS, P2) outside of brackets; to change species, use (N2, AF16, etc...) inside quotes of ID == ""
get.cellSize <- function(cell.type, species, Tmp){
  X <- cell.type[,,subset(ID.nums, ID==species & devTemp==Tmp)$num]
  X.mean <- cSize(mshape(X))
  X.boot <- replicate(1000, cellSize_func(X))
  X.CI <- quantile(X.boot, probs = c(0.025, .975))
  return(c(X.mean, X.CI))
}

type <- c(rep("ABa", 4), rep("ABp", 4), rep("EMS", 4), rep("P2",4))
species <- rep(c("N2", "AF16"),8)
temp <- rep(c(rep("20C", 2), rep("30C",2)),4)
cell.points <-   rbind(
  get.cellSize(ABa, "N2", "20"),
  get.cellSize(ABa, "AF16", "20"),
  get.cellSize(ABa, "N2", "30"),
  get.cellSize(ABa, "AF16", "30"),
  get.cellSize(ABp, "N2", "20"),
  get.cellSize(ABp, "AF16", "20"),
  get.cellSize(ABp, "N2", "30"),
  get.cellSize(ABp, "AF16", "30"),
  get.cellSize(EMS, "N2", "20"),
  get.cellSize(EMS, "AF16", "20"),
  get.cellSize(EMS, "N2", "30"),
  get.cellSize(EMS, "AF16", "30"),
  get.cellSize(P2, "N2", "20"),
  get.cellSize(P2, "AF16", "20"),
  get.cellSize(P2, "N2", "30"),
  get.cellSize(P2, "AF16", "30")
)
cells <- data.frame("type" = type, "species" = species, "temp" = temp, "mean" = cell.points[,1], "low" = cell.points[,2], "high" = cell.points[,3])
cells <- cells[order(cells$type, cells$species),]
rownames(cells) <- seq(1:nrow(cells))

par(mfrow=c(1,1), mar=c(4,4,4,0))
plot(0,0, main = "Cell Sizes by Species", xlim = c(0,length(unique(cells$type))*2), ylim = c(min(cells$low), max(cells$high)), xaxt = 'n', xlab = NA, ylab = "Centroid Size")
j <- 0
for (i in 1:nrow(cells)){
  if (cells$type[i] == "ABa"){ cell.pch  <- 23}
  if (cells$type[i] == "ABp"){ cell.pch  <- 22}
  if (cells$type[i] == "EMS"){ cell.pch  <- 21}
  if (cells$type[i] == "P2"){ cell.pch  <- 24}
  # if (cells$temp[i] == "20C"){ cell.bg  <- "white"}
  # if (cells$temp[i] == "30C"){ cell.bg  <- "gray"}
  # if (cells$species[i] == "N2") { cell.color <- "limegreen"}  #elegans
  # if (cells$species[i] == "AF16") { cell.color <- "blue"}     #briggsae
  if (cells$species[i] == "N2" & cells$temp[i] == "20C"){ cell.bg  <- "cyan"}
  if (cells$species[i] == "N2" & cells$temp[i] == "30C"){ cell.bg  <- "purple"}
  if (cells$species[i] == "AF16" & cells$temp[i] == "20C") { cell.bg <- "red"} 
  if (cells$species[i] == "AF16" & cells$temp[i] == "30C") { cell.bg <- "green"}    
  if (i %% length(unique(cells$type)) == 0) {k <- 1 }
  if (i %% length(unique(cells$type)) == 1) {k <- 0 }
  lines(x = c(j,j), y = c(cells$low[i], cells$high[i]), pch = cell.pch, col = "black", lwd = 2)
  points(j, cells$mean[i], pch = cell.pch, col = "black", cex = 2, bg = cell.bg)
  j <- j + 0.25 + k
}
legend(x = 0.25, y = 0.3, legend = unique(cells$type), pch = c(23,22,21,24), bty = 'n')
legend(x = 2, y = 0.3, legend = c("AF16-20˚C","AF16-30˚C","N2-20˚C","N2-30˚C"), pch = c(15,15,15,15), col = c("red","green","cyan","purple"), bty = 'n')
# legend(x = 2, y = 0.285, legend = unique(cells$temp), fill = c("white","gray"), bty = 'n')
par(new=T, mar=c(15,20,5,0))
GP3 <- gridPar(pt.bg = "red", link.col = "red", pt.size = 1, tar.pt.bg = "cyan", tar.link.col = "cyan") 
plotRefToTarget(lsmeans[,,1], lsmeans[,,3], links = cell.link.4, mag = 2, method = "points", gridPars = GP3)  
par(new=T, mar=c(15,20,5,0))
GP3 <- gridPar(pt.bg = "green", link.col = "green", pt.size = 1, tar.pt.bg = "purple", tar.link.col = "purple") 
plotRefToTarget(lsmeans[,,2], lsmeans[,,4], links = cell.link.4, mag = 2, method = "points", gridPars = GP3)  
par(mar=c(4,4,4,1))




### calculate Allometry-free values for shape dimorphism ---------------
dat_temp <- data.frame(two.d.cell, ID, devTemp)
csize <- gpa.cell$Csize
shape_residuals <- lm(as.matrix(dat_temp[,1:24]) ~ ID * devTemp, data=dat_temp)$resid

N2_shape_mean <- colMeans(shape_residuals[ID == "N2",])
AF16_shape_mean <- colMeans(shape_residuals[ID == "AF16",])
PD((N2_shape_mean - AF16_shape_mean))  

procD.lm(shape_residuals ~ ID * devTemp)
advanced.procD.lm(shape_residuals ~ ID + devTemp,
                  ~ ID * devTemp,
                  RRPP = T)


### calculate vector correlations ---------------------------------------

shape <- as.matrix(two.d.cell)
fitboot <- lm(shape ~ ID*devTemp)
estimates <- coef(fitboot)

# treatment groups 
AF16_20_shape <- estimates[1,]
N2_20_shape <- estimates[1,] + estimates[2,]
AF16_30_shape <- estimates[1,] + estimates[,3]
N2_30_shape <- estimates[1,] + estimates[2,] + estimates[3,] + estimates[4,]

ang.vec.abs((AF16_20_shape - N2_20_shape), (AF16_30_shape - N2_30_shape))   # effect of temperature
ang.vec.abs((AF16_20_shape - AF16_30_shape), (N2_20_shape - N2_30_shape))   # effect of species
ang.vec.abs(AF16_20_shape, AF16_30_shape)   # effect of temp on briggsae
ang.vec.abs(N2_20_shape, N2_30_shape)   # effect of temp on N2


### trajectory analysis -----------------------------------------------
# this analysis will be necessary for determining how much each species reacts to the change in temperature.
# performs a sort of vector-based analysis for differences in magnitude (path distance) and direction, etc...
# (to be used alongside traditional vector correlations)
TA <- trajectory.analysis(two.d.cell ~ ID * devTemp)
summary(TA) # all model parameters, p <= 0.003 | principal vector correlation, p = 0.004 | path distance, p = 0.873
plot(TA)

pc1.means <- by(pc.cell[,1], pc.cell[,27:28], mean)
pc2.means <- by(pc.cell[,2], pc.cell[,27:28], mean)

par(mfrow=c(1,2))
plot(pc.cell[,1], pc.cell[,2], col = as.factor(col.gp), pch=20, main = "Raw PC Trajectories", cex = 0.75)
lines(pc1.means[c(1,3)], pc2.means[c(1,3)], col = "turquoise", lwd= 2); lines(pc1.means[c(2,4)], pc2.means[c(2,4)], col = "maroon", lwd= 2)
points(pc1.means, pc2.means, bg = c("blue","black","green","red"), pch = c(21,21, 23,23), cex = 1.5)
legend("topleft", c("AF16-20","AF16-30","N2-20","N2-30") , pt.bg = c("blue","green","black","red"), pch = c(21,23,21,23))
plot(TA)
par(mfrow=c(1,1))

### Linear Discriminant Analysis (with Jacknifed Prediction) ------
library(MASS)
treat <- paste(ID, devTemp, sep="-")
fit <- lda(treat ~ two.d.cell, 
           na.action="na.omit", CV=TRUE)
fit # show results
# Assess the accuracy of the prediction
# percent correct for each category of G
ct <- table(treat, fit$class)
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))

plot(fit)
plot(fit, dimen=1, type="both") # fit from lda
partimat(treat ~ two.d.cell, method="lda")

source("scripts/DSD-NeuralNetAnalysis.R") #this has the strata_4var in it. just run by hand for now or copy/paste here later

#single LDA on all wings, no repetitions
#training and test set
wings_genotype_train_test <- strata_4var(average_wings[,c(1:71)], "genotype", tabw_g, 2/3, 1, 71) #variable position and table length are different for the averaged dataset
wings_wg_training <- wings_genotype_train_test[[1]]
wings_wg_test <- wings_genotype_train_test[[2]]
#model
linDiscrim_all <- lda(formula=genotype ~.,data = wings_wg_training, tol = 1.0e-8, CV=FALSE) #modified tolerance
all_lda_training_table <- table(actual = wings_wg_training$genotype,
                                predicted = predict(linDiscrim_all, newdata=wings_wg_training)$class)
all_lda_training_table
all_lda_test_table <- table(actual = wings_wg_test$genotype,
                            predicted = predict(linDiscrim_all, newdata=wings_wg_test)$class)
all_lda_test_table
#accuracy rate
100*sum(diag(all_lda_test_table)/sum(all_lda_test_table))

aLDA <-  data.frame(wings_wg_training, as.matrix(wings_wg_training[,2:37]) %*% linDiscrim_all$scaling ) #matrix-multiplying the scaling matrix by the LM's gives a vector of individual LD 'scores' for each LDA
genocols <- c( "#00008B99", "#B8860B99", "#00000099", "#00640099", "#8B000099")
genopch <- c( 21, 25, 23, 24, 22 )
par( mfrow= c(1, 2), mar=c(5,5,1,1))
plot( aLDA$LD1, aLDA$LD2, pch=genopch[aLDA$genotype], bg=genocols[aLDA$genotype], col=genocols[aLDA$genotype], xlab = "1st Discriminant Function", ylab = "2nd Discriminant Function", cex=1.5, main = "LDA Training Data", xlim = c(-50,-10), ylim = c(0,40))
# x = lda 1, y = lda2
legend("topright", bty = "n", pch = genopch, pt.bg=c("#00008B", "#B8860B", "#000000", "#006400", "#8B0000"), col=c("#00008B", "#B8860B", "#000000", "#006400", "#8B0000"), Names, xjust=0 )


#plotting LDA results, test set only
bLDA <-  data.frame(wings_wg_test, as.matrix(wings_wg_test[,2:37]) %*% linDiscrim_all$scaling ) #matrix-multiplying the scaling matrix by the LM's gives a vector of individual LD 'scores' for each LDA
genocols <- c( "#00008B99", "#B8860B99", "#00000099", "#00640099", "#8B000099")
genopch <- c( 21, 25, 23, 24, 22 )
# par( mfrow= c(1, 1), mar=c(5,5,1,1))
plot( bLDA$LD1, bLDA$LD2, pch=genopch[bLDA$genotype], col=genocols[bLDA$genotype], bg=genocols[bLDA$genotype], xlab = "1st Discriminant Function", ylab = "2nd Discriminant Function", cex=1.5, main= "LDA Test Data", xlim = c(-50,-10), ylim = c(0,40) ) # x = lda 1, y = lda2

legend("topright", bty = "n", pch = genopch, pt.bg=c("#00008B", "#B8860B", "#000000", "#006400", "#8B0000"), col=c("#00008B", "#B8860B", "#000000", "#006400", "#8B0000"), Names, xjust=0)



