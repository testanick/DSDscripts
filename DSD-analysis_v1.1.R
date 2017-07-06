library(geomorph)

#set working directory at images folder
#setwd("Documents/Rhevolution_4cell_images/")
setwd("Documents/DSD Shape Data/")
source("scripts/defineSliders.R")


#skip everything above and use this if you already have one .tps file with all landmarks on it
AF16.2 <- two.d.array(readland.tps("images/2-cell/AF16-2.tps", specID = "imageID"))
AF16.4 <- two.d.array(readland.tps("images/4-cell/AF16-4.tps", specID = "imageID"))
N2.2 <- two.d.array(readland.tps("images/2-cell/N2-2.tps", specID = "imageID"))
N2.4 <- two.d.array(readland.tps("images/4-cell/N2-4.tps", specID = "imageID"))

cells.4 <- arrayspecs( rbind(AF16.4, N2.4), p=18, k=2) 
cells.2 <- arrayspecs( rbind(AF16.2, N2.2), p=12, k=2) 



#looks like [,,22] has NAs
gpa.4 <- gpagen(cells.4, curves = curvepts.4, ProcD = T, max.iter = T)
gpa.2 <- gpagen(cells.2, curves = curvepts.2, ProcD = T, max.iter = T)

##turns data back to 2D
twoD.4<- two.d.array(gpa.4$coords)
twoD.2<- two.d.array(gpa.2$coords)

#pulls out data from tps files to allow for grouping
ID.rough <- strsplit(rownames(twoD.4), "_")
ID <- as.factor( matrix(unlist(ID.rough), ncol=2, byrow=T)[,1])





#### visualize the data (by species) #####
cell.link.2<-matrix(c(2,5,
                      5,6,
                      6,1,
                      1,7,
                      7,8,
                      8,3,
                      3,9,
                      9,10,
                      10,4,
                      4,11,
                      11,12,
                      12,2,
                      1,4
                      ),ncol=2,byrow=T)
cell.link.4<-matrix(c(4,11,
                    11,1,
                    1,12,
                    12,2,
                    2,13,
                    13,3,
                    3,14,
                    14,7,
                    7,15,
                    15,10,
                    10,16,
                    16,9,
                    9,17,
                    17,8,
                    8,18,
                    18,4,
                    1,5,
                    8,5,
                    5,6,
                    6,3,
                    6,10
                  ),ncol=2,byrow=T)


plotAllSpecimens(gpa.4$coords, mean=T, links = cell.link.4, label=T, plot.param = list(pt.bg=ID))
plotAllSpecimens(gpa.2$coords, mean=T, links = cell.link.2, label=T, plot.param = list(pt.bg=ID))


### plot PCA space by species ####
plotTangentSpace(gpa.4$coords, groups=ID, axis1 = 1, axis2 = 3)
plotTangentSpace(gpa.2$coords, groups=ID, axis1 = 1, axis2 = 2)



### plotting shape change from mean to each species ####
ref.4<-mshape(gpa.4$coords)
lsmeans.4<-arrayspecs((rowsum(two.d.array(gpa.4$coords), ID)/as.vector(table(ID))),18,2) #gp means

plotRefToTarget(lsmeans.4[,,2],lsmeans.4[,,1], links=cell.link.4, mag=2, method="points")  #plots elegans vs briggsae
plotRefToTarget(ref.4,lsmeans.4[,,1], links=cell.link.4, mag=2, method="points")  #plots briggsae vs mean
#plotRefToTarget(ref.4,lsmeans.4[,,1], links=cell.link.4, mag=3, method="TPS")  #uses thin plate spline

ref.2<-mshape(gpa.2$coords)
lsmeans.2<-arrayspecs((rowsum(two.d.array(gpa.2$coords), ID)/as.vector(table(ID))),12,2) #gp means
plotRefToTarget(lsmeans.2[,,2],lsmeans.2[,,1], links=cell.link.2, mag=2, method="points")  #magnified 2X


par(mfrow=c(1,2), mar=c(1,0,1,0) + 0.1)
plotRefToTarget(ref.2,lsmeans.2[,,1], links=cell.link.2, mag=3, method="vector", gridPars = gridPar(pt.size = 0.666))  
mtext(dimnames(lsmeans.2)[[3]][1], side = 3, line =  -5)
plotRefToTarget(ref.2,lsmeans.2[,,2], links=cell.link.2, mag=3, method="vector", gridPars = gridPar(pt.size = 0.666))  
mtext(dimnames(lsmeans.2)[[3]][2], side = 3, line =  -5)
par(mfrow=c(1,1))

#this is procrustes distance linear model that tests for differences in shape based on ID (i.e. species)
procD.lm(gpa.4$coords ~ ID, iter=999, RRPP=T)
procD.lm(gpa.2$coords ~ ID, iter=999, RRPP=T)


#renaming landmark data
Y.4 <- gpa.4$coords
Y.2 <- gpa.2$coords
#setting modules based on landmarks (in this case, P2 vs the rest of the embryo)
land.gps.4 <- rep('a',18); land.gps.4[c(3,14,7,15,10,6)]<-'b'
land.gps.2 <- rep('a',12); land.gps.2[c(1,4)]<-'b'

TREE <- "(elegans:4.2,briggsae:4.2);"
TREE
tree.ex <- read.tree(text = TREE)
str(tree.ex)
tree.ex


### plot movement of cell centroids across species #######





### running modularity test to see if above-defined modules differ from one another in shape ####
MT.4 <- modularity.test(Y.4, land.gps.4, CI=FALSE, iter=199)
MT.2 <- modularity.test(Y.2, land.gps.2, CI=FALSE, iter=199)
MT <- phylo.modularity(A = Y.gpa, partition.gp = land.gps, phy = tree.ex, CI=FALSE, iter=199)
summary(MT) # Test summary

mtext(paste("Anterior/Posterior Cell Border Modularity -- P: ", MT.2$P.value, sep=""), side = 1, line = -5)


#may want to look into using phylo.modularity() 
#find out how to get phylogenetic tree



#############################deprecated code

#populates vector "filelist" with names of all images
imagelist <- as.vector(list.files(pattern="*.jpg"))

#this was originally read in as a .nts file, but apparently geomorph only reads data in as a .tps format (despite the given examples running digitize2d for .nts)
tpslist <- paste(strsplit(filelist, ".jpg"), ".tps", sep="")

datalist <- paste(sub(".jpg","",filelist),".tps",sep="")


for(i in 1:length(imagelist)){
  digitize2d(filelist = imagelist[i], nlandmarks = 6, scale = NULL, tpsfile = tpslist[i])
}

cells <- data.frame()
for (i in 1:length(datalist)){
  cells <- rbind(cells, two.d.array(readland.tps(datalist[i])))
}

#after initial analysis, I found that the following samples had been landmarked incorrectly (counterclockwise vs clockwise)
flipped <- c(20,19,18,16,15,13,10,6,4)
for (i in flipped){
  cells[i, c(5,11)] <- cells[i, c(11,5)] #switches X3 with X6
  cells[i, c(6,12)] <- cells[i, c(12,6)] #switches Y3 with Y6
  cells[i, c(7,9)] <- cells[i, c(9,7)] #switches X4 with X5
  cells[i, c(8,10)] <- cells[i, c(10,8)] #switches Y4 with Y5
}

cell.land <- arrayspecs(cells, 10, 2)


##adds 'sex' to end of data frame
n <-nrow(twod.cell)
IDs <- paste(sub(".jpg","",filelist),"",sep="")[-22]
cc <- strsplit(IDs, '-')
ID <- as.factor(unlist(cc)[2*(1:length(IDs))-1])
cell.frame<-data.frame(twod.cell, ID)