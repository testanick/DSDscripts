

### Read in and process data ------------------------------------------------
N2.30.2 <- two.d.array(readland.tps("preliminary\ images/30C/2-cell/N2.tps", specID = "imageID"))
AF16.30.2 <- two.d.array(readland.tps("preliminary\ images/30C/2-cell/AF16.tps", specID = "imageID"))
AF16.2 <- two.d.array(readland.tps("images/20C/2-cell/AF16/AF16.tps", specID = "imageID"))
N2.2 <- two.d.array(readland.tps("images/20C/2-cell/N2/N2.tps", specID = "imageID"))


N2.30.4 <- two.d.array(readland.tps("preliminary\ images/30C/4-cell/N2.tps", specID = "imageID"))
AF16.30.4 <- two.d.array(readland.tps("preliminary\ images/30C/4-cell/AF16.tps", specID = "imageID"))
AF16.4 <- two.d.array(readland.tps("images/20C/4-cell/AF16/AF16.tps", specID = "imageID"))
N2.4 <- two.d.array(readland.tps("images/20C/4-cell/N2/N2.tps", specID = "imageID"))

cells.2 <- arrayspecs(rbind(AF16.2, N2.2), 12, 2)
cells.4 <- arrayspecs(rbind(AF16.4, N2.4), 18, 2)

cells.30.2 <- arrayspecs(rbind(AF16.30.2, N2.30.2), 12, 2)
cells.30.4 <- arrayspecs(rbind(AF16.30.4, N2.30.4), 18, 2)

 devTemp <- as.factor(c(rep("20C", length(cells.2[1,1,])), 
                        rep("30C", length(cells.30.2[1,1,]))))



cells.2 <- arrayspecs(rbind(AF16.2, N2.2, AF16.30.2, N2.30.2), 12, 2)
cells.4 <- arrayspecs(rbind(AF16.4, N2.4, AF16.30.4, N2.30.4), 18, 2)

# cells.2 <- arrayspecs(rbind(AF16.2, N2.2), 12, 2)
# cells.4 <- arrayspecs(rbind(AF16.4, N2.4), 18, 2)

### assign ID as species names based on raw image names -------------------
ID <- as.factor(unlist(strsplit(dimnames(cells.2)[[3]], "_"))[ c(TRUE,FALSE) ])

# run Generalized Procrustes Analysis ----------------------------------
gpa.2 <- gpagen(cells.2, ProcD = T, curves = curvepts.2)#, PrinAxes = T)
gpa.4 <- gpagen(cells.4, ProcD = T, curves = curvepts.4)#, PrinAxes = T)




### turn data back to 2D ------------------------------------------------
two.d.2 <- two.d.array(gpa.2$coords); rownames(two.d.2) <- ID
two.d.4 <- two.d.array(gpa.4$coords); rownames(two.d.4) <- ID




### calculate reference and target shapes ----------------------------
ref.2<-mshape(gpa.2$coords)
lsmeans.2<-arrayspecs((rowsum(two.d.array(gpa.2$coords), ID)/as.vector(table(ID))),12,2) #gp means

ref.4<-mshape(gpa.4$coords)
lsmeans.4<-arrayspecs((rowsum(two.d.array(gpa.4$coords), ID)/as.vector(table(ID))),18,2) #gp means

### calculate principle components of shape --------------------------
pcCell.2 <- data.frame(prcomp(two.d.2)$x[,1:20], ID)  #20 is the point at which ~100% of cumulative variance is explained
pcCell.4 <- data.frame(prcomp(two.d.4)$x[,1:30], ID) #30 is the point at which ~100% of cumulative variance is explained
