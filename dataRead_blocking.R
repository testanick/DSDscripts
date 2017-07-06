

### Read in and process data ------------------------------------------------
day1.2 <- two.d.array(readland.tps("images/Blocking/3-13/2-cell/day1_2-cell.TPS", specID = "imageID"))
day2.2 <- two.d.array(readland.tps("images/Blocking/3-17/2-cell/day2_2-cell.TPS", specID = "imageID"))
day3.2 <- two.d.array(readland.tps("images/Blocking/3-20/2-cell/day3_2-cell.TPS", specID = "imageID"))
day4.2 <- two.d.array(readland.tps("images/Blocking/3-24/2-cell/day4_2-cell.TPS", specID = "imageID"))

day1.4 <- two.d.array(readland.tps("images/Blocking/3-13/4-cell/day1_4-cell.TPS", specID = "imageID"))
day2.4 <- two.d.array(readland.tps("images/Blocking/3-17/4-cell/day2_4-cell.TPS", specID = "imageID"))
day3.4 <- two.d.array(readland.tps("images/Blocking/3-20/4-cell/day3_4-cell.TPS", specID = "imageID"))
day4.4 <- two.d.array(readland.tps("images/Blocking/3-24/4-cell/day4_4-cell.TPS", specID = "imageID"))

cells.2 <- arrayspecs(rbind(day1.2, day2.2, day3.2, day4.2), 12, 2)
cells.4 <- arrayspecs(rbind(day1.4, day2.4, day3.4, day4.4), 18, 2)

### assign ID as species names based on raw image names -------------------
day <- as.factor(unlist(strsplit(dimnames(cells.2)[[3]], "_|[.]"))[ c(TRUE,FALSE,FALSE) ])
ID <- as.factor(unlist(strsplit(dimnames(cells.2)[[3]], "_|[.]"))[ c(FALSE,TRUE,FALSE) ])


# run Generalized Procrustes Analysis ----------------------------------
gpa.2 <- gpagen(cells.2, ProcD = T, curves = curvepts.2)#, PrinAxes = T)
gpa.4 <- gpagen(cells.4, ProcD = T, curves = curvepts.4)#, PrinAxes = T)




### turn data back to 2D ------------------------------------------------
two.d.2 <- two.d.array(gpa.2$coords); rownames(two.d.2) <- day
two.d.4 <- two.d.array(gpa.4$coords); rownames(two.d.4) <- day




### calculate reference and target shapes ----------------------------
ref.2<-mshape(gpa.2$coords)
lsmeans.2<-arrayspecs((rowsum(two.d.array(gpa.2$coords), day)/as.vector(table(day))),12,2) #gp means

ref.4<-mshape(gpa.4$coords)
lsmeans.4<-arrayspecs((rowsum(two.d.array(gpa.4$coords), day)/as.vector(table(day))),18,2) #gp means

### calculate principle components of shape --------------------------
pcCell.2 <- data.frame(prcomp(two.d.2)$x[,1:20], day)  #20 is the point at which ~100% of cumulative variance is explained
pcCell.4 <- data.frame(prcomp(two.d.4)$x[,1:30], day) #30 is the point at which ~100% of cumulative variance is explained
