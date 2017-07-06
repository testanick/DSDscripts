

### Read in and process data ------------------------------------------------
cells <- two.d.array(readland.tps("images/40x/40x.tps", specID = "imageID"))
cells <- two.d.array(readland.tps("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/DSD Shape Data/images/40x/40x.tps", specID = "imageID"))

cell <- arrayspecs(cells, 18, 2)

# devTemp <- as.factor(c(rep("20C", length(cells.2[1,1,])), 
#                        rep("30C", length(cells.30.2[1,1,]))))


### assign ID as species names based on raw image names -------------------
ID <- as.factor(unlist(strsplit(dimnames(cell)[[3]], "_"))[ c(TRUE,FALSE,FALSE) ])
devTemp <- as.factor(unlist(strsplit(dimnames(cell)[[3]], "_"))[ c(FALSE,TRUE,FALSE) ])
IDnum_1 <- unlist(strsplit(dimnames(cell)[[3]], "_"))[ c(FALSE,FALSE,TRUE) ]
IDnum <-  as.factor(unlist(strsplit(unlist(strsplit(IDnum_1, "\\."))[c(TRUE,FALSE)], "-"))[c(F,F,F,T)])
date <- as.factor(unlist(strsplit(IDnum_1, "-"))[c(F,T,F,F)])

# run Generalized Procrustes Analysis ----------------------------------
gpa.cell <- gpagen(cell, ProcD = T, curves = curvepts.4) #, PrinAxes = T)


### turn data back to 2D ------------------------------------------------
two.d.cell <- two.d.array(gpa.cell$coords); rownames(two.d.cell) <- ID

### calculate reference and target shapes ----------------------------
ref <- mshape(gpa.cell$coords)
lsmeans <- arrayspecs((rowsum(two.d.array(gpa.cell$coords), paste(ID,devTemp,sep="_"))/as.vector(table(paste(ID,devTemp,sep="_")))),18,2) #gp means

### calculate principle components of shape --------------------------
pc.cell <- data.frame(prcomp(two.d.cell)$x[,1:26], ID, devTemp) #26 is the point at which >99.9% of cumulative variance is explained
