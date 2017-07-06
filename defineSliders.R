#defines sliders based on most recent 18-lm configuration of 4-cell embryos
#

curvepts.2<-matrix(c(
  2,5,6,
  5,6,1,
  1,7,8,
  7,8,3,
  3,9,10,
  9,10,4,
  4,11,12,
  11,12,2
), ncol=3, byrow=T)
colnames(curvepts.2) <- c("before", "slide", "after")

curvepts.4<-matrix(c(
  4,11,1,
  1,12,2,
  2,13,3,
  3,14,7,
  7,15,10,
  10,16,9,
  9,17,8,
  8,18,4,
  12,2,13,
  14,7,15,
  16,9,17,
  18,4,11
         ), ncol=3, byrow=T)
colnames(curvepts.4) <- c("before", "slide", "after")


cell.link.2 <- matrix(c(2,5,
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

cell.link.4 <- matrix(c(4,11,
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
                        5,8,
                        5,6,
                        6,3,
                        6,10
                        ),ncol=2,byrow=T)
