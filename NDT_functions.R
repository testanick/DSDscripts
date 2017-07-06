
# calculate tangent approximaton for tangent approximates Procrustes Distance (Euclidean Distance)
# This is just the magnitude of the vector!
PD <- function( x ) {
  sqrt( t( x ) %*% x )}
comment(PD) <- c("This just computes the Euclidean Distance (norm) for a vector")



####### When the vectors can be of arbitrary sign, use this which computes the magnitude of the vector correlation, and then computes the angle.
ang.vec.abs <- function(vec1, vec2){
  vec.cor <- abs((t(vec1) %*% vec2)/(PD(vec1)*PD(vec2)))
  vec.angle <- acos(vec.cor)*(180/pi)
  return(c(vector.cor=vec.cor, vec.angle=vec.angle))}
comment(ang.vec.abs) <- c(" This computes both the vector correlation, and angle, between two vectors.", " to compare to the Pearson correlation coefficient make sure to center and standardize vectors", "set it up to compute the absolute values of the vector correlation")



getShape <- function(land, gpa, id)
{
  temp2 <- data.frame(land, gpa$Csize)
  temp3 <- temp2[sample(nrow(temp2), replace=T),]
  temp <- subset(temp3, ID==id)
  if (nrow(temp)>=2)
  {  
    len <- ncol(temp)-2
    shape <- as.matrix(temp[,c(1:len)])
    size <- temp$gpa.Csize
    allometry.model <- lm(shape ~ size)
    Beta <- allometry.model$coefficients[2,] 
    shape <- ShapeScore(t(Beta), shape)
    return(mean(shape))
  } else{return(NA)}
}

getShapeByDay <- function(land, gpa, days)
{
  temp2 <- data.frame(land, gpa$Csize)
  temp3 <- temp2[sample(nrow(temp2), replace=T),]
  temp <- subset(temp3, day==days)
  if (nrow(temp)>=2)
  {  
    len <- ncol(temp)-2
    shape <- as.matrix(temp[,c(1:len)])
    size <- temp$gpa.Csize
    allometry.model <- lm(shape ~ size)
    Beta <- allometry.model$coefficients[2,] 
    shape <- ShapeScore(t(Beta), shape)
    return(mean(shape))
  } else{return(NA)}
}

# sexShapes <-data.frame()
# for (t in unique(wings2$treatment))
# {
#   male.shape <- getShape(t,"m")
#   female.shape <- getShape(t,"f")
#   temp.line <- cbind(t, male.shape, female.shape)
#   sexShapes <- rbind(sexShapes, temp.line)
# }

