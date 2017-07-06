
#######################################################


##### Plotting confusion matrices as pretty heatmaps
myHeatMap <- function( dim, names, data, xLab="Actual", yLab="Predicted" ){
	myCols <- colorRampPalette(c("#ffffff7f", "#ff00007f", "#00ff007f"), space="rgb")
	image( 1:dim, 1:dim, data[,dim:1], col=myCols(100), xaxt='n', yaxt='n', xlab=xLab, ylab=yLab )
	abline(h=0:dim+0.5, col="grey")
	abline(v=0:dim+0.5, col="grey")
	box( lwd=1 )
	text( 1:dim, rep(dim:1, each=dim), sub('^0$', '',data) )
	axis(1, at=1:dim, labels=names, cex.axis=0.8 )
	axis(2, at=1:dim, labels=rev(names), cex.axis=0.8, las=1 )
}




# this function works with the 'pairs' plotting function to add numerical correlation estimates to the plots
panel.cor <- function(x, y, digits=2) {
	r <- cor(x, y)
	par( usr =c(0, 1, 0, 1))
	Cor <- format(c(r, 0.123456789), digits=digits)[1]
 	text(0.5, 0.5, paste("r=", Cor), cex=1.5)
}


# Inverse log function
invlogit <- function( x ) { 1 / (1 + exp(-x) ) }


# geometric mean function
g_mean <- function( X ) {
    log_data <- log( X )
    gm <- exp( mean( log_data[ is.finite( log_data )] ))
    return( gm )
}

# mode (central tendency, not data-storage type) function

Mode <- function( x ) {
  ux <- unique( x )
  ux[ which.max( tabulate( match( x, ux))) ]
}


# 'apply'-able HPDinterval function
hpd <- function( X ) {
  HPDinterval( as.mcmc( X ), probs=c(0.025, 0.975) )[1:2]
}


# subset data by keeping only named columns
# cols <- c( "factorname", "variablename", etc... )
keepCols <- function(data, cols) {
	return( data[ , names(data) %in% cols] )
	}


# I got bored never remembering which 'pch' value I wanted - this reminds me...
show_pch <- function( ){
		X <- rep(1:5, 5)
		Y <- c( rep(5,5), rep(4,5), rep(3,5), rep(2,5), rep(1,5))
		par( mar=c(1, 1, 3, 1))
		plot( X, Y, type='n', xlab="", ylab="", main="PCH 21-25 are fillable", yaxt='n', xaxt='n', ylim=c(0,5) )
	for (i in 1:25){
			points( X[i], Y[i], pch=i, cex=3 )
			text( X[i], Y[i]-0.5, as.character(i) )
		}	}


# The obligatory homebrew standard error function
	stder <- function(x) sd(x)/sqrt(length(x))


# ...and the similarly obligatory coefficient of variation function
	CV <- function( x ) 100 * sd( x ) / mean( x )


# this function combines two vectors e.g. x1, x2..xn and y1, y2..yn to produce x1, y1, x2, y2...xn, yn
interleave <- function(v1, v2) {
	ord1 <- 2 * (1:length(v1))-1
	ord2 <- 2 * (1:length(v2))
	c(v1, v2)[order(c(ord1, ord2))]
	}

comment(interleave) <- "this function takes 2 vectors (a1, a2..an), (b1, b2..bn) and returns (a1, b1, a2, b2..an, bn)"




#ID's function to count minor alleles - if there are 2 or fewer instances then the locus is a statistical no-go...
minor.allele.function <- function(x) length(x[x==2])

comment(minor.allele.function)<-"counts the instances of the minor allele in 0/2 coded SNP data"


#MANOVA function for use in an apply-like statement - returns test stats only - will return NA if there are <3 instances of a minor case
manova.function3 <- function( test.SNP.data ) {
    D4ta <-  merge( test.SNP.data , line.PW.means.M )
    prob	<- minor.allele.function ( D4ta$SNPrecod )
    if (prob <= 2 ){
        return (c(Wilks=NA, p.value=NA))    }
    model <- manova( fmla, data=D4ta  )
    Wilks <- summary(model, test="Wilks")$stats[1,2]
    p     <- summary(model, test="Wilks")$stats[1,6]
    return(c(Wilks=Wilks, p.value=p))     }

comment(manova.function3)<-"performs a manova and returns test stats - will return NA if there are <=2 instances of the minor factor level"


ShapeScore <- function(Beta, Y) {
	## This computes the "shape score" of Drake and Klingenberg 2008
	## Beta is the coefficients for each procrustes coordinate from a multivariate multiple regression
	## Y is the raw procrustes coordinates
	s = Y %*% t(Beta) %*% ((Beta %*% t(Beta))^-0.5)
	return(shapeScore=s)
}


# univariate Rsquared and partial Rsquared calculators

Rsq <- function( model ){
	fitted.variance <- var(model$fitted)
	total.variance	<- var(model$fitted) + var(model$resid)
	fitted.variance / total.variance
}


PRsq <- function( model ){
	residual.variance <- var(model$resid)
	variables <- attr(terms(model), "term.labels")
		model.length <- length(variables)
		variable.name <- rep(NA, model.length )
		partial.Rsq <- rep(NA, model.length )
		univariate.Model.Rsq <- rep(NA, model.length )

	for (i in 1:model.length){
		variable.name[i] <- variables[i]
		drop <- parse( text=variables[i] )
		new.formula <- as.formula( paste( ".~.-", variables[i], sep=""))
		new.model <- update(model, new.formula )
		partial.Rsq[i] <- (var(new.model$resid) - residual.variance)/ var(new.model$resid)

		new.formula.univariate <- as.formula( paste( ".~", variables[i], sep=""))
		univariate.model <- update(model, new.formula.univariate)
		univariate.Model.Rsq[i] <- summary(univariate.model)$r.sq
		}

	R2 <- Rsq( model )
	adj.R2 <- summary(model)$adj.r

	partials <- data.frame(partial.Rsq, univariate.Model.Rsq )
	row.names(partials) <- variable.name

	list(FullModelRsquared=R2, FullModelAdjustedR2 = adj.R2, partials=partials	)
}



# multivariate Rsquared and partial Rsquared calculators

shapeRsq <- function( model ){
	fitted.variance <- sum(diag(var(model$fitted)))
	total.variance	<- sum(diag(var(model$fitted + model$resid)))
	fitted.variance / total.variance
}

comment(shapeRsq) <- "this function takes a miultivariate model object and returns the fitted variance / total variance: eqivalent to an Rsquared"


shapePRsq <- function( model ){
# Based on the derivation from page 269 of Kutner et. al. (Applied linear statistical models edition 5.)
	residual.variance <- var(model$resid)
	variables <- attr(terms(model), "term.labels")
		model.length <- length(variables)
		variable.name <- rep(NA, model.length )
		partial.Rsq <- rep(NA, model.length )

	for (i in 1:model.length){
		variable.name[i] <- variables[i]
		drop <- parse( text=variables[i] )
		new.formula <- as.formula( paste( ".~.-", variables[i], sep=""))
		new.model <- update(model, new.formula )
		partial.Rsq[i] <- (sum ( diag( var(new.model$resid))) - sum( diag( residual.variance)) ) / sum( diag( var(new.model$resid)))
		}
	R2 <- shapeRsq( model )
	list(Rsquared=R2, partials=data.frame( cbind( variable.name, partial.Rsq ))	)
}



####################################

# calculate tangent approximaton for tangent approximates Procrustes Distance (Euclidean Distance)
# This is just the magnitude of the vector!

PD <- function( x ) {
				sqrt( t( x ) %*% x )}
comment(PD) <- c("This just computes the Euclidean Distance (norm) for a vector")


###################################


# This function will compute the angle between two vectors.	BEWARE if the vectors can be of arbitrary sign, use the function below this one.

ang.vec <- function(vec1, vec2){
	vec.cor <- (t(vec1) %*% vec2)/(PD(vec1)*PD(vec2))
	vec.angle <- acos(vec.cor)*(180/pi)
	return(c(vector.cor=vec.cor, vec.angle=vec.angle))}
comment(ang.vec) <- c(" This computes both the vector correlation, and angle, between two vectors.", " to compare to the Pearson correlation coefficient make sure to center and standardize vectors", "DO NOT USE THIS IF VECTORS CAN BE OF ARBITRARY SIGN")



####### When the vectors can be of arbitrary sign, use this which computes the magnitude of the vector correlation, and then computes the angle.

ang.vec.abs <- function(vec1, vec2){
	vec.cor <- abs((t(vec1) %*% vec2)/(PD(vec1)*PD(vec2)))
	vec.angle <- acos(vec.cor)*(180/pi)
	return(c(vector.cor=vec.cor, vec.angle=vec.angle))}
comment(ang.vec) <- c(" This computes both the vector correlation, and angle, between two vectors.", " to compare to the Pearson correlation coefficient make sure to center and standardize vectors", "set it up to compute the absolute values of the vector correlation")



#################################

# PLS function from Claude 2010, with modification - updated 4-12-12 WP
PLS <- function( M1, M2 ) {
	p1 <- dim(M1)[2]
	# we allow p2 to be a simple vector, instead of a matrix, thus:-
	if ( is.vector(M1) == T ){
		stop( paste( "'M1' must be specified as a matrix" ) )
		 }
	if ( is.vector(M2) == T ) { p2 <- 1 }
		else { p2 = dim(M2)[2] }

	n <- dim(M1)[1]
	sM12 <- svd(var(cbind(M1,M2))[1:p1, (p1+1):(p1+p2)])
	vM12 <- var(cbind(M1,M2))[1:p1, (p1+1):(p1+p2)]
	vM21 <- var(cbind(M1,M2))[(p1+1):(p1+p2), 1:p1]
	v11 <- var(M1)
	v22 <- var(M2)
	D <- sM12$d
	F1 <- sM12$u
	F2 <- sM12$v
	Rv <- sum(diag(vM12%*%vM21))/sqrt(sum(diag(v11%*%v11))*
	sum(diag(v22%*%v22)))
	return( list(Rv_coef=Rv, Singular_Values=D, block1_PLS_vecs=F1, block2_PLS_vecs=F2) )
}

comment(PLS) <- "M1 should be a matrix of data. M2 may be a matrix or vector of data"


############################################

# Klingenberg's "shape score"
# now let us compute Klingenberg's "shape score" as per Drake and Klingenberg

ShapeScore <- function(Beta, Y) {
	## This computes the "shape score" of Drake and Klingenberg 2008
	## Beta is the coefficients for each procrustes coordinate from a multivariate multiple regression.
	## i.e. shape.model.1 <- lm(cbind(Y1, Y2) ~ X1 + X2, data=data.set )
	## Beta ...   beta.model <- t(coef(shape.model.1)[k,]) # where k is the covariate of focus
	##
	## Y is the raw procrustes coordinates (Y1, Y2, ...)
	# Since Klingenberg uses the generalized inverses, we have tested this for the 2p-4 PCs and have gotten identical results.
	s = Y %*% t(Beta) %*% ((Beta %*% t(Beta))^-0.5)
	return(shapeScore=s)
}

comment(ShapeScore) <- c("This computes the 'shape score' of Drake and Klingenberg 2008",
  "Beta is the coefficients for each procrustes coordinate from a multivariate multiple regression. ",
  "i.e. shape.model.1 <- lm(cbind(Y1, Y2) ~ X1 + X2, data=data.set )",
  "Beta ...   beta.model <- t(coef(shape.model.1)[k,]) # where k is the covariate of focus",
  "Y is the raw procrustes coordinates (Y1, Y2, ...)",
  "Klingenberg uses the generalized inverses, we have tested this for the 2p-4 PCs and have gotten identical results.")


#################################

# Back-Transformation from a vector of coefficients calculated from Principal Components to a vector of coefficients applicable to the Proc. coordinates from which the PC's were calculated
# Main purpose - allow easy production of wingplots from analysis output based on principal components that must be necessarily of lower dimensionality than the coordinates

BackTrans <- function( PC_coefs, PC_object, real_dimensions=58 ){
	if ( class(PC_object) == "prcomp" ){
		rotation_of_real_dimensions <- PC_object$rotation[, 1:real_dimensions ]
		} else {
		rotation_of_real_dimensions <- PC_object[, 1:real_dimensions ]
		}
	coefs_in_LM_space <- rotation_of_real_dimensions %*% PC_coefs
	coefs_in_LM_space
}

comment(BackTrans) <- c("Backtransform a reduced-dimensionality vector of PC coefs to a vector LM coefs. Default dimensionality of PC's is 58")


##################################
KirkpatrickDimensionality <- function(x){
	## This implements Kirkpatrick's (2009) simple measures of evolvability and dimensionality
	## The input is a Genetic covariance matrix for mean standardized traits.
	EigenValues <- svd(x)$d
	EffDim <- sum(EigenValues/EigenValues[1])
	MaxEvol <- sqrt(EigenValues[1])
	TotalGenetVar <- sum(EigenValues)
    return( c(Dimensionality = EffDim, MaximumEvolvability = MaxEvol, TotalGenetVar = TotalGenetVar  ))
}
comment(KirkpatrickDimensionality) <- c("This implements Kirkpatrick's (2009) simple measures of evolvability and dimensionality", "The input is a Genetic covariance matrix for mean standardized traits.")




######################################
# Function to computer the CI's for vector correlations and angles.
###### This still needs to be generalized for either MCMC or bootstrapping. Then functionalized.

##### Currently Using the posterior distributions to make the comparisons

# Grab the resampled data or posterior MCMC simulations
# data.1 <- model.legs.2.wt.18$VCV #SET 1
# data.2 <- model.legs.3.wt.18$VCV #SET2

# # Make sure the resampling/ posterior sims are the same length (number of iterations)
# N1 <- nrow(data.1)
# N2 <- nrow(data.2)
# N1==N2 # check to make sure the data sets are the same size. Make this into a function, and if it fails, either give an error or just use min(N1,N2) for length

# # Set up empty matrix to store the results for the vector correlations and angles for each of three eigenvectors
# # Make a variable in the function for the number of vectors to compare (and change ncol)
# angles.vectors <- matrix(NA, ncol=6, nrow=N1)

# # just setting up names. This can easily be extended using paste or interaction.
# colnames(angles.vectors) <- c("VecCor1", "Angle1","VecCor2", "Angle2", "VecCor3", "Angle3")

# for (i in 1:N1) {
	# # First extract the sample of the covariance matrix from the posterior simulation
	# cov.1 <- matrix(data.1[i,1:9], nrow=3, ncol=3)
	# cov.2 <- matrix(data.2[i,1:9], nrow=3, ncol=3)
    # # compute the vector correlations and angles. I could use svd(), also needs to be generalized to more than 3 vectors.
    # angles.vectors[i,] <- c(  ang.vec.abs(eigen(cov.1)$vectors[,1], eigen(cov.2)$vectors[,1]),
                              # ang.vec.abs(eigen(cov.1)$vectors[,2], eigen(cov.2)$vectors[,2]),
                              # ang.vec.abs(eigen(cov.1)$vectors[,3], eigen(cov.2)$vectors[,3]) )
# }



########################### WHAT STILL NEEDS TO BE IMPLEMENTED
# PROCRUSTES ANOVA
# VARIANCE ACCOUNTED FOR VIA PROCRUSTES DISTANCE
