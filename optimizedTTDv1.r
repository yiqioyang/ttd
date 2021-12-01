# The algorithm is written in R. The lines starting with or followed by "#" are comments for functions, command lines, or delimiters of separate functions.
###### Load packages 
###### Author: Qingyuan Yang and Marcus Bursik
###### License: GPL. Use at your own risk.
###### Requires: to run the examples, the appropriate CSV files are needed.
library(gstat)
library(sp)
library(rgdal)
library(RandomFields)
library(maptools)
library(raster)
library(deldir)
library(ggplot2)
library(segmented)
library(poweRlaw)
library(segmented)
library(geoR)
library(combinat)


###### Workflow #########################################################################################################################################
######(1) Calculate the distances, downwind distances of sample sites with respect to the input       <------ function "prep"				######
######    source vent location and wind direction.	 															######
######																								######
######(2) Apply segmented/simple linear regression to the datasets: the log-thickness is the linear   <------ function "sample.reg"			######
######    combination of distance and downwind distance plus a constant.													######
######																								######
######(3) Record the fitted values and residuals from the fitting in step (2)                         <------ implemented by command lines		######
######																								######
######(4) Form a grid for mapping of the trend and the final result		            		          <------ function "form.grid"			######
######																								######
######(5) Form the trend model 											                              <------ function "prep2" and "trendmodel"	######
######																								######
######(6) Form and fit the variogram function of the residuals recorded in step (3)	                  <------ function "var.read"			######
######																								######
######(7) Implemnt ordinary kriging using the residuals in step (3) and the variogram models          <------ function "prediction"			######
######    from step (6) to form the local thickness variation. We take the sum of 			                    							######
######    the trend model in step (5) and the local variation to form the final thickness distribution								######
######																								######
###### End of workflow ##################################################################################################################################



###### Functions ########################################################################################################################################
prep <- function(sv, td, ang){	
											
### See step (1) in the workflow above
### Input variables: 									  ### Data type			  ### Column naes
###		sv: the source vent location							dataframe				"x" and "y"
###		td: the tephra thickness dataset						dataframe				"x", "y", "rz" (real thickness), "z"(log-transformed thickness)
###		ang: the inferred wind direction (from north clockwise)		a single value			     	   -

### Intermediate variables:
###		pts: the coordinates of sample sites
###		nr: the number of observations
###		sv.mtx: the matrix prepared for the computation of distance and downwind distance
###		disp.mtx: the unit vector pointing along the wind direction for the computation of downwind distance
###		dist, dd, cd: computed distances, downwind distances, and crosswind distances from the sample sites to
###		the source vent with the inferred wind direction

### Output:
###		The output is a dataframe that contains the original thickness dataset and the calculated distances, downwind
###		distances, and the crosswind distances of the sample sites.


### Method:
###		This function uses a simple dot product to calculate the downwind distance

	pts <- as.matrix(td[,c(1,2)])									### Extract the coordinates of sample sites
	sv <- as.matrix(sv)
	nr <- nrow(pts)	


	sv.mtx <- matrix(ncol=2, nrow=nr)								### Repeated rows of the coordinates of the source vent
	sv.mtx[,1] <- rep(sv[1], nr)
	sv.mtx[,2] <- rep(sv[2], nr)


	disp.mtx <- c(sin(ang*pi/180),cos(ang*pi/180))						### Form the unit vector pointing to the wind direction

	
	pts.mtx <- pts-sv.mtx										### Calculate the distances, downwind distances, and crosswind distances
	pts.sq <- (pts.mtx)^2										### of the sample sites
	dist<- (pts.sq[,1]+pts.sq[,2])^0.5
	dd <- pts.mtx %*% disp.mtx
	cd <-(dist^2-dd^2)^0.5
      

	output <- as.data.frame(cbind(pts, dist, dd, cd, td[, c(3,4)]))			### Reorganize the data for output
	colnames(output) <- c("x", "y", "dist", "dd", "cd", "rz", "z")
		

	return(output)
}

###-----------------

prep2 <- function(sv, td, ang){

### This function is written for the computation of the trend model, which requires distances and 
### downwind distances of all the cells within the region of interest. Therefore this function works
### exactly the same way as function "prep", except that the trend thicknesses of these cells are to
### be computed. This function is called in the function "trendmodel" for further calculation of trend thickness.


### Input variables: 									  ### Data type			  ### Column names
###		sv: the source vent location			     				dataframe				"x" and "y"
###		td: the coordinates of the cells of the grid				dataframe or matrix		does not matter, as they are transfered into matrix in this function
###		ang: the inferred wind direction (from north clockwise)		a single value			     -


### Output variables:
### 		The output is a dataframe that records the coordinates, distances, downwind distances, and crosswind distances
###		of the cells of the grid.

	pts <- as.matrix(td[,c(1,2)])
	sv <- as.matrix(sv)
	nr <- nrow(pts)	


	sv.mtx <- matrix(ncol=2, nrow=nr)
	sv.mtx[,1] <- rep(sv[1], nr)
	sv.mtx[,2] <- rep(sv[2], nr)
	disp.mtx <- c(sin(ang*pi/180),cos(ang*pi/180))

	
	pts.mtx <- pts-sv.mtx
	pts.sq <- (pts.mtx)^2
	dist<- (pts.sq[,1]+pts.sq[,2])^0.5
	dd <- pts.mtx %*% disp.mtx
	cd <-(dist^2-dd^2)^0.5
      

	output <- as.data.frame(cbind(pts, dist, dd, cd))
	colnames(output) <- c("x", "y", "dist", "dd", "cd")
		

	return(output)
}

###-----------------

sample.reg <- function(td, segdist, segdd){

### See step (2) of the workflow above
### This function provides two alternative fitting schemes for the fitting of thickness with distance and downwind distance:
### 	(1) Use linear regression with the "segmented" fitting function to simulate
###	    the different decay rates of tephra thickness within and outside the plume corner.

###	(2) Use simple linear regression for the fitting.

###	This method will tend to use the segmented fitting option as it is better reconciled with the physics of tephra settling.
###   If error occurs, meaning that the data are inappropriate for option(1), a "tryCatch" function is set up to switch to 
###	option(2).

###	Input variables:
###		td: the output from function "prep"
###		segdist: the inferred breakpoint for distance, this can be estimated by plotting the thickness against distance of the observations.
###		segdd: the inferred breakpoint for downwind distance.

### 	Intermediate functions:
###		segfit: segmented fitting function
###		simplefit: simple linear regression

### 	Output variables:
###		An object that can be regarded as the fitted model. It contains the fitted values, residuals and the
###		coefficients of the fitting. It can be used to predict the trend thickness for the grid and also store the residuals 
###		for ordinary kriging.

###	*** NOTES ***
###		The segemented fitting runs iterative processes. Sometimes even with the same "segdist" and "segdd" inputs, it still
###		comes up with different breakpoints and consequently slightly different coefficients. Most of the time, such variation
###		in the models is negligible, but we still suggest users to be careful when dealing with the segmented fitting.
###		The goodness of fit can be examined by printing or plotting the residuals, which are stored in the output.
### 		If option (1) is being used, the output also contains the values of the breakpoints for the two variables.
		
	segfit <-function(segdist, segdd){															### Segmented fitting function
		reg.mod <- lm(z~ dist + dd,data=td)
			seg.lmmodel <- segmented(reg.mod, seg.Z = ~dist + dd, psi = list(dist = c(segdist), dd = c(segdd)))
      	return(seg.lmmodel)
	}


	simplefit <-function(){																	### Simple linear regression function
		reg.mod <- lm(z ~ dist + dd ,data=td)
      	return(reg.mod)
	}

	
	result = tryCatch(segfit(segdist, segdd),error=function(e){simplefit()})									### TryCatch structure to automatically switch to simple linear regression
	

  return(result)																			### when error occurs
  }

###-----------------

form.grid <- function(td, cellsize, displacement.coef, expanding.coef){
### See step (4) in the workflow above
### This function first finds the bounding box of the sample sites. It can be used as the extent of the grid.
### However, sometimes one needs a larger extent. To do that, "displacement.coef" and "expanding.coef" are set.
### "displacement.coef" determines if the lowerleft corner of the bounding box should be further shifted to 
### the lower left. "expanding.coef" determines if the length and width of the bouding box should be expanded.
### These two values are all ratios. For example, for a 500*500 squared bounding box, if the 'displacement.coef"
### is set to be 1, then the lowerleft corner (if it was 0,0), will be moved to -500, -500. If the "expanding.coef"
### is set to be 2, then the size of the new bounding box will be 1000*1000.

### Input variables:
### td: the output from function "prep"
### cellsize: cellsize (a single value, as they are squared cells)
### displacement.coef: explained above (a single value)
### expanding.coef: explained above (a single value)


### Output variables:
### An empty grid


  td.sp <- SpatialPoints(cbind(td$x, td$y))
  bbox <- bbox(td.sp)
  cs <- c(cellsize, cellsize)
  cc <- as.numeric(ceiling(bbox[,1] + (cs / 2) - diff(t(bbox)) * displacement.coef))
  cd <- ceiling(diff(t(bbox)) / cs) * expanding.coef
  td.temp.grid <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd)
  td.grid <- SpatialGrid(td.temp.grid)
  return(td.grid)
} 

###-----------------

trendmodel <- function(td, sv, grid, ang, ftng.rslt){
### See step (5) in the workflow above
### This function calculates the distances, downwind distances of the cells of the grid, and makes the
### prediction of trend thickness using the model from function "sample.reg".

### Input variables:
###		td: the output from function "prep"
###		sv: the source vent coordinate
###		grid: the output from function "form.grid"
###		ang: the inferred wind direction (from north clockwise)
###		ftng.rslt: the output from function "sample.reg"


### Output variables:
### A dataframe that contains the coordinates, distances, downwind and crosswind distances, trend thicknesses and 
### trend thicknesses in logarithm form for all the cells within the grid.

	td.grid <- SpatialPoints(as(grid, "SpatialPixels"))
	td.grid.coord <- coordinates(td.grid)


		grided<- prep2(sv, td.grid.coord, ang)

	
		gridz <- predict(ftng.rslt, data.frame(cbind(dist = grided$dist, dd = grided$dd)))
		gridrz <- 10^gridz

		
	output <- as.data.frame(cbind(grided, gridrz, gridz))
	colnames(output) <- c("x", "y", "dist", "dd", "cd", "zr", "z")

	
	return(output)
}

###-----------------

var.read <- function(td, width.val, cutoff.val, var.par, var.model){
### This function reads the locations of the sample sites and the residuals from the fitting in function "sample.reg" to 
### form and fit a variogram model. This function serves as a preparation for ordinary kriging.

### Input variables:
###		td: the output from function "prep"
###		width.val & cutoff.val: the parameters to calculate the experimental variogram. "width.val" denotes the
###		averaged lag for the computation, and "cutoff.val" is the maximum lag that is taken into consideration for the 
###		calculation.
###		var.par: the parameters required for fitting variogram models. It includes three values, namely sill, range, and nugget.
###		Detailed explanations of these parameters can be found in any geostatistics reference.
###		var.model: a text variable which describes the shape of different variogram models, type the command line "show.vgms()" to view a list
###		of variogram models.

### Output variables:
###		The output is a list, the first object of the list is the variogram model, which is the main output.
###		The second object is the plot of the experimental variogram and the fitted variogram curve.
###		The third  object is the values of the experimental variogram

  td.sp <- SpatialPointsDataFrame(cbind(td$x, td$y), td)
  td.var <- variogram(res ~ 1, td.sp, width=width.val, cutoff = cutoff.val)



  td.vgm <- vgm(psill=var.par[1], model=var.model, range=var.par[2], nugget=var.par[3])
  td.fit <- fit.variogram(td.var, td.vgm, fit.sills = c(TRUE, TRUE)) 
  return(list(td.fit,plot(td.var,td.fit),td.var))
}

###-----------------

prediction <- function(td, max, min, maxdist, grid, trend, variog){
### This function first implements ordinary kriging on the residuals to construct the local variation map, and then
### takes the sum of trend and local variation to present the final result.

### Input variables:
###		td: the output from function "prep"
###		max:
###		min:
###		maxdist:
###			The "max" and "min" variables describe the maximum and minimum number of observations taken into account 
###			for the kriging estimation at a single cell. "maxdist" is used to define a maximum radius beyond which the 
###			sample sites are assumed non-related to a particular cell.
###		grid: the output from function "form.grid"
###		trend: the output from function "trendmodel"
###		variog: the first object of the output from function "var.read"

### Output variables:
###		The output is a dataframe that contains the columns of the output of function "trendmodel". In addition to that,
###		it also has the predicted thickness and its logarithm, and the variance of kriging estimation.
###		Here the variance only represents the variance from kriging the log-residuals of the observations and thus cannot
###		be used to characterize the total variance of this method.


	td.sp <- SpatialPointsDataFrame(cbind(td$x, td$y), td)  
	td.temp.k <- krige(res ~ 1, td.sp, nmax = max, nmin = min, maxdist = maxdist, grid, variog)
	k.result <- td.temp.k$var1.pred
	k.result[is.na(k.result)] <- 0
	k.variance <- td.temp.k$var1.var
	k.rz <- 10^(k.result + trend$z)


	output <- cbind(trend, k.result, k.variance, k.rz)


return(output)
}


###### End of functions ####################################################################################################



###### Demo ################################################################################################################

setwd("C:\\Users\\working directory")										### Set working directory, "\" has to be replaced by "\\".  Change to "/" for *nix, mac.
td <- read.table("thickness dataset.csv", header=TRUE, sep=",")						### Read the thickness dataset
sv <- read.table("source vent.csv", header=TRUE, sep=",")							### Read source vent location

td <- prep(sv, td, 20)													### Calculate the distances, downwind distances, and crosswind distances of the sample sites,
																### here the wind direction is 20 degree from north clockwise.

model <- sample.reg(td, 3000,3000)											### Form either the segmented or the simple linear regression model of thickness against
																### distance and downwind distance. The breakpoints for the two variables are set as
																### 3000 and 3000 meters.

td$fit <- model$fitted.values												### Expand "td" by storing the fitted values and
td$res <- model$residuals												### residuals from the fitting in the last command line

grid <- form.grid(td, 500, 0.1, 1.2)											### Form the grid.

trend <- trendmodel(td, sv, grid, 20, model)									### Form the trend. Note the wind direction has to be set again for
																### the calculation of grid cells, and the two wind directions have to be the same.
	
var <- var.read(td, 300,10000, c(0.04,5000,0.01), "Cir")							### Form the variogram model. The parameters are explained in the 
																### description of function "var.read"

result <- prediction(td, 50,20, 5000, grid,trend, var[[1]])							### Calculate the local variation, and take its sum and trend thickness for final result.


result.sp <- SpatialPixelsDataFrame(cbind(trend$x, trend$y), result)					### Turn the final result into "SpatialPixelsDataFrame" for visualization.
spplot(result.sp["k.rz"])	or	spplot(result.sp["zr"])							### Plot the final result or the trend model.

levels <- c(0,250,500,1000,2000,3000,5000,10000,20000)								### Set up the levels for plotting the isopach map.
isomap <- ContourLines2SLDF(contourLines(as.image.SpatialGridDataFrame(result.sp["k.rz"]), levels = levels))
																### Form the isopach map

plot(isomap)															### Plot the isopach map
																### Note: "isomap" is a "SpatialLinesDataFrame" and can be used in GIS programs for further mapping purposes.
																### Geo-referencing can be done in R or GIS programs.

																### ***Notes on infer the extent of tephra fall deposits:
																### (1) Add one unit thickness to all the observations and log-transform them.
																### (2) Apply these transformed log-thickness to the method
																### (3) Specify the "levels <- c(1)" and apply it to "ContourLines2SLDF" function to find the unit thickness isopach
																### (4) The resultant contour is the inferred extent.
																### ***
																### The resultant contour is strongly influenced by the trend model,
																### as normally fewer observations are made on the distal and thin parts of a tephra fall deposit


###### End of demo ##########################################################################################################

























