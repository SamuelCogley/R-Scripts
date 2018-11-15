# Jordan Brown and Rahul Isaacs
# 12/5/14
# Pit Vipers

#source("RatDiffusion.R",print.eval=TRUE)

rm(list=ls())

# initialize rat
par(mar = rep(0,4))
m = 5;
n = 5;
diffusionRate=0.1;

x = floor(runif(1,min=1,max=n+1))
y = floor(runif(1,min=1,max=m+1))

rat = c(x,y)

hotSites = matrix(rat,ncol=2)
t=10;

# Diffusion Simulation
diffusionSim = function(m,n,diffusionRate,hostSites,t) {
  # Declare global variables hot and ambient
  # Initialize grid  
  bar = initBar(m,n,hotSites)
  bar  
  # Perform simulation  
  grids<-array(0, dim=c(m,n,t+1))
  grids[,,1]<-bar  
  for(i in 2:(t+1)) {    
    #     Extend matrix
    barExtended = reflectingLat(bar)
    #     Apply spread of heat function to each grid point
    bar = applyDiffusionExtended(diffusionRate,barExtended)
    #     reapply hot spots
    bar = applyHotCold(bar,hotSites)
    #    save new Matrix
    grids[,,i]<-bar    
  }
  return(grids)
}


# initialize the grid
initBar <-function(m,n,hotSites){
  utils::globalVariables(c("AMBIENT"))
  bar <-AMBIENT*(matrix(c(rep(1,m*n)), nrow = m))
  bar = applyHotCold(bar, hotSites)
  return(bar)}


# APPLYHOTCOLD return bar with hot and cold sites
applyHotCold = function(bar, hotSites, coldSites) {
  utils::globalVariables(c("HOT"))
  newbar = bar
  for(k in 1:length(hotSites[,1])){

    newbar[hotSites[k,1],hotSites[k,2]]<-HOT
  }
    return(newbar)
}


# REFLECTINGLAT returns extended lattice with reflecting boundary
# conditions
reflectingLat = function(lat){
  latNS = rbind(lat[1,],lat,lat[nrow(lat),])
  extLat = cbind(latNS[,1], latNS, latNS[,ncol(latNS)])
  return (extLat)}


# APPLYEXTENDED - Function to apply 
# diffusion(diffusionRate, site, N, NE, E, SE, S, SW, W, NW)
# site of matrix latExt and to return the resulting matrix
applyDiffusionExtended <- function(diffusionRate,latExt){
  m = length(latExt[,1]) - 2
  n = length(latExt[2,]) - 2
  newLat = matrix(rep(0,m*n), nrow=m)
  # calculate one column at a time because R is column-oriented
  for(j in 2:(n+1)){
    for(i in 2:(m+1)){
      site = latExt[i, j]
      N = latExt[i - 1, j]
      NE = latExt[i -1, j + 1]
      E = latExt[i, j + 1]
      SE = latExt[i + 1, j + 1]
      S  = latExt[i + 1, j]
      SW = latExt[i + 1, j - 1]
      W = latExt[i, j - 1]
      NW = latExt[i - 1, j - 1]      
      newLat[i - 1, j - 1] = diffusion(diffusionRate, site,N, NE, E, SE, S, SW, W, NW)
    }
  }
  return(newLat)
  }


# DIFFUSION new value at cell due to diffusion
diffusion <- function(diffusionRate, site, N, NE, E, SE, S, SW, W, NW){
  return ((1 - 8*diffusionRate) * site + diffusionRate*(N+NE+E+SE+S+SW+W+NW))
  }

# animate grids in gray with HOT being black, COLD more white
animDiffusionGray<- function(graphList){  
  utils::globalVariables(c("HOT"))
  lengthGraphList = dim(graphList)[3]	
  # set up grayscale map
  map = gray(0:HOT / HOT)	
  m = dim(graphList)[1]
  n = dim(graphList)[2]
  # determine window size
  # 1.6 is approximately the golden ratio; used so cell pictured as square
  fraction = n/m;
  dev.new(width = 2 * fraction, height = 2 * 1.6)
  for( k in 1:lengthGraphList) {
    dev.hold()
    g = t(graphList[,,k])  # transpose because of following comment
    #"image interprets the z matrix as a table of f(x[i], y[j]) values, so that the x axis corresponds to row number and the y axis to column number, with column 1 at the bottom, i.e. a 90 degree counter-clockwise rotation of the conventional printed layout of a matrix."
    # first row's image is on bottom, unlike figures in text
    image(HOT - g + 1, col = map, axes = FALSE)
    box()
    Sys.sleep(0.1)
    dev.flush()
  }	
}
# move rat
moveRat <- function(rat){
  x = rat[1,1]
  y = rat[1,2]
  xtst = runif(1)
  if(x==0){
    if(xtst<0.5){
      x=x+1
    }
  }
  else if(x==n){
    if(xtst<0.5){
      x=x-1
    }
  }
  else{
    if(xtst<1/3){
      x=x+1
    }
    else
    {
      if(xtst>2/3){
        x=x-1
      }
    }
  }
  ytst = runif(1)
  if(y==0){
    if(ytst<0.5){
      y=y+1
    }
  }
  else if(y==m){
    if(ytst<0.5){
      y=y-1
    }
  }
  else{
    if(ytst<1/3){
      y=y+1
    }
    else
    {
      if(ytst>2/3){
        y=y-1}
    }
  }
  newRat = c(x,y)    
  return newRat
}
# hotSites = matrix(rat,ncol=2)

# animate grids in color with hotter being more red, colder more blue
animDiffusionColor  <- function(graphList){  
  utils::globalVariables(c("HOT"))
  lengthGraphList = dim(graphList)[3]
  # set up color map
  map <- 1:(HOT + 1)
  for (i in 0:HOT) {
    amt <- i/HOT
    map[i + 1] <- rgb(1-amt, 0, amt)  # red-green-blue amounts
  }	 
  m = dim(graphList)[1]
  n = dim(graphList)[2]	
  # determine window size
  # 1.6 is approximately the golden ratio; used so cell pictured as square
  fraction = n/m
  dev.new(width = 2 * fraction, height = 2 * 1.6)	
  for( k in 1:lengthGraphList) {
    dev.hold()
    g = t(graphList[,,k])  # transpose because of following comment
    #"image interprets the z matrix as a table of f(x[i], y[j]) values, so that the x axis corresponds to row number and the y axis to column number, with column 1 at the bottom, i.e. a 90 degree counter-clockwise rotation of the conventional printed layout of a matrix."
    # first row's image is on bottom, unlike figures in text
    image(HOT - g + 1, col = map, axes = FALSE)
    box()
    Sys.sleep(0.1)
    dev.flush()
    cat("t = ", k,"\n")
  }	
}

#Declare global variables hot and ambient.
utils::globalVariables(c("HOT","AMBIENT","m","n"))
AMBIENT = 20.0 # degrees Celsius;
HOT = 37.0 # degrees Celsius;

########################
# Run Rat

grids = diffusionSim(m,n,diffusionRate, hotSites,t);
animDiffusionGray(grids)
#animDiffusionColor(grids)
