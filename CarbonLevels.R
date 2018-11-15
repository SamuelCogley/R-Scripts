rm(list=ls())   # Clear workspace environment

#Rate constants
#Units: 1/year
#Given rate/Initial source carbon
f_terrPhoto = 0.1467
f_marinePhoto = 0.0533
f_terrResp = .0916
f_marineResp = .05
f_carbonDissolve = .1334
f_evap = .125
f_upwell = .00071
f_marineDeath = .005
f_plantDeath = .0916
f_plantDecay = .0367

#Unit: gigatons carbon/year
f_downwell = 23

#Integration constants 
dt = 0.01 # years; the time step
tf = 100 # year; the number of years the models simulates
numIter = tf/dt # Number of iterations

#Result vectors

#Flux vectors
terrPhoto = vector()
marinePhoto = vector()
terrResp = vector()
marineResp = vector()
carbonDissolve = vector()
evap = vector()
upwell = vector()
marineDeath = vector()
plantDeath = vector()
plantDecay = vector()
downwell = vector()

t = vector() #time

#source vectors
atmosphere = vector()
terrBio = vector()
oceanSurface = vector()
deepOcean = vector()
soil = vector()

#Initial carbon levels in reservoirs
#Gigatons of carbon
atmosphere[1] = 750 
terrBio[1] = 600
oceanSurface[1] = 800
deepOcean[1] = 38000
soil[1] = 1500

#Initial Conditions: source * rate

terrPhoto[1] = atmosphere * f_terrPhoto
marinePhoto[1] = atmosphere * f_marinePhoto
terrResp[1] = terrBio * f_terrResp
marineResp[1] = oceanSurface * f_marineResp
carbonDissolve[1] = atmosphere * f_carbonDissolve
evap[1] = oceanSurface * f_evap
upwell[1] = deepOcean * f_upwell
marineDeath[1] = oceanSurface * f_marineDeath
plantDeath[1] = terrBio * f_plantDeath
plantDecay[1] = soil * f_plantDecay
downwell[1] = f_downwell


t[1] = 0



for(i in 2:numIter ){

  t[i] = (i-1)*dt
  #flux equations
  
  terrPhoto[i] = terrPhoto[i-1] + dt * f_terrPhoto * atmosphere[i-1]
  marinePhoto[i] = marinePhoto[i-1] + dt * f_marinePhoto * atmosphere[i-1]
  terrResp[i] = terrResp[i-1] + dt * f_terrResp * terrBio[i-1]
  marineResp[i] = marineResp[i-1] + dt * f_marineResp * oceanSurface[i-1]
  carbonDissolve[i] = carbonDissolve[i-1] + dt *f_carbonDissolve * atmosphere[i-1]
  evap[i]= evap[i-1] + dt * f_evap * oceanSurface[i-1]
  upwell[i] = upwell[i-1] + dt * f_upwell * deepOcean[i-1]
  marineDeath[i] = marineDeath[i-1] + dt * f_marineDeath * oceanSurface[i-1]
  plantDecay[i] = plantDecay[i-1] + dt * f_plantDecay * soil[i-1]
  plantDeath[i] = plantDeath[i-1] + dt * f_plantDeath * terrBio[i-1]
  downwell[i] = downwell[i-1] + dt * f_downwell
  
  #source equations
  
  atmosphere[i] = (atmosphere[i-1] + (terrResp[i] +
                                       plantDecay[i] +
                                       evap[i]+
                                       marineResp[i]) -
                             (carbonDissolve[i] +
                              marinePhoto[i] +
                              terrPhoto[i]))
  
  terrBio[i] = (terrBio[i-1] + (terrPhoto[i])
                            - (terrResp[i] + plantDeath[i]))
  oceanSurface[i] = (oceanSurface[i-1] + (marinePhoto[i] +
                                           carbonDissolve[i] +
                                           upwell[i])
                                  - (marineResp[i] +
                                       evap[i] +
                                       downwell[i] +
                                       marineDeath[i]))
  deepOcean[i] = (deepOcean[i-1] + (downwell[i] +
                                     marineDeath[i])
                                - (upwell[i]))
  soil[i] = soil[i-1] + plantDeath[i] - plantDecay[i]

}

#Plot graphs of all reservoirs
plot( 
  x = t, 
  y = atmosphere, 
  type = "l", 
  xlab = "years", 
  ylab = "carbon (gigatons of Carbon)",  
  main = "Carbon in Atmosphere Over Time")

plot( 
  x = t, 
  y = terrBio, 
  type = "l", 
  xlab = "years", 
  ylab = "carbon (gigatons of Carbon)",  
  main = "Carbon in Terrestrial Biosphere Over Time")

plot( 
  x = t, 
  y = oceanSurface, 
  type = "l", 
  xlab = "years", 
  ylab = "carbon (gigatons of Carbon)",  
  main = "Carbon in Ocean Surface Over Time")

plot( 
  x = t, 
  y = deepOcean, 
  type = "l", 
  xlab = "years", 
  ylab = "carbon (gigatons of Carbon)",  
  main = "Carbon in Deep Ocean Over Time")

plot( 
  x = t, 
  y = soil, 
  type = "l", 
  xlab = "years", 
  ylab = "carbon (gigatons of Carbon)",  
  main = "Carbon in Soil Over Time")