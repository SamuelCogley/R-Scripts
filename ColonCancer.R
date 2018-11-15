rm(list=ls())   # Clear workspace environment

#RATE CONSTANTS: unknown units
cat_activator = 0.01 #Interaction between B_cat  & ActiveReceptor
deg_constant = 0.01 # Constant for degradation of BetaCatenin
TF = .00059 # Transcription Factor, interaction between B_cat_t & TCF
TF_Prom = .005 #Constant for Cell Division Promotion
TF_Inhib = .005 #Constant for Cell Division Inhibition
WNT_FZL = .00475 #Const for interaction between WNT and FZL


#INTEGRATION CONSTANTS: unknown units
dt = 0.01 #  the time step, set arbitrarily
finalTime =  75# Set arbitrarily
numIter = finalTime/dt # Number of iterations

#RESULT VECTORS, no units
activeReceptor = vector()
betaCat = vector()
betaCatT = vector()
betaCatTCF = vector()
cellDivProm = vector()
cellDivInhib = vector()
FZL = vector()
targetGeneExp = vector()
TCF = vector()
WNT = vector()
t = vector()

#INITIAL VECTOR VALUES, no units
activeReceptor[1] = 0
betaCat[1] = 20
betaCatT[1] = 0
betaCatTCF[1] = 6
cellDivProm[1] = 20
cellDivInhib[1] = 20
targetGeneExp[1] = 10
FZL[1] = 10
TCF[1] = 10
WNT[1] = 5.616*0
t[1] = 0

#MAIN LOOP
for(i in 2:numIter ){
  t[i] = (i-1)*dt
  
  WNT[i] = WNT[i-1] - ( WNT_FZL * WNT[i-1] * FZL[i-1]) * dt
  FZL[i] = FZL[i-1] - ( WNT_FZL * WNT[i-1] * FZL[i-1]) * dt
  
  activeReceptor[i] = activeReceptor[i-1] + 
    ( ( WNT_FZL * WNT[i-1] * FZL[i-1] ) -
    ( cat_activator * activeReceptor[i-1] * betaCat[i-1] ) ) * dt
  
  betaCat[i] = betaCat[i-1] -
    ( deg_constant * betaCat[i-1] + 
        cat_activator * activeReceptor[i-1] * betaCat[i-1]) * dt
  
  betaCatT[i] = betaCatT[i-1] + 
    ( (cat_activator * activeReceptor[i-1] * betaCat[i-1]) -
        (TF * betaCatT[i-1] * TCF[i-1]) ) * dt
  TCF[i] = TCF[i-1] - ( TF * betaCatT[i-1] * TCF[i-1]) * dt
  
  betaCatTCF[i] = betaCatTCF[i-1] +
    ( (TF * betaCatT[i-1] * TCF[i-1]) -
      ( ( TF_Prom * betaCatTCF[i-1] * targetGeneExp[i-1] ) +
          (TF_Inhib * targetGeneExp[i-1])/(betaCatTCF[i-1] + 0.01) )) * dt
  targetGeneExp[i] = targetGeneExp[i-1] - ( TF_Prom * betaCatTCF[i-1] * targetGeneExp[i-1] + 
                                              (TF_Inhib * targetGeneExp[i-1])/(betaCatTCF[i-1] + 0.01)) * dt
  
  
  
  cellDivProm[i] = cellDivProm[i-1] + (TF_Prom * betaCatTCF[i-1] * targetGeneExp[i-1]) * dt
  cellDivInhib[i] = cellDivInhib[i-1] + ((( TF_Inhib * targetGeneExp[i-1] )/(betaCatTCF[i-1] + 0.01)) )* dt



  
}

#PLOTS
plot(t,
     WNT,
     type="l",
     xlab="Time",
     ylab="Amount",
     lty=1,
     col="blue",
     #main = "Inhibitors/Promoters vs Time",
     ylim=c(19,26)
     ) 

lines(t,
      cellDivProm,
      lty=1,
      col="blue")
lines(t,
      cellDivInhib,
      lty=1,
      col="red")

legend(0, 26, 
       c("Promotion", "Inhibition"), 
       lty = c(1,1,1),
       col = c("blue", "red", "blue"),
       cex = .65) 

# 
# 
# plot(t,
#      cellDivProm,
#      type="l",
#      xlab="Time",
#      ylab="Cell Divisison Promoters",
#      lty="dashed",
#      col="blue",
#      main = "Cell Divisison Promoters vs Time")
# plot(t,
#      cellDivInhib,
#      type="l",
#      xlab="Time",
#      ylab="WNT amout",
#      lty="dashed",
#      col="blue",
#      main = "Cell Divisison Inhibitors vs Time")
# plot(t,
#      betaCat,
#      type="l",
#      xlab="Time",
#      ylab="Beta Cat amout",
#      lty="dashed",
#      col="blue",
#      main = "Beta Cat vs Time")
# plot(t,
#      betaCatT,
#      type="l",
#      xlab="Time",
#      ylab="betaCatT",
#      lty="dashed",
#      col="blue",
#      main = "betaCatT vs Time")
# plot(t,
#      betaCatTCF,
#      type="l",
#      xlab="Time",
#      ylab="betaCatTCF",
#      lty="dashed",
#      col="blue",
#      main = "betaCatTCF vs Time")