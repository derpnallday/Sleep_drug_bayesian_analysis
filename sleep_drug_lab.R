#below is code from the yMet-Xnom1grp-Mrobust-Example using sleep_drug.csv


#===============================================================================
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Load The data file 
myDataFrame = read.csv( file="sleep_drug.csv" )

# We are loading the data for each drug group

Drug1 = myDataFrame$L.Hyoscyamine
Drug2 = myDataFrame$L.Hyoscine
Drug3 = myDataFrame$R.Hyoscine



#-------------------------------------------------------------------------------
#explore the distribution of the data
par(mfrow=c(2,2))

hist(Drug1, 
     breaks=16, 
     main="Plot of L Hysocymine",
     xlab="extra hours of sleep",
     xlim=c(-2,6))

hist(Drug2, 
     breaks=16, 
     main="Plot of L Hyoscine",
     xlab="extra hours of sleep",
     xlim=c(-2,6))

hist(Drug3, 
     breaks=16, 
     main="Plot of Hyoscyamine R",
     xlab="extra hours of sleep",
     xlim=c(-2,6))

#------------------------------------------------------------------------------- 
# Optional: Specify filename root and graphical format for saving output.
# Otherwise specify as NULL or leave saveName and saveType arguments 
# out of function calls.
fileNameRoot1 = "Drug1-Jags-" 
fileNameRoot2 = "Drug2-Jags-" 
fileNameRoot3 = "Drug3-Jags-"

graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-Xnom1grp-Mrobust.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chains for 3 groups:
mcmcCoda1 = genMCMC( data=Drug1 , numSavedSteps=20000 , saveName=fileNameRoot1 )
mcmcCoda2 = genMCMC( data=Drug2 , numSavedSteps=20000 , saveName=fileNameRoot2 )
mcmcCoda3 = genMCMC( data=Drug3 , numSavedSteps=20000 , saveName=fileNameRoot3 )


#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda1) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda1 , parName=parName , 
            saveName=fileNameRoot1 , saveType=graphFileType )
}

parameterNames = varnames(mcmcCoda2) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda2 , parName=parName , 
            saveName=fileNameRoot2 , saveType=graphFileType )
}

parameterNames = varnames(mcmcCoda3) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda3 , parName=parName , 
            saveName=fileNameRoot3 , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of 1:
summaryInfo = smryMCMC( mcmcCoda1 , 
                        compValMu=3 , ropeMu=c(2.5,3.5) ,
                        compValSigma=1 , ropeSigma=c(0.5,1.5) ,
                        compValEff=0.0 , ropeEff=c(-0.1,0.1) ,
                        saveName=fileNameRoot1 )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda1 , data=Drug1 , 
          compValMu=2 , ropeMu=c(1.5,2.5) ,
          compValSigma=1 , ropeSigma=c(0.5,1.5) ,
          compValEff=0.2 , ropeEff=c(0.1,0.3) ,
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot1 , saveType=graphFileType )



# Get summary statistics of 2:
summaryInfo = smryMCMC( mcmcCoda2 , 
                        compValMu=0 , ropeMu=c(1,2) ,
                        compValSigma=2 , ropeSigma=c(5,6) ,
                        compValEff=0.0 , ropeEff=c(-0.1,0.1) ,
                        saveName=fileNameRoot2 )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda2 , data=Drug2 , 
          compValMu=2 , ropeMu=c(1.5,2.5) ,
          compValSigma=1 , ropeSigma=c(0.5,1.5) ,
          compValEff=0.2 , ropeEff=c(0.1,0.3) ,
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot2 , saveType=graphFileType )


# Get summary statistics of 3:
summaryInfo = smryMCMC( mcmcCoda3 , 
                        compValMu=0 , ropeMu=c(1,2) ,
                        compValSigma=2 , ropeSigma=c(5,6) ,
                        compValEff=0.0 , ropeEff=c(-0.1,0.1) ,
                        saveName=fileNameRoot3 )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda3 , data=Drug3 , 
          compValMu=2 , ropeMu=c(1.5,2.5) ,
          compValSigma=1 , ropeSigma=c(0.5,1.5) ,
          compValEff=0.2 , ropeEff=c(0.1,0.3) ,
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot3 , saveType=graphFileType )
#===============================================================================



