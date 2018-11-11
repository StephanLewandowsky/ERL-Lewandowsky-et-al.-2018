# check yrw and col vals
#HiatusHistoryModels.R: Code to do all model-based analyses for second ERL paper 2017 - 2018.
library(fArma)
library(timsac)
library(lattice)
library(stats)  #putting this ahead of dplyr so filter is masked by dplyr
library(car)
library(diagram)
library(latticeExtra) #for layer
library(sp) #for layer sp.points
library(doParallel)
library(foreach)
library(effects)
library(gplots)
library(rgl)
library(nlme)
library(grid)
library(tidyverse)

rm(list=ls())
graphics.off()
setwd("C:/Users/Lewan/Documents/Research Projects/Climate Change/_climate with James/HiatusHistory")
source("hiatushistoryfuncs.r") #define all required functions
source("hiatushistorymodelfuncs.r")


#######################################################################################
#function to read all available versions of a data set, 
#   no matter how many there are.
#   Alphabetical order of data files must correspond to their chronological order
#   Return all versions of data sets in a list.
readallversions <- function(datadir,ds) {
  wd <- paste(datadir,ds,sep="/")
  tbr<- paste (wd,list.files(wd,pattern="*.temp"),sep="/")
  r  <- vector("list",length(tbr))       #build a list of all versions of the data 
  for (i in c(1:length(tbr))) {          #now read each file and return result in a list
    temp <- read.table(tbr[i])
    r[[i]] <- data.frame(t=temp[,1],anom=temp[,2])
  }
  return(r)
}

#function to extract model ID so different models can be considered separately
readmodelIDs <- function (datadir,ds) {
  wd  <- paste(datadir,ds,sep="/")
  fns <- list.files(wd,pattern="*.temp")
  return( as.vector(sapply(fns,FUN=function(x) strsplit(x,"_")[[1]][2])) )
}

# function to ensure that year always runs to late December, even if provided as integer
mkyr <- function(y) {return(ifelse(abs(floor(y)-y) > 1e-8, y, y+.99))}
# function to verify that a dataset is annualized to guard against silly errors
verifyannual <- function(ds) {try(if(abs(floor(max(ds$t))-max(ds$t)) > 1e-8) stop("DS not annualized")) }

#function rebaselines anomalies to the desired baseline period 
#  and cleans the data to fall within the desired window of use
cleandata <- function(rawvariants) {
  cln <- rawvariants
  for (i in c(1:length(cln))) {
    cln[[i]] <- dplyr::filter(cln[[i]],t>=yrsused[1] & t<=mkyr(yrsused[2]))
    blmean <- mean(dplyr::filter(cln[[i]],t>=baselnpd[1] & t<=mkyr(baselnpd[2]))$anom)
    cln[[i]]$anom <- cln[[i]]$anom - blmean
  }
  return(cln)
}

#function annualizes data set for convenience
annualize <- function(ds) {
  ca <- aggregate(anom ~ floor(t), data=ds, FUN=mean)
  names(ca)[1] <- "t"
  return(ca)
}

#http://stackoverflow.com/questions/21011672/automatically-add-variable-names-to-elements-of-a-list
listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}

#######################################################################################
#major  driver variables for the entire analysis. Some are packaged into lists below.
datadir   <- "./HiatusHistory TemperatureDataV1.0"
modeldir  <- "./HiatusHistory AllModels"
outputdir <- "./HiatusHistory Model Output"

#now decide how to adjust forcings--and update file name accordingly
forceadjf <- "schmidtadj.dat"
adjname   <- "schmadj"

# forceadjf <- "huberadjback.dat" #<--not used for paper (13/1/18)
# adjname   <- "hubackadj"

#forceadjf <- "huberadjvolc.dat"
#adjname   <- "hubvolcadj"

yrsused  <- c(1880,2016.999)  #inclusive
yrslscap <- c(1965,2016.999)  #inclusive for landscaping
baselnpd <- c(1981,2010.999)  #inclusive: note this may be different from ERL I c(1981,2010.999).
boftrnd  <- c(1998,1970)      #[1]=beginning of trend for Cowtan retro graph,
                              #[2]=beginning of global warming trend for continuous slopes
mintrnd  <- 10                #minimum trend length for histocond figures
maxyrlscap <-24               #maximum years back in time for landscaping
nMC        <- 100            #number of Monte Carlo realizations

#obtain data sets by checking directories in the data directory
datasets <- list.dirs(datadir,full.names=FALSE,recursive=FALSE) 
names(datasets) <- datasets
datasets <- datasets[-which(datasets=="NOAA")]  #eliminate NOAA because of masking issue

#now obtain model variants by checking directories in the appropriate model directory
cmipversions <- list.dirs(modeldir,full.names=FALSE,recursive=FALSE) 
names(cmipversions) <- cmipversions

#manually add first-available-date
ds1strel <- rep(1879.999,length(datasets))   #default availability from early on
names(ds1strel) <- datasets
ds1strel["CW"]  <- 2013.80                   #November 2013 (for GT test)
ds1strel["BERKELEY"]  <- 2013.90             #December 2013 (for GT test) [Actually March 2014 but that wouldnt be meaningful]
#and had3 ends ...
had3endyr <- 2012
#now keep track of changepoints (as per Cahill et al.)
cps <- rep(boftrnd[2],length(datasets))      #just in case some explicitly missing (e.g. BERKELEY)
names(cps) <- datasets
cps["GISTEMP"] <- 1970 #SD=3
cps["CW"]      <- 1974 #SD=2
cps["HADCRUT"] <- 1974 #SD=2


#model versions: online publication dates of the relevant papers
blenddate <- 2015.625  #Blending TAS/TOS became available
fadjdate  <- 2014.125  #Schmidt et al forcing adjustments (use same for Huber for now)

#plotting variables such as color and so on. Some packaged into list below.
plcol     <- c("red","purple4","darkblue","dark green")
plthincol <-c("mistyrose","plum","lightblue1","palegreen3")
plchar    <- c(21,22,23,24)
names(plcol) <- names(plthincol) <- names(plchar) <- datasets


#######################################################################################
#now read, clean, and rebaseline all available variants of each data set
mcresults <- alldata  <- vector("list",length(datasets))
names(mcresults) <- names(alldata) <- datasets      #permit convenient indexing by name
for (i in c(1:length(datasets))) {
  rawvariants   <- readallversions(datadir,datasets[i])
  alldata[[i]]  <- cleandata(rawvariants)    
}

#read and duplicate after adjustment versions
allmodels  <- vector("list",length(cmipversions))
for (i in c(1:length(cmipversions))) {
  cmipraw <- readallversions(modeldir,cmipversions[i])
  allmodels[[i]]  <- cleandata(cmipraw)
}
modelIDs <- readmodelIDs(modeldir,cmipversions[1])
#first create adjusted models
allmodels <- c(allmodels,makeadj(allmodels,forceadjf))
names(allmodels) <- c(cmipversions,paste(cmipversions,adjname,sep="_"))
#now create a concatenated mask, unadjusted and not blended, for had3 and had4
allmodels <- c(allmodels,list(concatmasks(allmodels,had3endyr)))                      
names(allmodels)[length(allmodels)] <- "concatmasked_had34_tas_hemi"
cmipversions <- names(allmodels)
names(cmipversions) <- cmipversions
################################################################################################################################################################
############ anything above here is preface for general setup. Segments below can be executed directly after the preface has been run ##########################
################################################################################################################################################################


#######################################################################################
#get the literature corpus and analyze the subset of papers that deal with model-obs
#discrepancy
flit <- getlit() 
lit <- subset(flit,M.O==2)
# x11()
# hist(lit$Start)
# x11()
# hist(lit$duration)
summary(lit$Start)
summary(lit$duration)

#use the unblended and unmodified original CMIP5 for all
usemmm4lit <- 1  #if this is 1, then use multi model mean, otherwise the ensemble members
if (usemmm4lit) {
  cmip4lit <- list(data.frame(t=seq(yrsused[1],floor(yrsused[2])),anom=cmipenvelope(allmodels[["global_tas"]])$ensmean))
  fnmarker <- "_mmm"
} else {
  fnmarker <- "_ensmems"
  cmip4lit <- allmodels[["global_tas"]]
}

#now process each paper
pslopes <- prefslopes <- qofslope <- NULL
for (i in 1:dim(lit)[1]) {
  print(c(i,lit$ds[i,]))
  nds <- unique(lit$ds[i,!is.na(lit$ds[i,])])
  for (ds4histo in nds[nds != "NOAA"]) {
    apaper     <- getmodstat4paper(cmip4lit,ds4histo,lit$Start[i],lit$End[i])
    pslopes    <- c(pslopes,apaper$pslope)
    prefslopes <- c(prefslopes,apaper$refslopes)
    qofslope   <- c(qofslope,apaper$q)
  }
}
#plot the results of literature analysis
y1 <- hist(pslopes, plot=FALSE,breaks="FD")
y1$density <- y1$density/sum(y1$density)
x11()
plot( y1, col=rgb(0,0,1,1/4), freq=FALSE,  xlab="K/decade",ylab="Density",main=NULL,
      ylim=c(0,.45),xlim=c(-.2,.7))
os <- hist(prefslopes, plot=FALSE,breaks="FD")
os$density <- os$density/sum(os$density)
plot( os, col=rgb(1,0,1,1/4), freq=FALSE,  add=TRUE)
abline(v=0,lty="dashed")

savePlot(filename = paste(outputdir,"/modlitHistoSlopes",fnmarker,".pdf",sep=""), type = "pdf", device = dev.cur(), restoreConsole = TRUE)


#######################################################################################
#first plot the latest annualized means with model projections (no historical conditioning)
p2f <- 0
for (whichversion in cmipversions[-grep("masked_had3",cmipversions)]) {
    if (length(grep("had",whichversion))) {
      plotcmip(whichversion,"HADCRUT")
    } else {
      plotcmip(whichversion,datasets[-grep("HADCRUT",datasets)])
    }
}



#######################################################################################
#now plot historically conditioned temperatures for various model variants
#   conditioning here is for both data (always, though broken lines provide 
#   counterfactual) and models (only if requested by name of version being passed)
p2f         <- 1 #1 = to file, 0 = screen only
tofr        <- 1 #broken trend = 1, continuous trend = 2
trendlength <- 15
startyr     <- 1970  #use this when multiple data sets other than HADCRUT are plotted
endyr       <- 2016

#unmasked results (so no hadCRUT obs): historical 
plmodsanter(alldata,c("GISTEMP"),allmodels,"historical",adjname,
            tofr,plcol,trendlength,startyr,endyr,p2f)
#unmasked results: retrospective from latest model version (blended and adjusted)
plmodsanter(alldata,c("BERKELEY","CW","GISTEMP"),allmodels,paste("global_blend_",adjname,sep=""),adjname,
            tofr,plcol,trendlength,startyr,endyr,p2f)

#masked results (hadCRUT only): historical, using concatenated masks for model
plmodsanter(alldata,c("HADCRUT"),allmodels,"histo_concat",adjname,
            tofr,plcol,trendlength,cps["HADCRUT"],endyr,p2f)  #note different AGW onset

#masked results (hadCRUT only): retrospective from latest masked and updated model
plmodsanter(alldata,"HADCRUT",allmodels,paste("masked_had4_blend_hemi_",adjname,sep=""),adjname,
            tofr,plcol,trendlength,cps["HADCRUT"],endyr,p2f)  #note different AGW onset

#ignoring adjustments for now
# plmodsanter(alldata,"HADCRUT",allmodels,paste("concatmasked_had34_tas_hemi","",sep=""),adjname,
#             tofr,plcol,trendlength,cps["HADCRUT"],endyr,p2f)  #note different AGW onset
# plmodsanter(alldata,"HADCRUT",allmodels,paste("masked_had4_blend_hemi","",sep=""),adjname,
#             tofr,plcol,trendlength,cps["HADCRUT"],endyr,p2f)  #note different AGW onset



#######################################################################################
#Monte Carlo analysis: compute residuals model vs. observations, then simulate multiple
#  realizations of the same structure.
#  Because continuous time series are modeled, historical conditioning can only be up to 
#  'end year' being considered
p2f         <- 1
tofr        <- 1  #broken trend = 1, continuous trend = 2
modwt       <- 1 #if 1 then average computed within each model first, if 0 then average computed across all runs
trendlength <- 15

#There are 3 endyears that are of interest: 2012 (end of Had3), 2013 (old forcings), and 2016 (latest)
endyrs4MC <- c(2012,2012,2012, 2013,2013,2013, 2016,2016,2016)
dataptrs <- c("HADCRUT","GISTEMP","CW","HADCRUT","GISTEMP","CW","HADCRUT","GISTEMP","CW")
for (ptr in 1:length(endyrs4MC)) {  
  thisds <- dataptrs[ptr]
  endyr  <- endyrs4MC[ptr]
  tofm   <- gettofm(thisds,endyr,adjname)
  #call Monte Carlo to generate realizations with the correct noise structure (but no mmm added in)
  mcresults[thisds] <- mobsmontecarlowrapper(allmodels,modelIDs,alldata[[thisds]],endyr,nMC,tofm,modwt) 
  #mc is done on residuals alone, with no multi-model mean. Need to add back in mmm--global variable is set in function
  mcresults[thisds]  <- list(lapply(mcresults[[thisds]], FUN=function(x) data.frame(t=x$t,anom=x$anom + globalmmm$anom)))
  #plot realization envelope, mmm, and data (data historically conditioned to end year)
  plotmobsmc(alldata,
             annualize(getmmm (allmodels,tofm,endyr,modelIDs,modwt))$anom,  #compute MMM on the fly
             mcresults[[thisds]],thisds,endyr,tofm)
  #plot trends... data are historically conditioned (but realizations are not)
  plmodsanter(alldata,thisds,mcresults,thisds,adjname,
              tofr,plcol,trendlength,cps[thisds],endyr,p2f,1)  #last arg=1 flags Monte Carlo, so no MMM but spaghetti sample
} #end of condition 




#######################################################################################
#Model Monte Carlo: now do landscaping analysis
tofr            <- 1   #1==broken, 2==continuous trends
modwt           <- 1 #if 1 then average computed within each model first, if 0 then average computed across all runs
trendlength     <- 15
nMC4ls          <- 1000
pausegrid       <- matrix(NA,110,3)  #manually computed max, value is way higher than needed
earliestvantage <- 2007
present         <- 2016
minyr           <- 10  #check if this should be 9 (or yrw+1 in landscape column values below)
yrw             <-c(minyr:19) 
k <- 0 
for (p2 in c(earliestvantage:present)) {
  for (yback in yrw) {
    if (p2-yback+1 >= 1998) {
      k <- k+1
      pausegrid[k,] <- c(p2-yback+1,mkyr(p2),yback)
    }
  }
}
pausegrid <- na.omit(pausegrid) 
ncells <-dim(pausegrid)[1]

gridds <- c("GISTEMP","HADCRUT") 
gridmcresults <- vector("list",length(gridds))
names(gridmcresults) <- gridds
for (i in gridds) {                     # gridds do all the historically-available data sets
  print(i)
  gridmcresults[[i]] <- vector("list",ncells)
  for (k in 1:ncells) {                 #run through grid
    print(k)                            #provide sign of life....
    pauseydec  <- pausegrid[k,1:2]
    tpres      <- floor(pauseydec[2])  #this ensures historical conditioning
    realizatns <- mobsmontecarlowrapper(allmodels,modelIDs,alldata[[i]],tpres,nMC4ls,
                                        gettofm(i,tpres,adjname), #obtain correct type of model on the fly
                                        modwt)
    gridmcresults[[i]][[k]] <- examineMC(alldata[[i]],cps[i],realizatns,pauseydec,tpres,nMC4ls,tofr) 
  }
}
save(gridmcresults,file=paste("modelgridmc_noMMM_",adjname,"_",tofr,".RData",sep=""))

#direct entry here possible without performing MC --- tofr defined above
#load(file=paste("modelgridmc_noMMM_",adjname,"_",tofr,".RData",sep=""))

#unpack results for landscaping plot
pfromallreg <- 3   #this determines which type of p values to use: 1=allreg,2=pauseonly,3=agwonly
for (i in gridds) {  #do all the historically-available data sets gridds
  pdev <- magdev <- matrix(NA,10,10)
  for (k in 1:ncells) { #run through grid and fill matrix in landscape-conducive order
    row    <-  (floor(pausegrid[k,2]) - earliestvantage + 1)
    column <-  (floor(pausegrid[k,2]) - pausegrid[k,1] - 8)
    
    pdev[row,column]   <- switch(pfromallreg,
                                 gridmcresults[[i]][[k]]$pcpause,
                                 gridmcresults[[i]][[k]]$pcpauseponly,
                                 gridmcresults[[i]][[k]]$pcpauseagwonly)
    magdev[row,column] <- pdev[row,column]   
    
    #print stuff for debugging
    pauseyears <- pausegrid[k,1:2]
    tpres      <- floor(pauseyears[2])  
     print(c(round(trunc(pauseyears),0),magdev[row,column],
             round(lm(anom~t,data=annualize(subset(alldata[[i]][[length(alldata[[i]])]],t>=pauseyears[1] & t<=tpres+.999)))$coefficients[2]*10,4)))
  }

  #now do landscaping plot
  minmaxscale <- c(0,100) #c(-.3,.3) #symmetrical to force white being at zero
  mags4p <- magdev  #create copy that can be altered for printing
  mags4p[mags4p < minmaxscale[1]] <- minmaxscale[1]
  mags4p[mags4p > minmaxscale[2]] <- minmaxscale[2]
  
  pvals4p <- pdev
  print(min(pvals4p,na.rm=TRUE))
  
  xysig <- which(pvals4p <= 5, arr.ind = TRUE)     #row & column of significance
  if (length(xysig) == 0) {xysig<-matrix(NA,1,2)}  #avoid failure in printing points
  xy4num <- which(!is.na(pvals4p), arr.ind = TRUE) #row & column of printing numbers

  x11()  #if embedded in function, levelplot needs a print call
  print(
    levelplot(
      main=i,
      mags4p,
      row.values = c(earliestvantage:present),
      column.values = yrw , #+ 1, #should this +1 be here?
      xlab = "Vantage year",
      ylab = "Years included",
      at = do.breaks(minmaxscale, 100),
      col.regions = colorRampPalette(c("red", "yellow", "yellow3", "white", "green", "green2", "green4"), 
                                     space = "rgb"),
        
      #add monte carlo numbers and 'significance' indicators
      panel = function(...) {
        panel.levelplot(...)
        panel.text(
          xy4num[, 1] + earliestvantage - 1,
          xy4num[, 2] + minyr - 1,
          round((pvals4p[xy4num]),0)
        )
        grid.points(
          xysig[, 1] + earliestvantage - 1,
          xysig[, 2] + minyr - 1,
          pch = 1,
          gp = gpar(col = "yellow",lwd=1.5),
          size = unit(1.9, "char"))
      }
    )
  )
     savePlot(filename = paste(outputdir,"/modelmclandsc_reg",tofr,"_data", i,"_typep",pfromallreg,"_",adjname,".pdf",sep=""),type = "pdf",
                 device = dev.cur(),restoreConsole = TRUE)
} #end of data sets


#debugging bits down here....
ds <- "CW"
for (k in 1:ncells) {                 #run through grid
  pauseyears <- pausegrid[k,1:2]
  tpres      <- floor(pauseyears[2])  
  tl <- floor(pauseyears[2]-pauseyears[1]+1)
  mymmm<-annualize(getmmm (allmodels,gettofm(ds,2016,adjname),tpres,modelIDs,0))
  print(c(trunc(pauseyears),tl,
          round(lm(anom~t,data=subset(mymmm,t>=pauseyears[1] & t<=tpres+.999))$coefficients[2]*10,4),
          gettrends(floor(pauseyears[2]),startyr,tl,getobs(alldata[[ds]],tpres,1)),
          round(lm(anom~t,data=subset(annualize(getobs(alldata[[ds]],tpres,1)),t>=pauseyears[1] & t<=tpres+.999))$coefficients[2]*10,4)))
}  

#######################################################################################
#demonstrate multiple-testing problem using synthetic data alone
# Note: thisds is the last one unless loop above is skipped
tofr <- 2    #1==broken, 2==continuous trends
compsim <- illustratemt(mcresults[[thisds]],trendlength,thisds,tofr)
save(compsim,mcresults,thisds,trendlength,tofm,file=paste("multipletestingsim_",as.character(tofr),".RData",sep="")) #or load this again to reconvene here

x11()
pts2plot <- c(seq(.05,1,.05))
hexa <- c(-0x0031L,-0x0032L,-0x0044L,-0x00BCL, -0x00BDL)
hexapex <- c(.75,.75,.75,1,1)

#set up random control
ltran4p <- sapply(pts2plot,FUN=function(x) sum(compsim$ltranslope<x)/length(compsim$ltloslope))
plot(pts2plot,ltran4p,las=1,xlim=c(0,1),ylim=c(0,1),xlab="Proportion expected",ylab="Proportion observed",type="l",lty="dashed")
points(pts2plot,ltran4p,pch=21,bg="red")

#plot the cherry-picked results
for (j in 1:3){ #skip to avoid clutter
  ltlowish4p <- sapply(pts2plot,FUN=function(x) sum(compsim$ltlowishslope[,j]<x)/length(compsim$ltloslope))
  #ltlowishcor4p <- sapply(pts2plot,FUN=function(x) sum(compsim$ltlowishslopecor[,j]<x)/length(compsim$ltloslope))
  lines(pts2plot,ltlowish4p,lty="dotdash")
  points(pts2plot,ltlowish4p,pch=21,bg="white",cex=2.9)
  points(pts2plot,ltlowish4p,pch=hexa[j],bg="white",cex=hexapex[j])
  
  # lines(pts2plot,ltlowishcor4p,col="green")
  # points(pts2plot,ltlowishcor4p,pch=21,bg="green",cex=2.9)
  # points(pts2plot,ltlowishcor4p,pch=hexa[j],bg="green",cex=hexapex[j])
}
abline(c(0,1),lty="solid")

legend(.7,.3,c("Random","Lowest","Second lowest","First decile (Q10)"), 
       col=c("black","black","black","black"),
       #lty=c("dashed","dotdash","dotdash","dotdash"),
       pch=c(21,-0x0031L,-0x0032L,-0x0044L),
       pt.bg=c("red","black","black","black","black"))
savePlot(filename = paste(outputdir,"/typeIinflate_",as.character(tofr),".pdf",sep=""),type = "pdf",device = dev.cur(),restoreConsole = TRUE)

# remaining plot options
# ltlo4p <- sapply(pts2plot,FUN=function(x) sum(compsim$ltloslope<x)/length(compsim$ltloslope))
# 
# ltlocor4p <- sapply(pts2plot,FUN=function(x) sum(compsim$ltloslopecor<x)/length(compsim$ltloslope))
# lines(pts2plot,ltlocor4p,col="red")
# points(pts2plot,ltlocor4p,pch=22,bg="red",cex=2)