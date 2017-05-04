##############################################################################
################                                               ###############
################             GAM Randomise Wrapper             ###############
################           Angel Garcia de la Garza            ###############
################              angelgar@upenn.edu               ###############
################                 05/02/2017                    ###############
##############################################################################

suppressMessages(require(optparse))

##############################################################################
################                 Option List                   ###############
##############################################################################


option_list = list(
  make_option(c("-c", "--covariates"), action="store", default=NA, type='character',
              help="Full path to RDS covariate file.  
              Please include bblid and scanid in this file as a column each"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="Full path to output directory"), 
  make_option(c("-p", "--imagepaths"), action="store", default=NA, type='character',
              help="Name of the variable in the covariate file that contains the path to the images to be analyzed"), 
  make_option(c("-m", "--mask"), action="store", default=NA, type='character',
              help="Full path to mask"), 
  make_option(c("-i", "--inclusion"), action="store", default=NA, type='character',
              help="Name of inclusion variable on dataset. By default 1 means include. This will subset your rds file"),
  make_option(c("-u", "--subjId"), action="store", default=NA, type='character',
              help="subjID name on the covariates dataset"), 
  make_option(c("-f", "--ffull"), action="store", default=NA, type='character',
              help="Formula for full model, should only include the right hand side of the formula.
              Example: ~ stai_stai_tr+sex+s(age)+s(age,by=sex)"),
  make_option(c("-r", "--freduced"), action="store", default=NA, type='character',
              help="Formula for covariates of the reduced model, should only include the right hand side of the formula.
              Example: ~ stai_stai_tr+sex+s(age)"),
  make_option(c("-k", "--splits"), action="store", default=10, type='numeric',
              help="number of splits to divide the data in, default is 10. To minimize data usage"),
  make_option(c("-d", "--skipfourD"), action="store", default=FALSE, type='logical',
              help="Option to skip creation of fourdD image and look for it in the Analysis Directory.
              4D image must be labeled as 'fourd.nii.gz'. Will also skip smoothing step.
              Default (FALSE) means to not skip"),
  make_option(c("-t", "--thresh"), action="store", default=0.01, type='numeric',
              help="P-value threshold for cluster correction"),
  make_option(c("-n", "--nsim"), action="store", default=500, type='numeric',
              help="Number of simulations, default is 500"),
  make_option(c("-e", "--execute"), action="store", default=FALSE, type='logical',
              help="Whether to run the command default is to only print out command ")
  )

opt = parse_args(OptionParser(option_list=option_list))

for (i in 1:length(opt)){
  if (is.na(opt)[i] == T) {
    cat('User did not specify all arguments.\n')
    cat('Use gam_voxelwise.R -h for an expanded usage menu.\n')
    quit()
  }
}


print("##############################################################################")
print("################  Generalized Additive Model Randomise Script  ###############")
print("################            Angel Garcia de la Garza           ###############")
print("################              angelgar@upenn.edu               ###############")
print("################                 Version 3.0.1                 ###############")
print("##############################################################################")

##############################################################################
################                  Load Libraries               ###############
##############################################################################

print("Loading Libraries")

suppressMessages(require(ggplot2))
suppressMessages(require(base))
suppressMessages(require(reshape2))
suppressMessages(require(nlme))
suppressMessages(require(lme4))
suppressMessages(require(gamm4))
suppressMessages(require(stats))
suppressMessages(require(knitr))
suppressMessages(require(mgcv))
suppressMessages(require(plyr))
suppressMessages(require(oro.nifti))
suppressMessages(require(parallel))
suppressMessages(require(optparse))
suppressMessages(require(fslr))
suppressMessages(require(voxel))



##############################################################################
################              Declare Variables               ###############
##############################################################################

print("Reading Arguments")

subjDataName <- opt$covariates
OutDirRoot <- opt$output
namePaths <- opt$imagepaths
maskName  <- opt$mask
inclusionName <- opt$inclusion
subjID <- opt$subjId
fullFormula <- opt$ffull
redFormula <- opt$freduced
splits <- opt$splits
skipFourD <- opt$skipfourD
thresh <- opt$thresh
nsim <- opt$nsim
runCommand <- opt$run


##############################################################################
################            Load subject data                  ###############
##############################################################################

print("Loading covariates file")
subjData<-readRDS(subjDataName) ##Read Data
subset <- which(subjData[inclusionName] == 1) ##Find subset for analysis
subjData <- subjData[subset, ] #subset data


##############################################################################
################    Create Analysis Directory                  ###############
##############################################################################


print("Creating Analysis Directory")
OutDir <- paste0(OutDirRoot, "/n",dim(subjData)[1],"_",namePaths,"_",inclusionName,"_smooth",as.character(smooth))
dir.create(OutDir)
setwd(OutDir)

##############################################################################
################     Create and output fourd image             ###############
##############################################################################


if (!skipFourD) {
  print("Merging images and saving out a fourd image")
  subjList <- as.character(subjData[,grep(namePaths, names(subjData))])
  length.subj <- length(subjList)
  k <- splits
  break.subj <- ceiling(length.subj / k)
  
  subMergeNames <- "foo"
  for (i in 1:k) {
    if (i == k) {
      out <- paste0("fourd_",i,".nii.gz")
      fslmerge(subjList[(1 + (i-1)*break.subj):length.subj], direction="t", outfile=out, drop_dim=F)
      subMergeNames <- c(subMergeNames, out)
    } else {
      out <- paste0("fourd_",i,".nii.gz")
      fslmerge(subjList[(1 + (i-1)*break.subj):((i)*break.subj)], direction="t", outfile=out, drop_dim=F)
      subMergeNames <- c(subMergeNames, out)
    }
  }
  
  subMergeNames <- subMergeNames[-1]
  fslmerge(subMergeNames, direction="t", outfile="fourd.nii.gz")
  
  
  system('rm -f fourd_*.nii.gz')
  
  
} else {
  print("Skipping fourd image creation; Script will looking for file names fourd.nii.gz under first level directory")
}


system(paste0("scp ", maskName," ",OutDir, "/mask.nii.gz"), wait=T)
print("mask succesfully copied")


##############################################################################
################        Make Output Directory                  ###############
##############################################################################


print("Creating output directory")
outName <- gsub("~", "", fullFormula)
outName <- gsub(" ", "", outName)
outName <- gsub("\\+","_",outName)
outName <- gsub("\\(","",outName)
outName <- gsub("\\)","",outName)
outName <- gsub(",","",outName)
outName <- gsub("\\.","",outName)
outName <- gsub("=","",outName)
outName <- gsub("\\*","and",outName)
outName <- gsub(":","and",outName)

print("Creating output directory")
outNameRed <- gsub("~", "", redFormula)
outNameRed <- gsub(" ", "", outNameRed)
outNameRed <- gsub("\\+","_",outNameRed)
outNameRed <- gsub("\\(","",outNameRed)
outNameRed <- gsub("\\)","",outNameRed)
outNameRed <- gsub(",","",outNameRed)
outNameRed <- gsub("\\.","",outNameRed)
outNameRed <- gsub("=","",outNameRed)
outNameRed <- gsub("\\*","and",outNameRed)
outNameRed <- gsub(":","and",outNameRed)

outsubDir <- paste0("n",dim(subjData)[1],"gam_randomise_full_",outName,"_reduced_",outNameRed)

outsubDir<-paste(OutDir,outsubDir,sep="/")

logDir<-paste(OutDir,"logs",sep="/")

#Will return a warning if logDir and outsubDir already exist
dir.create(logDir)
dir.create(outsubDir)

system(paste('rm -f', file.path(outsubDir, '*')))

##############################################################################
################            Output summary files               ###############
##############################################################################

setwd(outsubDir)

write.table(subjData[, namePaths], paste0(namePaths,".csv"), row.names = F, col.names=FALSE)
write.table(subjData[, subjID], paste0(subjID,".csv"), row.names = F, col.names=FALSE)

print("Succesfully wrote paths and id files")


##############################################################################
################              Echo Arguments                   ###############
##############################################################################


system( paste0("echo Arguments are: >> ", outsubDir,"/logs.txt"))
system( paste0("echo Covariates file is: ", subjDataName,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Output directory is: ", OutDir,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Path name variable in covarites file is: ", namePaths,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Mask path is: ", maskName,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Inclusion variable name is: ", inclusionName,">> ", outsubDir,"/logs.txt"))
system( paste0("echo ID variable name is: ", subjID,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Formula for full model is: ", fullFormula,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Formula for reduced model is: ", redFormula,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Number of splits is: ", splits,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Skip 4D image creation is: ", skipFourdD,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Threshold for correction is: ", thresh,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Number of simulations is: ", nsim,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Run command is: ", outName,">> ", runCommand,"/logs.txt"))

###cleanup logdir
system(paste('rm -f', file.path(logDir, '*')))



##############################################################################
################              Create Design Matrix             ###############
##############################################################################

subjData$dummy <- rnorm(dim(subjData)[1])
# model matrices
X = model.matrix(gam(update.formula(fullFormula, "dummy ~ .") , data=subjData))
Xred = model.matrix(gam(update.formula(redFormula, "dummy ~ .") , data=subjData))


##############################################################################
################          Output Design and Contrasts          ###############
##############################################################################


## DESIGN AND CONTRASTS ##
# design file
n = nrow(X)
p = ncol(X)
p2 = p - ncol(Xred)
matfile = file.path(outsubDir, 'design.mat')
cat('/NumWaves\t', ncol(X), '\n/NumPoints\t', nrow(X), '\n/PPheights\t', paste(apply(X, 2, function(x) abs(diff(range(x))) ), collapse='\t'), '\n\n/Matrix\n', sep='', file=matfile)
write.table(X, append=TRUE, file=matfile, row.names=FALSE, col.names=FALSE)
  
# contrast file
confile1 = file.path(outsubDir, 'design.con') # for f-test
cons = matrix(0, nrow=p2, ncol=ncol(X))
cons[ cbind(1:(p2), which(! colnames(X) %in% colnames(Xred) ) ) ] = 1
cat('/ContrastName1\t temp\n/ContrastName2\t\n/NumWaves\t', ncol(X), '\n/NumPoints\t', nrow(cons), '\n/PPheights\t', paste(rep(1,ncol(cons)), collapse='\t'), '\n/RequiredEffect\t1\t1\n\n/Matrix\n', sep='', file=confile1)
write.table(cons, append=TRUE, file=confile1, row.names=FALSE, col.names=FALSE)
  
# fts file
ftsfile = file.path(outsubDir, 'design.fts')
fts = matrix(1, nrow=1, ncol=nrow(cons)) # ftest of all contrasts
cat('/NumWaves\t', nrow(cons), '\n/NumContrasts\t', 1, '\n\n/Matrix\n', sep='', file=ftsfile)
write.table(fts, append=TRUE, file=ftsfile, row.names=FALSE, col.names=FALSE)
  

##############################################################################
################          Generate and Run command             ###############
##############################################################################

##change mergednifti
##change mask file 

# t distribution is two tailed, F is one tailed. -x outputs voxelwise statistics -N outputs null distribution text files
# F-test
fcmd = paste('randomise -i', mergednifti, '-m', maskfile, '-o', file.path(outsubDir, 'randomise'), '-d', matfile, '-t', confile1, '-f', ftsfile, '--fonly -F', qf( (1-thresh),df1=p2, df2=(n-p) ), '-x -N -n', nsim, '--uncorrp' )


##Change run
if(run){
  system(fcmd)
}

print(fcmd)

print("Write t-maps and p-maps is done")
