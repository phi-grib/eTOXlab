# load packages
suppressPackageStartupMessages(library(randomForest))

#set seed
set.seed(1992)
# read arguments from command line
arg <- commandArgs( trailingOnly = T )
vpath <- arg[1]

# define files
fRdata = paste(vpath,"/modelR.Rdata",sep="")
finfo  = paste( vpath,'/infoModelR.csv',sep="")

# clean results of previous calculations
file.remove(fRdata)
file.remove(finfo)

# load X and Y
Y <- as.factor(read.table(paste( vpath,'/tmpY.csv',sep=""))[,1])
X <- read.table(paste( vpath,'/tmpX.csv',sep="") ,sep=",",header=F)

# tune Random Forest  model and return best model
mdl <- tuneRF(X, Y ,  plot=F ,doBest=TRUE , 
            improve = 0.05 , 
            stepFactor=2,
            ntreeTry=50)
            
# alternatively build a single Random Forest model
# mdl <- randomForest(X,Y,ntree=500,...)

# save model
save(mdl, file= fRdata )

# save confusion matrix
df <- data.frame( TP=mdl$confusion[2,2],
                  TN=mdl$confusion[1,1],
                  FP=mdl$confusion[1,2],
                  FN=mdl$confusion[2,1] )

write.table(df,file= finfo,row.names=F,sep=",",quote=F)

