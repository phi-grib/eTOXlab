# load packages
suppressPackageStartupMessages(library(randomForest))

# read arguments from command line
arg <- commandArgs( trailingOnly = T )
vpath   <- arg[1]
tmpfile <- arg[2]

# load model 
load( file= paste(vpath,"/modelR.Rdata",sep="") )

# load descriptors to predict
toPredict <- read.table(tmpfile ,sep=",",header=F)[,1]

# predict, get label, and get numeric value
pred <- as.numeric(mdl$classes[predict(mdl,toPredict)])

# save result
write.table(pred,file= tmpfile ,row.names=F,sep=",",quote=F,col.names=F)
