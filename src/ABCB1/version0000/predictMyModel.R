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
#pred <- as.numeric(mdl$classes[predict(mdl,toPredict)])

pred <- predict(mdl,toPredict,type="prob")[,"1"]
pred <- ifelse( pred < mdl$threshold, 0,1)


# save result
write.table(pred,file= tmpfile ,row.names=F,sep=",",quote=F,col.names=F)
