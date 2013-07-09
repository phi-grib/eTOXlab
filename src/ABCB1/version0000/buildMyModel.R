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
#mdl <- tuneRF(X, Y ,  plot=F ,doBest=TRUE , 
#            improve = 0.05 ,  stepFactor=2, ntreeTry=50)
            
# alternatively build a single Random Forest model
mdl  <- randomForest(X,Y,ntree=500)

getSensSpec <- function(thr, pred, actual){
    pred <- sapply( pred , function(x){ ifelse(x<thr, 0,1)})
    tn = sum( (pred == actual ) & actual == 0 )
    tp = sum( (pred == actual ) & actual == 1 )  
    fp = sum( (pred != actual ) & actual == 0 )
    fn = sum( (pred != actual ) & actual == 1 )
    spec = tn / (tn+fp)
    sens = tp / (tp+fn)
    abs(spec-sens)
}
prob <- predict(mdl,X,type="prob")

ths <- seq(0,1,length.out=30)
nY <- as.numeric(levels(Y)[Y])

sesp <- sapply( ths, getSensSpec, prob[,"1"], nY)

mdl$threshold <- ths[which(min(sesp)==sesp)[1]]

# save model
save(mdl, file= fRdata )

pred <- sapply( prob[,"1"] , function(x){ ifelse(x<mdl$threshold, 0,1)})
tn = sum( (pred == nY ) & nY == 0 )
tp = sum( (pred == nY ) & nY == 1 )  
fp = sum( (pred != nY ) & nY == 0 )
fn = sum( (pred != nY ) & nY == 1 )
    
df <- data.frame( TP = tp, TN = tn, FP = fp, FN = fn )

write.table(df,file= finfo,row.names=F,sep=",",quote=F)

