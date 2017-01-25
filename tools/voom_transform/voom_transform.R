#!/usr/bin/env Rscript

options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
library("limma")

spec <- matrix(c( 
    "expressionmatrix", "e", 1, "character",
    "transformedmatrix", "t", 1, "character"
    ), byrow=TRUE, ncol=4)
opt <- getopt(spec)

em <-  read.delim(opt$expressionmatrix,header=T,row.names=1,stringsAsFactors=F,check.names=FALSE,na.strings=c(""))

vem <- voom(em)
write.table(file=opt$transformedmatrix,vem$E,sep="\t",row.names=TRUE,col.names=NA,quote=F)
