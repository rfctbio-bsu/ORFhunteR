### Loading of required packages.
##  Listing of required packages.
reqPkgs <- c("Biostrings", "BSgenome.Hsapiens.UCSC.hg38", "data.table",
             "parallel", "Peptides", "Rcpp", "randomForest", "rtracklayer",
             "stringr")
##  Checking for installed required packages.
insPkgs <- reqPkgs[reqPkgs %in% installed.packages()[, 1]]
if (length(x=reqPkgs) > length(x=insPkgs)){
    if (length(x=reqPkgs) - length(x=insPkgs) == 1){
        stop(paste(paste("The following package has been not installed: ",
                         toString(x=paste0(reqPkgs[!reqPkgs %in% insPkgs])),
                         ".",
                         sep=""),
                   "Please, install this package before analysis."))
    }else{
        stop(paste(paste("The following packages have been not installed: ",
                         toString(x=paste0(reqPkgs[!reqPkgs %in% insPkgs])),
                         ".",
                         sep=""),
                   "Please, install these packages before analysis."))
    }
}else{
    print("All required packages are installed.")
}
##  Loading of required packages.
suppressMessages(expr=library(package=Biostrings))
suppressMessages(expr=library(package=BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(expr=library(package=data.table))
suppressMessages(expr=library(package=parallel))
suppressMessages(expr=library(package=Peptides))
suppressMessages(expr=library(package=Rcpp))
suppressMessages(expr=library(package=randomForest))
suppressMessages(expr=library(package=rtracklayer))
suppressMessages(expr=library(package=stringr))
#suppressMessages(expr=library(package=genome, character.only=TRUE))
##  Setting of C++ environment.
setwd("D:/Vasily Grinev")
Rcpp::sourceCpp("getBaoMetrics.cpp")
Rcpp::sourceCpp("getCorrelationFactors.cpp")


true_orfs = read.table(file = "true_orfs.txt",
                  sep = "\t",
                  header = TRUE,
                  quote = "\"",
                  as.is = TRUE)
rownames(anno) <- anno[, 1]
anno <- as.matrix(anno)
input <- as.matrix(prob_orfs)
input <- cbind(input, "", "")
index = input[, 1] %in% anno[, 1]
input[index, 6:7] <- anno[input[index, 1], 2:3]
input <- data.frame(input)
input[, 2:7] <- apply(input[, 2:7], 2, as.numeric)
input$V6 <- input$V6 - input$start
input$V7 <- input$V7 - input$end
input
