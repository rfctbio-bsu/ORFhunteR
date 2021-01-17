#' Classify the pseudo and true ORF candidates derived from RNA molecules
#' @description Clussify the pseudo and true ORF candidates.
#' @param ORFLncRNAs the oject of the class list of non-coding pseudo ORFs.
#' @param ORFmRNAs the oject of the class list of coding true ORFs.
#' @param pLearn probability threshold for the "traing" class of ORFs.
#'     Default value is 0.75.
#' @param nTrees nuber of the decision trees in randomForest.
#'     Default value is 500.
#' @param modelRF character string giving the name of RDS-file to store
#'     the classification model. NULL by default.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @param showAccuracy logic TRUE or FALSE. Use TRUE for print accuracy metrics. 
#'     Default value is FALSE
#' @return The classificator object of the class randomForest.
#' @author Mikalai M. Yatskou
#' @examples
#' \dontrun{
#' clt <- classifyORFsCandidates(ORFLncRNAs, ORFmRNAs)
#' }
#' @export
#' @import Biostrings BSgenome.Hsapiens.UCSC.hg38 Peptides
#' @import data.table stringr randomForest rtracklayer
#' @import xfun parallel
#' @importFrom Rcpp sourceCpp
#' @importFrom stats aggregate predict 
#' @importFrom utils read.table
#' @importFrom graphics plot
#' @useDynLib ORFhunteR

classifyORFsCandidates <- function(ORFLncRNAs, 
                                   ORFmRNAs,
                                   pLearn = 0.75, 
                                   nTrees = 500,
                                   modelRF = NULL,
                                   workDir = NULL,
                                   showAccuracy = FALSE){
    ### Vectorizing ORFs.
    ORFLncRNAs <- DNAStringSet(x = do.call(what=rbind, args=ORFLncRNAs))
    ORFLncRNAsVctd <- vectorizeORFs(ORFLncRNAs)
    ORFmRNAs <- DNAStringSet(do.call(rbind, ORFmRNAs))
    ORFmRNAsVctd <- vectorizeORFs(x = ORFmRNAs)
    
    #****************************************
    ### Making the ORF dataset for classification.
    NL <- dim(x=ORFLncRNAsVctd)[1]
    Nm <- dim(x=ORFmRNAsVctd)[1]
    dataset <- rbind(ORFLncRNAsVctd, ORFmRNAsVctd)
    Labels <- c(rep(x=2, times=NL), rep(x=3, times=Nm))
    dataset <- cbind(dataset, Labels)
    
    #****************************************
    ### Calssification of ORFs by randomForest.
    D <- dim(x=dataset)
    NRNAs <- D[1]
    NFtrs <- D[2]
    ##  Forming the learning and testing datasets.
    LearnSize <- pLearn * D[1]
    # set.seed(1)
    IndexesLearn <- sample(x=1:D[1], size=LearnSize)
    dataLearn <- dataset[IndexesLearn, ]
    dataTestCls <- dataset[-IndexesLearn, ]
    dataTest <- dataTestCls[, -NFtrs]
    dataTemp <- dataLearn
    dataTemp$Labels <- NULL
    clt <- randomForest(dataTemp, factor(dataLearn$Labels))
    ##  Learning dataset.
    e <- clt$confusion
    TP <- e[1, 1] #True positive
    FN <- e[2, 1] #False negative
    FP <- e[1, 2] #False positive
    TN <- e[2, 2] #True negative
    A <- (TP + TN)/(TP + TN + FP + FN)
    P <- TP/(TP + FP)
    R <- TP/(TP + FN)
    F1 <- (2 * TP)/(2 * TP + FP + FN) 
    if(showAccuracy){
        message("Training dataset\n")
        message("Accuracy:", A, "\nPrecision rate:", P,
            "\nRecall:", R, "\nThe score F1:", F1, "\n")
    }
    ##  Test dataset.
    e <- table(predict(object = clt, newdata = dataTest), 
               factor(x = dataTestCls$Labels))
    TP <- e[1, 1] 
    FN <- e[2, 1]
    FP <- e[1, 2]
    TN <- e[2, 2]
    A <- (TP + TN)/(TP+ TP + FP + FN)
    P <- TP/(TP + FP)
    R <- TP/(TP + FN)
    F1 <- (2 * TP)/(2 * TP + FP + FN) 
    if(showAccuracy){
        message("Test dataset\n")
        message("Accuracy:", A, "\nPrecision rate:", P,
            "\nRecall:", R, "\nThe score F1:", F1, "\n")
    }
    #****************************************
    ### Saving the classification model into the RDS-file.
    ##  Full path to the file.
    if (!is.null(x = modelRF)){
        if (!is.null(x=workDir)){
            saveRDS(object=clt, file=paste(workDir, modelRF, sep="/"))
        } else{
            saveRDS(object=clt, file=modelRF)
        }
    }
    ### Returning a final object of class randomForest.
    return(clt)
}
