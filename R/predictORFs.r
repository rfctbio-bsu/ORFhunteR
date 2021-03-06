#' Predict the true ORFs in mRNA molecules
#' @description Predict the true ORFs in mRNA molecules.
#' @param tr character string giving the name of file with experimental
#'     transcripts. Allowed file formats are "fasta", "fa", "gtf" or "gff".
#' @param genome character string giving the name of BSgenome data package with
#'     full genome sequences. Default value is "BSgenome.Hsapiens.UCSC.hg38".
#' @param prThr probability threshold for the "winning" class of ORFs.
#'     Default value is 0.5.
#' @param model character string giving the connection or full path to the file
#'     from which the classification model is read. Use default NULL value to 
#'     use our default model.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return The coordinates of ORFs in mRNA molecules of interest.
#' @author Mikalai M. Yatskou
#' @examples
#' \dontrun{
#' tr_path <- system.file("extdata",
#'                        "Set.trans_sequences.fasta",
#'                        package="ORFhunteR")
#' model <- "http://www.sstcenter.com/download/ORFhunteR/classRFmodel_1.rds"
#' ORFs <- predictORF(tr=tr_path, model=model)
#' }
#' @export

predictORF <- function(tr,
                       genome="BSgenome.Hsapiens.UCSC.hg38",
                       model=NULL,
                       prThr=0.5,
                       workDir=NULL){
    ### Loading of the experimental transcripts as a list of character strings.
    exp_trans <- loadTrExper(tr=tr, genome=genome, workDir=workDir)
    ### Identification of all possible ORFs in experimental transcripts.
    cl <- makeCluster(spec=detectCores() - 1)
    clusterExport(cl=cl,
                varlist=c("findORFs", "DNAStringSet", "matchPDict", "DNAString",
                        "aggregate", "codonStartStop", "c", "outer", "round",
                        "length", "t", "lower.tri", "cbind", "list", "min",
                        "order", "substring", "do.call", "colnames", "return",
                        "names", "as.character", "sort", "unlist", "start"))
    all_orfs <- parLapply(X=exp_trans,
                          fun=function(y){orf <- findORFs(x=y)},
                          cl=cl)
    stopCluster(cl=cl)
    all_orfs <- cbind(transcript_id=rep(x=names(x=exp_trans),
                                        times=lengths(all_orfs)/4),
                      do.call(what=rbind, args=all_orfs))
    all_orfs <- all_orfs[as.numeric(all_orfs[, "length"]) >= 42, ]
    ### Classification of the all identified ORFs.
    ##  Annotation of the all identified ORFs with sequence features.
    prob_orfs <- vectorizeORFs(x=DNAStringSet(x=all_orfs[, "orf.sequence"]))
    ##  Import of the classification model.
    if(is.null(x=model)){
        model <- download_model_file()
    } else {
        model <- readRDS(file=file(description=model))
    }
    ##  Classification.
    prob_orfs <- predict(object=model, newdata=prob_orfs, type="prob")
    prob_orfs <- data.table::data.table(cbind(all_orfs[, c("transcript_id",
                                               "start",
                                               "end",
                                               "length")],
                                  prob=prob_orfs[, "3"]))
    ### Data aggregation anf filtration.
    prob <- NULL
    true_orfs <- prob_orfs[prob_orfs[, .I[prob == max(prob)],
                                   by="transcript_id"]$V1]
    true_orfs <- true_orfs[true_orfs[, .I[length == max(length)],
                                     by="transcript_id"]$V1]
    true_orfs <- true_orfs[true_orfs[, .I[start == min(start)],
                                     by="transcript_id"]$V1]
    true_orfs <- data.frame(true_orfs)
    if (length(x=exp_trans) > length(x=true_orfs$transcript_id)){
        mis <- names(x=exp_trans)[!names(x=exp_trans) %in% 
                                      true_orfs$transcript_id]
        mat <- matrix(0L, nrow=length(x=mis), ncol=4)
        colnames(x=mat) <- c("start", "end", "length", "prob")
        true_orfs <- rbind(true_orfs, cbind(transcript_id=mis, mat))
    }
    true_orfs <- true_orfs[order(x=true_orfs$transcript_id), ]
    true_orfs[, c("start", "end", "length", "prob")] <-
        apply(X=true_orfs[, c("start", "end", "length", "prob")],
              MARGIN=2,
              FUN=as.numeric)
    true_orfs <- true_orfs[true_orfs[, "prob"] >= prThr, ]
    rownames(x=true_orfs) <- NULL
    ### Returning a final object of class data.frame.
    return(true_orfs)
}

download_model_file <- function() {
    path <- system.file("extdata", "cl_model.rds", package = "ORFhunteR")
    model <- readRDS(path)
    if(length(model) < 1){
        model <- readRDS(url("http://www.sstcenter.com/download/ORFhunteR/classRFmodel_1.rds", "r"))
        saveRDS(model, path)
    }
    return(model)
}
