#' Predict the ORFs coordinates
#' @description Predict the ORFs coordinates and write it in txt file in your directory.
#' @param f.TrExper character string giving the name of the data file.
#' @param genome character string giving the name of BSgenome data package
#'               with full genome sequences. Default value is "BSgenome.Hsapiens.UCSC.hg38".
#' @param d.work character string giving the path to and name of work directory. Concatenated with filename via "/". Default value is NULL.
#' @author Mikalai M. Yatskou
#' @return The ORF coordinates of the RNA molecules.
#' @examples
#' \dontrun{
#' ORFs_coordinates <- predictORF(f.TrExper = "mRNAs_first200.fasta", d.work = NULL)
#' }
#' @export
#'
#' @import Biostrings BSgenome.Hsapiens.UCSC.hg38 Peptides
#' @import data.table stringr randomForest rtracklayer
#' @import xfun
#' @importFrom Rcpp sourceCpp
#' @importFrom stats aggregate predict
#' @importFrom utils read.table write.table
#' @useDynLib ORFhunteR


predictORF <- function(f.TrExper, genome = "BSgenome.Hsapiens.UCSC.hg38", d.work = NULL){

  t0 <- Sys.time()
  seq.set_all <- loaderTrExper(f.TrExper = f.TrExper,
                               genome = genome,
                               d.work = d.work)

  N_parce <- 5
  N_mol_total <- length(seq.set_all)
  Parced_Parts <-  N_mol_total/N_parce
  ORFs_coords_mRNAs <- NULL

  j_total <- 0
  N_mol_slctd <- 0

  for(l in 1 : N_parce){
    t1 <- Sys.time()
    seq.set <- seq.set_all[round((l-1)*Parced_Parts+1):round(l*Parced_Parts)]

    ##  Identification of possible ORFs.
    res = lapply(X = seq.set,
                 FUN = function(y){ orf = finderORFs(x = y, d.work = d.work) })

    t2 <- Sys.time()
    cat("Time 1 = ",t2-t1,'\n')

    t3 <- Sys.time()
    #clt_path <- system.file("extdata", "Classif_RF_ORFslncRNA_vs_ORFmRNAs.rds", package = "ORFhunteR")
    #clt <- readRDS(clt_path)
    clt <- readRDS(url("http://www.sstcenter.com/download/ORFhunteR/Classif_RF_ORFslncRNA_vs_ORFmRNAs.rds", "r"))
    N_mol <- length(res)

    for (i in 1 : N_mol) {
      a <- res[[i]]
      x <- a[[4]]
      x_vctd <- Vectorising_RNAs_V(x)
      ORFs_lbls <- predict(clt,x_vctd, type="prob")
      ind <- which.max(as.vector(ORFs_lbls[,2]))
      Estd_start <- a[[1]][ind]
      Estd_stop <- a[[2]][ind]
      Estd_length <- a[[3]][ind]

      ORFs_coords = data.frame(Estd_start, Estd_stop, Estd_length)
      row.names(ORFs_coords) <- names(res[i])
      ORFs_coords_mRNAs <- rbind(ORFs_coords_mRNAs, ORFs_coords)
    }

    t4 <- Sys.time()
    cat("Time 2 = ",t4-t3,'\n')
    cat("mRNAs = ",N_mol,'\n')

    N_mol_slctd <- N_mol_slctd + N_mol
    cat("Total mRNAs = ",N_mol_slctd,'\n')
    cat("Iteration := ",l,'\n','\n')

  }

  print("********************" )
  cat("Total mRNAs = ",N_mol_slctd,'\n')
  cat("Total_time = ",t4-t0,'\n')

  ORFs_coords_mRNAs <- data.frame(row.names(ORFs_coords_mRNAs),
                                  ORFs_coords_mRNAs[,1],ORFs_coords_mRNAs[,2],
                                  ORFs_coords_mRNAs[,3])
  colnames(ORFs_coords_mRNAs) <- c('transcript_id','start','end', 'length')

  write.table(ORFs_coords_mRNAs,'ORF_predicted_coordinates.txt', 
              row.names = FALSE, quote = FALSE)
  return(ORFs_coords_mRNAs)
}
