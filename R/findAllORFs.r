#' Find all variants of open reading frames in nucleotide sequence(-s)
#' @description Identify all possible variants of open reading frames in
#'     nucleotide sequence(-s) of interest.
#' @param x nucleotide sequence(-s) of interest which can be specified in the
#'     following ways (see examples below):
#'     i)   as character string of single nucleotide sequence;
#'     ii)  as character vector of several nucleotide sequences;
#'     iii) as R system file "Set.trans_sequences.fasta" with test set of
#'          nucleotide sequences;
#'     iv)  as character string giving the name of file with custom sequences
#'          of interest. Usually it is a set of sequences of RNA molecules.
#'          This file must be located in working directory. Allowed file formats
#'          are "fasta" or "fa";
#'     v)   as ready-to-use object of class DNAStringSet stored in the
#'          computer's RAM.
#' @param codStart character string with type of start codon (according to
#'     Chothani S. P. et al. Mol Cell. 2022 Aug 4;82(15):2885-2899.e8): "ATG",
#'     "CTG", "GTG", "AGG", "TTG", "AAG", "ACG", "ATC", "ATT" or "ATA".
#'     Default value is "ATG".
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return data frame with six following fields:
#'     i) transcript_id - ID of sequence;
#'     ii) orf_id       - ID of open reading frame;
#'     iii) start       - start coordinate of open reading frame in a sequence;
#'     iv) end          - end coordinate of open reading frame in a sequence;
#'     v) length        - length of open reading frame;
#'     vi) orf.sequence - sequence of open reading frame.
#' @author Vasily V. Grinev
# Last updated: June 12, 2024.
#' @examples
#' ### Input of nucleotide sequence of interest as character string:
#' tr <- "ACGCCTAGATGGTAGGGCCCATTTGCGAAGGACGAAGTGCCCAGTGAGTTAGCCAC"
#' orf <- findAllORFs(x=tr, codStart="ATG", workDir=NULL)
#' ### Input of several nucleotide sequences of interest as character vector:
#' tr <- c("ACGCCTAGATGGTAGGGCCCATTTGCGAAGGACGAAGTGCCCAGTGAGTTAGCCAC",
#'         "GGGCCATGAGGGTAGCGCACATGGATGGACGATGTGCCCAGTGGGTTAGCCAC")
#' orf <- findAllORFs(x=tr, codStart="ATG", workDir=NULL)
#' ### Input of R system file with test set of nucleotide sequences:
#' tr <- system.file("extdata",
#'                   "Set.trans_sequences.fasta",
#'                   package="ORFhunteR")
#' orf <- findAllORFs(x=tr, codStart="ATG", workDir=NULL)
#' ### Input of custom file with sequences of interest:
#' # tr <- "path/Set.trans_sequences.fasta" 
#' # orf <- findAllORFs(x=tr, codStart="ATG", workDir=NULL)
#' ### Usage of alternative start codon:
#' tr <- "ACGCCTAGCTGGATTTGCCCATATGCGATTTACGAAGTGCCCAGTGAGTTAGCCAC"
#' orf <- findAllORFs(x=tr, codStart="CTG", workDir=NULL)
#' @export

findAllORFs <- function(x, codStart="ATG", workDir=NULL){
  ### Loading of the nucleotide sequence(-s) of interest as an object of
  #   class DNAStringSet.
  if (class(x=x) == "character" & unique(x=tools::file_ext(x=x) == "")){
    seqs <- DNAStringSet(x=x)
    names(x=seqs) <- paste("seq",
                           formatC(x=seq(from=1,
                                         to=length(x=seqs)),
                                   width=nchar(x=length(x=seqs)),
                                   flag="0"),
                           sep="")
  }
  if (class(x=x) == "character" & unique(x=gsub(pattern=".+extdata.+",
                                                replacement="extdata",
                                                x=x) == "extdata")){
    seqs <- readDNAStringSet(filepath=x)
    seqs <- seqs[order(x=names(x=seqs)), ]
  }
  if (class(x=x) == "character" & unique(x=tools::file_ext(x=x)
                                         %in% c("fa", "fasta"))){
    if (!is.null(workDir)){
     x <- paste(workDir, x, sep="/")
    }
    seqs <- readDNAStringSet(filepath=x)
    seqs <- seqs[order(x=names(x=seqs)), ]
  }
  if (class(x=x) == "DNAStringSet"){
    seqs <- x
    seqs <- seqs[order(x=names(x=seqs)), ]
  }
  ### Identification of all possible locations of start and stop codons.
  tr_start <- unlist(x=vmatchPattern(pattern=codStart, subject=seqs))
  tr_start <- data.frame(tr_start)
  tr_start <- tr_start[, c(4, 1:2)]
  colnames(x=tr_start) <- c("seqnames", "start", "end")
  tr_stop <- unlist(x=c(vmatchPattern(pattern="TAA", subject=seqs),
                        vmatchPattern(pattern="TAG", subject=seqs),
                        vmatchPattern(pattern="TGA", subject=seqs)))
  tr_stop <- data.frame(tr_stop)
  tr_stop <- tr_stop[, c(4, 1:2)]
  colnames(x=tr_stop) <- c("seqnames", "start", "end")
  if(nrow(tr_start) == 0 | nrow(tr_stop) == 0){
    message("No start of stop codons found \n")
    return(NULL)
  }
  tr_start <- tr_start[tr_start$seqnames %in% tr_stop$seqnames, ]
  tr_start <- tr_start[order(x=tr_start$seqnames), ]
  rownames(x=tr_start) <- NULL
  tr_stop <- tr_stop[tr_stop$seqnames %in% tr_start$seqnames, ]
  tr_stop <- tr_stop[order(x=tr_stop$seqnames), ]
  rownames(x=tr_stop) <- NULL
  tr_start <- split(x=tr_start, f=tr_start$seqnames)
  
  tr_stop <- split(x=tr_stop, f=tr_stop$seqnames)
  seqs <- seqs[names(x=seqs) %in% names(x=tr_start)]
  tr_start <- tr_start[names(x=tr_start) %in% names(x=seqs)]
  tr_stop <- tr_stop[names(x=tr_stop) %in% names(x=seqs)]
  ### Identification of all possible variants of open reading frames.
  ORFs <- list()
  for (i in seq.int(1,length(tr_start))){
    inFrame <- outer(X=tr_stop[[i]]$start,
                     Y=tr_start[[i]]$start,
                     FUN="-")/3
    inFrame <- which(x=round(x=inFrame) == inFrame & inFrame > 0,
                     arr.ind=TRUE)
    if (nrow(x=inFrame) > 0){
      inFrame[, 1] <- tr_stop[[i]]$start[inFrame[, 1]]
      inFrame[, 2] <- tr_start[[i]]$start[inFrame[, 2]]
      inFrame <- aggregate(x=inFrame[, 1],
                           by=list(inFrame[, 2]),
                           FUN=min)
      inFrame <- aggregate(x=inFrame[, 1],
                           by=list(inFrame[, 2]),
                           FUN=min)
      inFrame <- inFrame[, c(2, 1)]
      inFrame <- inFrame[order(inFrame[, 1]), ]
      if (nrow(x=inFrame) > 0){
        inFrame[, 2] <- inFrame[, 2] + 2
        inFrame <- cbind(names(x=tr_stop[i]),
                         inFrame,
                         substring(text=as.character(x=seqs[[i]]),
                                   first=inFrame[, "x"],
                                   last=inFrame[, "Group.1"]))
        colnames(x=inFrame) <- c("seqnames",
                                 "start", "end",
                                 "orf.sequence")
        rownames(x=inFrame) <- NULL
        ORFs[[i]] <- inFrame
      }
    }
  }
  ORFs <- do.call(what=rbind, args=ORFs)
  orf_id <- unlist(x=sapply(X=table(ORFs$seqnames),
                            FUN=function(y){paste(".ORF",
                                                  formatC(x=seq(from=1,
                                                                to=y),
                                                          width=3,
                                                          flag="0"),
                                                  sep="")}))
  ORFs$orf_id <- paste(ORFs$seqnames, as.vector(x=orf_id), sep="")
  ORFs$length <- ORFs$end - ORFs$start + 1
  ORFs <- ORFs[, c(1, 5, 2:3, 6, 4)]
  colnames(ORFs) <- c("transcript_id", "orf_id",
                      "start", "end", "length",
                      "orf.sequence")
  ### Returning the final object.
  return(ORFs)
}
