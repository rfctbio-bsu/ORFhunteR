#' Translate any nucleotide sequences to peptide sequences
#' @description Translation of any nucleotide sequences to peptide sequences.
#' @param x character string giving the path to FASTA file with nucleotide
#'     sequences of interest.
#' @param codones location of tab-delimited TXT file with genetic code. Default
#'     value is "internal" for system extdata. Alternatively, the name of
#'     tab-delimited TXT file located in work directory should be specified.
#' @param aaSymbol type of amino acid symbols (one- or three-letter coding).
#'     Default value is 1.
#' @param pep_length minimal length of in silico translated peptide/protein to
#'     be saved. Default value is 8.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return list of three AAStringSet objects with peptide sequences.
#' @author Vasily V. Grinev, Katsiaryna V. Huzava
#' Last updated: August 19, 2024.
#' @examples
#' seqs <- system.file("extdata",
#'                    "Set.trans_sequences.for.translation.fasta",
#'                    package="ORFhunteR")
#' res <- translate(x=seqs,
#'                  codones="internal",
#'                  aaSymbol=1,
#'                  pep_length=8,
#'                  workDir=NULL)
#' @export

translate <- function(x,
                      codones="internal",
                      aaSymbol=1, 
                      pep_length=8, 
                      workDir=NULL){
  ### Loading of the nucleotide sequence(-s) of interest as an object of
  #   class DNAStringSet.
  if (!is.null(x=workDir)){
    seqs <- paste(workDir, x, sep="/")
  }
  ##  Retrieving the file extension.
  frt <- tools::file_ext(x=x)
  ##  Validation of file format.
  if (!frt %in% c("fa", "fasta")){
    message("Invalid file format")
    return(NULL)
  }
  ##  Full path to the file.
  ##  Loading of the sequence(-s) of interest.
  seqs0 <- readDNAStringSet(filepath=seqs)
  seqs1 <- subseq(x=seqs0, start=2, end=width(x=seqs0))
  seqs2 <- subseq(x=seqs0, start=3, end=width(x=seqs0))
  ### Calling the codon table.
  if (codones == "internal"){
    codon_table_path <- system.file("extdata",
                                    "codon.table.txt",
                                    package="ORFhunteR")
  }else{
    if (!is.null(x=workDir)){
      codon_table_path <- paste(workDir, codones, sep = "/")
    } else{
      codon_table_path <- codones
    }
  }
  codon_table <- read.table(file=codon_table_path,
                            sep="\t",
                            header=TRUE,
                            quote="\"",
                            as.is=TRUE)
  ### Nucleotide sequences to proteins translation.
  aa0 <- lapply(X=seqs0,
                FUN=function(y){y <- as.character(x=y)
                prt <- sapply(X=seq(from=1,
                                    to=nchar(x=y) - 2,
                                    by=3),
                              FUN=function(z){substr(x=y,
                                                     start=z,
                                                     stop=z + 2)
                              }
                )
                prt <- paste(codon_table[match(x=prt,
                                               table=codon_table[, 1]),
                                         if (aaSymbol == 1){
                                           4
                                         }else{
                                           3
                                         }],
                             collapse="")
                }
  )
  aa0 <- AAStringSet(x=unlist(x=aa0), use.names=TRUE)
  aa0 <- AAStringSet(x=gsub(pattern="\\*.*", replacement="", x=aa0))
  aa0 <- aa0[width(x=aa0) >= pep_length, ]
  aa1 <- lapply(X=seqs1,
                FUN=function(y){y <- as.character(x=y)
                prt <- sapply(X=seq(from=1,
                                    to=nchar(x=y) - 2,
                                    by=3),
                              FUN=function(z){substr(x=y,
                                                     start=z,
                                                     stop=z + 2)
                              }
                )
                prt <- paste(codon_table[match(x=prt,
                                               table=codon_table[, 1]),
                                         if (aaSymbol == 1){
                                           4
                                         }else{
                                           3
                                         }],
                             collapse="")
                }
  )
  aa1 <- AAStringSet(x=unlist(x=aa1), use.names=TRUE)
  aa1 <- AAStringSet(x=gsub(pattern="\\*.*", replacement="", x=aa1))
  aa1 <- aa1[width(x=aa1) >= pep_length, ]
  aa2 <- lapply(X=seqs2,
                FUN=function(y){y <- as.character(x=y)
                prt <- sapply(X=seq(from=1,
                                    to=nchar(x=y) - 2,
                                    by=3),
                              FUN=function(z){substr(x=y,
                                                     start=z,
                                                     stop=z + 2)
                              }
                )
                prt <- paste(codon_table[match(x=prt,
                                               table=codon_table[, 1]),
                                         if (aaSymbol == 1){
                                           4
                                         }else{
                                           3
                                         }],
                             collapse="")
                }
  )
  aa2 <- AAStringSet(x=unlist(x=aa2), use.names=TRUE)
  aa2 <- AAStringSet(x=gsub(pattern="\\*.*", replacement="", x=aa2))
  aa2 <- aa2[width(x=aa2) >= pep_length, ]
  ### Returning a final object.
  return(list(aa0=aa0, aa1=aa1, aa2=aa2))
}
