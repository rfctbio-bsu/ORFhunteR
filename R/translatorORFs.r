#' Translate open reading frames to proteins
#' @description Translate the identified open reading frames to proteins.
#' @param f.seq_orf character string giving the name of FASTA file with sequences of identified open
#'              reading frames.
#' @param aa.symbol type of amino acid symbols (one- or three-letter coding). Default value is 1.
#' @param f.prts character string giving the name of output FASTA file with sequences of in sillico
#'               translated proteins. Default value is NULL.
#' @param d.work character string giving the path to and name of work directory. Default value is NULL.
#' @return AAStringSet object with sequences of obtained proteins.
#' @author Vasily V. Grinev
#' @examples
#' f.seq_orf_path <- system.file("extdata", "Set.trans_ORFs.sequences.fasta", package = "ORFhunteR")
#' p.seq <- translatorORFs(f.seq_orf = f.seq_orf_path,
#'                         aa.symbol = 1)
#' @export

translatorORFs <- function(f.seq_orf, aa.symbol = 1, f.prts = NULL, 
                           d.work = NULL){
  ### Loading of required auxiliary library.
  #   This code was successfully tested with library Biostrings v.2.52.0.
  #suppressMessages(expr = library(package = Biostrings))
  ### Internal calling of codon table.
  #codon.table = read.table(file = system.file("extdata", "codon.table.txt", 
  #package = "ORFhunteR"),
  codon.table_path <- system.file("extdata", "codon.table.txt", 
                                  package = "ORFhunteR")
  codon.table <- read.table(file = codon.table_path,
                           sep = "\t",
                           header = TRUE,
                           quote = "\"",
                           as.is = TRUE)
  ### Loading of sequences of open reading frames as an object of class 
  #DNAStringSet.
  if(!is.null(d.work)){
    f.seq_orf <- paste(d.work, f.seq_orf, sep = "/")
    f.prts <- paste(d.work, f.prts, sep = "/")
  }
  seq.orf <- readDNAStringSet(filepath = f.seq_orf)

  ### Open reading frames to proteins translation.
  seq.prts <- lapply(X = seq.orf,
                    FUN = function(y){y = as.character(x = y)
                    seq.protein = sapply(X = seq(from = 1,
                                                 to = nchar(x = y) - 3,
                                                 by = 3),
                                         FUN = function(z){substr(x = y,
                                                                start = z,
                                                                stop = z + 2)})
                    seq.protein = as.vector(x = seq.protein)
                    seq.protein = paste(codon.table[match(x = seq.protein,
                                                      table = codon.table[, 1]),
                                                if (aa.symbol == 1){4}else{3}],
                                        collapse = "")})
  seq.prts <- AAStringSet(x = unlist(x = seq.prts), use.names = TRUE)
  ### Write an XStringSet object to a file.
  if (!is.null(f.prts)){
    writeXStringSet(x = seq.prts, filepath = f.prts)
  }
  ### Returning a final object.
  return(seq.prts)
}
