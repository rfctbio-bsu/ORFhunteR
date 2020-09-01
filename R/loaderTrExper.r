#' Load a set of transcripts
#' @description Load a set of experimental transcripts.
#' @param f.TrExper character string giving the name of file with experimental transcripts. Allowed
#'               file formats are "fasta", "fa", "gtf" or "gff".
#' @param genome character string giving the name of BSgenome data package
#'               with full genome sequences. Default value is "BSgenome.Hsapiens.UCSC.hg38".
#' @param d.work character string giving the path to and name of work directory. Concatenated with filename via "/". Default value is NULL.
#' @return list of loaded transcript sequences.
#' @author Vasily V. Grinev
#' @examples
#' f.TrExper_path <- system.file("extdata", "Set.trans_sequences.gtf", package = "ORFhunteR")
#' tr.seq = loaderTrExper(f.TrExper = f.TrExper_path,
#'                        genome = "BSgenome.Hsapiens.UCSC.hg38",
#'                        d.work = NULL)
#' @export

loaderTrExper <- function(f.TrExper, genome = "BSgenome.Hsapiens.UCSC.hg38", 
                          d.work = NULL){
  ### Validation of file format.
  ##  Loading of required auxiliary library tools.
  #suppressMessages(expr = library(package = tools))
  ##  Retrieving the file extension.
  frt <- tools::file_ext(x = f.TrExper)
  ##  Validation of file format.
  if (!frt %in% c("fasta", "fa", "gtf", "gff")){
    stop("Invalid file format")
  }
  ### Loading of experimental transcripts in FASTA/FA format.
  if(!is.null(d.work)){
    f.TrExper <- paste(d.work, f.TrExper, sep = "/")
  }
  if (frt == "fasta" || frt == "fa"){
    ##  Loading of required auxiliary library Biostrings.
    #suppressMessages(expr = library(package = Biostrings))
    ##  Loading of experimental transcripts as a list of character strings.
    seq.set <- as.list(x = as.character
                       (x = readDNAStringSet(filepath = f.TrExper,
                                                            format = frt)))
  }
  ### Loading of experimental transcripts in GTF/GFF format.
  if (frt == "gtf" || frt == "gff"){
    ##  Loading of required auxiliary libraries.
    #suppressMessages(expr = library(package = rtracklayer))
    #suppressMessages(expr = library(package = genome, character.only = TRUE))
    ##  Loading of experimental transcripts as a list of character strings.
    seq.set <- import(con = f.TrExper, format = frt)
    seq.set <- seq.set[seq.set$type == "exon", ]
    seq.set <- split(x = seq.set, f = seq.set$transcript_id)
    seq.set <- lapply(X = seq.set,
                     FUN = function(y){as.character(x = unlist(
                       x = getSeq(x = get(x = genome), names = y)))})
  }
  return(seq.set)
}
