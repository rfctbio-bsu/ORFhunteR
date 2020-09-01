#' Filter of human Ensembl transcripts with inconsistent open reading frames
#' @description Filter of human Ensembl transcripts with inconsistent open reading frames.
#' @param f.fasta character string giving the name of multi-FASTA file.
#' @param f.mtRNAs character string giving the name of TXT file in tab-delimited format containing
#'                  the list of mitochondrial transcripts. This file should include one mandatory
#'                  field transcript_id.
#' @param f.orf_coords character string giving the name of TXT file in tab-delimited format containing
#'                  coordinates of open reading frames for transcripts to be analyzed.
#'                  This file should include three mandatory fields:
#'                  i) transcript_id - transcript ID;
#'                  ii) start        - start coordinate of ORF;
#'                  iii) end         - end coordinate of ORF.
#' @param d.work character string giving the path to and name of work directory. Concatenated with filename via "/". Default value is NULL.
#' @return DNAStringSet object
#' @author Vasily V. Grinev
#' @examples
#' f.fasta_path <- system.file("extdata", "Ensembl_release_85_GRCh38.p7_mRNAs_ex.fasta",
#'                            package = "ORFhunteR")
#' f.mtRNAs_path <- system.file("extdata", "Ensembl_release_85_GRCh38.p7_mitochondrial_RNAs.txt",
#'                            package = "ORFhunteR")
#' f.orf_coords_path <- system.file("extdata",
#'                            "Ensembl_release_85_GRCh38.p7_CDS_coordinates_of_mRNAs.txt",
#'                            package = "ORFhunteR")
#' res <- filterEnsemblORFs(f.fasta = f.fasta_path,
#'                          f.mtRNAs = f.mtRNAs_path,
#'                          f.orf_coords = f.orf_coords_path,
#'                          d.work = NULL)
#' # writeXStringSet(x = res, filepath = "filtered_mRNAs.fasta")
#' @export

filterEnsemblORFs <- function(f.fasta, f.mtRNAs, f.orf_coords, d.work = NULL){
  ##  Loading of required auxiliary library.
  #   This code was successfully tested with library Biostrings v.2.52.0.
  # suppressMessages(expr = library(package = Biostrings))
  ##  Loading of FASTA file in an XStringSet object.
  if(!is.null(d.work)){
    f.fasta <- paste(d.work, f.fasta, sep = "/")
    f.mtRNAs <- paste(d.work, f.mtRNAs, sep = "/")
    f.orf_coords <- paste(d.work, f.orf_coords, sep = "/")
  }
  tr <- readDNAStringSet(filepath = f.fasta,
                         format = "fasta",
                         use.names = TRUE)
  ##  Removing of mitochondrial transcripts.
  m.tr <- read.table(file = f.mtRNAs,
                    sep = "\t",
                    header = TRUE,
                    quote = "\"",
                    as.is = TRUE)$transcript_id
  tr <- tr[!names(x = tr) %in% m.tr, ]
  ##  Removing transcripts lacking a canonical start and/or 
  #stop codons inside the sequence.
  orf_coords <- read.table(file = f.orf_coords,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          as.is = TRUE)
  orf_coords <- orf_coords[!orf_coords$transcript_id %in% m.tr, ]
  orf_coords <- orf_coords[order(orf_coords$transcript_id), ]
  tr <- tr[names(x = tr) %in% orf_coords$transcript_id, ]
  tr <- tr[order(names(x = tr)), ]
  codon.start <- substr(x = tr, start = orf_coords$start, 
                        stop = orf_coords$start + 2)
  codon.start <- attr(x = codon.start[codon.start == "ATG"], which = "names")
  codon.stop <- substr(x = tr, start = orf_coords$end - 2, 
                       stop = orf_coords$end)
  codon.stop <- attr(x = codon.stop[codon.stop %in% c("TAA", "TAG", "TGA")], 
                     which = "names")
  tr <- tr[names(x = tr) %in% codon.start[codon.start %in% codon.stop], ]
  ##  Returning the final object.
  return(tr)
}
