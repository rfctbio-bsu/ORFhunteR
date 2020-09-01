#' Annotate open reading frames
#' @description Annotate the identified open reading frames in nucleotide sequences of interest.
#' @param f.orf  character string giving the name of tab-delimited TXT file with coordinates of open reading frames. This file should include three mandatory fields:
#'             i) transcript_id - transcript ID;
#'             ii) start        - start coordinate of open reading frame in a transcript;
#'             iii) end         - end coordinate of open reading frame in a transcript.
#' @param f.fasta character string giving the name of FASTA file with nucleotide sequences of interest.
#' @param f.gtf character string giving the name of GTF/GFF file with transcripts of interest.
#'             Default value is NULL.
#' @param f.prts character string giving the name of FASTA file with sequences of in sillico translated proteins.
#' @param d.work character string giving the path to and name of work directory. Concatenated with filename via "/". Default value is NULL.
#' @return data.frame object with columns:
#'  \item{transcript_id}{transcript ID.}
#'  \item{f_utr.length}{length of 5'UTRs.}
#'  \item{start.codon}{type of start codon.}
#'  \item{orf.start}{start coordinate of open reading frames.}
#'  \item{orf.stop}{stop coordinate of open reading frames.}
#'  \item{stop.codon}{type of stop codon.}
#'  \item{stop.status}{PTC status of stop codon.}
#'  \item{orf.length}{length of open reading frames.}
#'  \item{t_utr.length}{length of 3'UTRs.}
#'  \item{MW}{molecular weight.}
#'  \item{pI}{isoelectic point of a protein sequence.}
#'  \item{indexPPI}{potential protein interaction index.}
#' @author Vasily V. Grinev
#' @examples
#' f.orf_path <- system.file("extdata", "Set.trans_ORFs.coordinates.txt", 
#'                           package = "ORFhunteR")
#' f.fasta_path <- system.file("extdata", "Set.trans_sequences.fasta", 
#'                           package = "ORFhunteR")
#' f.gtf_path <- system.file("extdata", "Set.trans_sequences.gtf", 
#'                           package = "ORFhunteR")
#' f.prts_path <- system.file("extdata", "Set.trans_proteins.sequences.fasta", 
#'                           package = "ORFhunteR")
#' anno = annotatorORFs(f.orf = f.orf_path,
#'                      f.fasta = f.fasta_path,
#'                      f.gtf = f.gtf_path,
#'                      f.prts = f.prts_path,
#'                      d.work = NULL)
#' @export

annotatorORFs <- function(f.orf, f.fasta, f.gtf = NULL, f.prts, d.work = NULL){
  ### Loading of required auxiliary library.
  #   This code was successfully tested with the libraries Biostrings v.2.52.0 
  #and Peptides v.2.4.1.
  #suppressMessages(expr = library(package = Biostrings))
  #suppressMessages(expr = library(package = Peptides))
  # Loading of open reading frames coordinates as an object of class data frame.
  if(!is.null(d.work)){
    f.orf <- paste(d.work, f.orf, sep = "/")
    f.fasta <- paste(d.work, f.fasta, sep = "/")
    f.prts <- paste(d.work, f.prts, sep = "/")
  }
  orf <- read.table(file = f.orf,
                   #sep = "\t",
                   header = TRUE,
                   quote = "\"",
                   as.is = TRUE)
  orf = orf[order(orf$transcript_id), ]
  ### Length of 5'UTRs.
  a.orf = orf[, 1:2]
  a.orf$start = a.orf$start - 1
  colnames(a.orf) = c("transcript_id", "f_utr.length")
  ### Type of start codon.
  ##  Loading of nucleotide sequences of interest as an object of class 
  #DNAStringSet.
  seq.set = readDNAStringSet(filepath = f.fasta)
  seq.set = seq.set[order(names(seq.set))]
  ##  Extraction of start codons.
  a.orf$start.codon = substr(x = seq.set, start = orf$start, 
                             stop = orf$start + 2)
  ### Coordinates of open reading frames.
  a.orf$orf.start = orf$start
  a.orf$orf.stop = orf$end - 2
  ### Type of stop codon.
  a.orf$stop.codon = substr(x = seq.set, start = a.orf$orf.stop, 
                            stop = a.orf$orf.stop + 2)
  ### PTC status of stop codon.
  if (!is.null(f.gtf)){
    ptc = finderPTC(f.orf = f.orf, f.trans = f.gtf, d.work = d.work)[, c(1, 5)]
    ptc = ptc[order(ptc$transcript_id), ]
    a.orf$stop.status = ptc$stop.status
    ptc = NULL
  }else{
    a.orf$stop.status = "ND"
  }
  ### Length of open reading frames.
  a.orf$orf.length = orf$end - orf$start + 1
  ### Length of 3'UTRs.
  a.orf$t_utr.length = width(seq.set) - orf$end
  ### Calculation the molecular weight of a protein sequence.
  ##  Loading of amino acid sequences of interest as an object of class 
  #AAStringSet.
  seq.prts = readAAStringSet(filepath = f.prts)
  ##  Molecular weight calculation.
  a.orf$MW = unlist(x = lapply(X = seq.prts,
                               FUN = function(y){mw(seq = y)}))/1000
  a.orf$MW = round(x = as.vector(x = a.orf$MW), digits = 2)
  ### Calculation the isoelectic point of a protein sequence.
  a.orf$pI = round(x = unlist(x = lapply(X = seq.prts, 
                                         FUN = function(y){pI(seq = y)})), 
                   digits = 2)
  ### Calculation the potential protein interaction index.
  a.orf$indexPPI = unlist(x = lapply(X = seq.prts, 
                                     FUN = function(y){boman(seq = y)}))
  a.orf$indexPPI = round(x = a.orf$indexPPI, digits = 2)
  ### Returning a final object.
  return(a.orf)
}
