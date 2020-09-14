#' Extract the sequences of identified open reading frames
#' @description Extract the sequences of identified open reading frames.
#' @param f.orf character string giving the name of tab-delimited TXT file with coordinates of open
#'             reading frame(-s). This file should include three mandatory fields:
#'             i) transcript_id - transcript ID;
#'             ii) start        - start coordinate of open reading frame in a transcript;
#'             iii) end         - end coordinate of open reading frame in a transcript.
#' @param f.trans character string giving the name of file with transcript(-s) of interest.
#' @param frt format of file with transcript(-s) of interest. Valid formats are "gtf", "gff",
#'             "fasta" and "fa". Default value is "fasta".
#' @param f.seq character string giving the name of output FASTA file with sequences of identified
#'             open reading frames. Default value is NULL.
#' @param d.work character string giving the path to and name of work directory. Concatenated with filename via "/". Default value is NULL.
#' @return DNAStringSet object with sequences of extracted open reading frames.
#' @author Vasily V. Grinev
#' @examples
#' f.orf_path <- system.file("extdata", "Set.trans_ORFs.coordinates.txt", package = "ORFhunteR")
#' f.trans_path <- system.file("extdata", "Set.trans_sequences.fasta", package = "ORFhunteR")
#' seqORF = sequenceORFs(f.orf = f.orf_path,
#'                       f.trans = f.trans_path,
#'                       frt = "fasta",
#'                       f.seq = NULL,
#'                       d.work = NULL)
#' @export

sequenceORFs <- function(f.orf, f.trans, frt = "fasta", f.seq = NULL, 
                         d.work = NULL){
  ### Validation of file format.
  if (!frt %in% c("gtf", "gff", "fasta", "fa")){
    stop("Invalid file format")
  }
  ### Loading of required auxiliary libraries.
  #   This code was successfully tested with libraries Biostrings v.2.52.0, 
  #rtracklayer v.1.44.3 and
  #   BSgenome.Hsapiens.UCSC.hg38 v.1.4.1.
  #suppressMessages(expr = library(package = Biostrings))
  #suppressMessages(expr = library(package = rtracklayer))
  #suppressMessages(expr = library(package = BSgenome.Hsapiens.UCSC.hg38))
  ##Loading of open reading frames coordinates as an object of class data frame.
  if(!is.null(d.work)){
    f.orf <- paste(d.work, f.orf, sep = "/")
    f.trans <- paste(d.work, f.trans, sep = "/")
    f.seq <- paste(d.work, f.seq, sep = "/")
  }
  seq.orf <- read.table(file = f.orf,
                       header = TRUE,
                       quote = "\"",
                       as.is = TRUE)
  seq.orf = seq.orf[order(seq.orf$transcript_id), ]
  ### Loading of nucleotide sequences of interest as an object of class 
  #DNAStringSet.
  if (frt %in% c("fasta", "fa")){
    seq.set = readDNAStringSet(filepath = f.trans)
    seq.set = seq.set[order(names(seq.set))]
  }else{
    ##  Loading of the GTF/GFF annotations as a GRanges object.
    gtf = sort(x = import(con = paste(d.work, f.gtf, sep = "/"), format = "gtf"))
    gtf = gtf[gtf$type == "exon", ]
    ##  Converting a GRanges object into transcript-splitted GRangesList object.
    list.gtf = split(x = gtf, f = gtf$transcript_id)
    ##  Retrieving of transcript sequences as a list of DNAString instances. 
    #This step has not been optimized and may take a long time. Be patient!
    list.seq = lapply(X = list.gtf, FUN = function(x){
      if (as.character(strand(x)@values) == "-"){
      unlist(x = rev(x = getSeq(Hsapiens, x)))
    }else{unlist(x = getSeq(Hsapiens, x))}})
    ##  Consolidation of the DNAString objects into a DNAStringSet instance.
    seq.set = DNAStringSet(x = list.seq, use.names = TRUE)
    seq.set = seq.set[order(names(seq.set))]
  }
  ### Extraction of sequences of identified open reading frames.
  seq.orf <- DNAStringSet(x = substr(x = seq.set, start = seq.orf$start, 
                                     stop = seq.orf$end), use.names = TRUE)
  if (!is.null(f.seq)){
    ### Write an XStringSet object to a file.
    writeXStringSet(x = seq.orf, filepath = f.seq)
  }
  ### Return the final object of class DNAStringSet.
  return(seq.orf)
}
