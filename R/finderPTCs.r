#' Identify the premature termination codons in nucleotide sequences
#' @description Identify the premature termination codons in nucleotide sequences of interest.
#' @param f.orf character string giving the name of TXT file in tab-delimited format with
#'             coordinates of open reading frames.
#'             This file should include four mandatory fields:
#'             i) transcript_id - transcript ID;
#'             ii) start        - start coordinate of open reading frame in a transcript;
#'             iii) end         - end coordinate of open reading frame in a transcript;
#'             iv) length       - length of open reading frame.
#' @param f.trans character string giving the name of file with transcript(-s) of interest.
#' @param frt format of file with transcript(-s) of interest. Default value is "gtf".
#' @param d.work character string giving the path to and name of work directory. Concatenated with filename via "/". Default value is NULL.
#' @return data.frame object with start and stop positions, length and stop status of codons for each transcript ID.
#' @author Vasily V. Grinev
#' @examples
#' f.orf_path <- system.file("extdata", "Set.trans_ORFs.coordinates.txt", package = "ORFhunteR")
#' f.trans_path <- system.file("extdata", "Set.trans_sequences.gtf", package = "ORFhunteR")
#' ptc <- finderPTC(f.orf = f.orf_path,
#'                  f.trans = f.trans_path,
#'                  frt = "gtf",
#'                  d.work = NULL)
#' @export

finderPTC <- function(f.orf, f.trans, frt = "gtf", d.work = NULL){
  ### Validation of file format.
  if (!frt %in% c("gtf", "gff")){
    stop("Invalid file format")
  }
  ##Loading of open reading frames coordinates as an object of class data frame.
  if(!is.null(d.work)){
    f.orf <- paste(d.work, f.orf, sep = "/")
    f.trans <- paste(d.work, f.trans, sep = "/")
  }
  orf <- read.table(file = f.orf,
                   #sep = "\t",
                   header = TRUE,
                   quote = "\"",
                   as.is = TRUE)
  orf$tr_length = 0
  orf$l.ex_width = 0
  orf = as.matrix(orf)
  ### Loading of experimental transcripts in GTF/GFF format.
  ##  Loading of required auxiliary libraries.
  #suppressMessages(expr = library(package = rtracklayer))
  ##  Loading of experimental transcripts as a list of character strings.
  seq.set = rtracklayer::import(con = f.trans, format = frt)
  seq.set = sort(x = seq.set[seq.set$type == "exon", ])
  seq.set = split(x = seq.set, f = seq.set$transcript_id)
  ### Identification of PTC(-s).
  l.ex_width = lapply(X = seq.set,
        FUN = function(y){if (as.character(rtracklayer::strand(y))[1] == "-"){
                        return(c(sum(width(y)), width(y[1])))
                      }else{
                        return(c(sum(width(y)), width(y[length(x = y)])))
                      }})
  l.ex_width = do.call(what = rbind, args = l.ex_width)
  idx = orf[, 1] %in% rownames(l.ex_width)
  orf[idx, 5:6] = l.ex_width[orf[idx, 1], 1:2]
  orf = data.frame(orf, stringsAsFactors = FALSE)
  orf[, 2:6] = apply(X = orf[, 2:6], MARGIN = 2, FUN = as.numeric)
  orf$stop.status = "mature"

  premature_controll <- orf[, 3] - (orf[, 5] - orf[, 6] + 1) <= -50
  if(any(premature_controll)){
    orf[premature_controll, ]$stop.status = "premature"
  } 
  orf = orf[, c(1:4, 7)]
  ### Returning a final object.
  return(orf)
}
