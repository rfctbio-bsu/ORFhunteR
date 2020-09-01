#' Identify all possible variants of open reading frames in a nucleotide sequence
#' @description Identify all possible variants of open reading frames in a nucleotide sequence of interest.
#' @param x character string giving the nucleotide sequence of interest.
#' @param cod.start character string with type of start codon: "ATG", "GTG", "TTG" or "CTG".
#'                  Default value is "ATG".
#' @param f.orf character string giving the name of output TXT file in tab-delimited format with
#'              coordinates of identified open reading frames. Default value is NULL.
#' @param d.work character string giving the path to and name of work directory. Concatenated with filename via "/". Default value is NULL.
#' @return matrix with start and stop positions, length and sequence of identified variants of open reading frames.
#' @author Vasily V. Grinev
#' @examples
#' x = "AAAATGGCTGCGTAATGCAAAATGGCTGCGAATGCAAAATGGCTGCGAATGCCGGCACGTTGCTACGT"
#' orf = finderORFs(x = x,
#'                 cod.start = "ATG",
#'                 f.orf = NULL,
#'                 d.work = NULL)
#' @export

finderORFs <- function(x, cod.start = "ATG", f.orf = NULL, d.work = NULL){
  ### Calculation coordinates of defined start codon and all stop codons 
  #in the sequence of interest.
  cod.pos = codonStartStop(x)
  START = cod.pos[[1]][cod.pos[[2]] == cod.start]
  STOP = cod.pos[[1]][cod.pos[[2]] == "TAA" | cod.pos[[2]] == "TAG" | 
                        cod.pos[[2]] == "TGA"]
  ### Development a matrix of calculated coordinates.
  in.frame = outer(X = STOP, Y = START, FUN = "-")
  L = in.frame[lower.tri(x = in.frame, diag = TRUE)]/3
  L = round(x = L) == L & L > 0
  if (length(L[L == "TRUE"]) > 0){
    ### Calculation of the in-frame start codons.
    starts = t(x = in.frame)
    starts[] = START
    starts = t(x = starts)
    starts = starts[lower.tri(x = starts, diag = TRUE)]
    starts = starts[L]
    ### Calculation of the in-frame stop codons.
    stops = in.frame
    stops[] = STOP
    stops = stops[lower.tri(x = stops, diag = TRUE)]
    stops = stops[L]
    ### Calculation the coordinates of open reading frames.
    ORFs = cbind(starts, stops)
    ORFs = aggregate(x = ORFs[, 2], by = list(ORFs[, 1]), FUN = min)
    ORFs = aggregate(x = ORFs[, 1], by = list(ORFs[, 2]), FUN = min)
    ORFs = ORFs[order(ORFs[, 2]), ]
    seq.ORFs = substring(x, first = ORFs[, 2], last = ORFs[, 1] + 2)
    ### Final object.
    ORFs = list(ORFs[, 2], ORFs[, 1] + 2, ORFs[, 1] - ORFs[, 2] + 3, seq.ORFs)
    if (!is.null(f.orf)){
      ORFs = do.call(what = cbind, args = ORFs)
      colnames(ORFs) = c("start", "end", "length", "orf.sequence")
      write.table(x = ORFs,
                  file = paste(d.work, f.orf, sep = "/"),
                  sep = "\t",
                  quote = FALSE,
                  col.names = TRUE,
                  row.names = FALSE)
    }
    return(ORFs)
  }
}
