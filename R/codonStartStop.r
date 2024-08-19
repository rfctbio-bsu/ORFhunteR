#' Identify all potential start and stop codons in a nucleotide sequence
#' @description This function scans a nucleotide sequence of interest
#'     in the search of start codon ATG or non-canonical start codons "CTG", 
#'     "GTG", "AGG", "TTG", "AAG", "ACG", "ATC", "ATT" or "ATA" (according to 
#'     Chothani S. P. et al. Mol Cell. 2022 Aug 4;82(15):2885-2899.e8) as well
#'     as stop codons TAA, TAG and TGA.
#' @param x character string giving the nucleotide sequence.
#' @return list of potential start and stop codons with their coordinates.
#' @author Vasily V. Grinev
#' Last updated: June 12, 2024
#' @examples
#' codons <- codonStartStop(x="AAAATGGCATGGTAAGTCAAAATGGCATGGTAAGTCAAAATGGCGG")
#' @export

codonStartStop <- function(x){
  ### Calculation of codon positions.
  codons <- DNAStringSet(x=c("ATG", "CTG", "GTG", "AGG", "TTG", "AAG", "ACG",
                             "ATC", "ATT", "ATA", "TAA", "TAG", "TGA"))
  names(x=codons) <- as.character(x=codons)
  codonPositions <- sort(x=unlist(x=matchPDict(pdict=codons,
                                               subject=DNAString(x=x))))
  codonPositions <- list(start(x=codonPositions), names(x=codonPositions))
  ### Returning a final object of class list.
  return(codonPositions)
}
