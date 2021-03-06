% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findORFs.r
\name{findORFs}
\alias{findORFs}
\title{Identify all variants of open reading frames in a nucleotide sequence}
\usage{
findORFs(x, codStart = "ATG")
}
\arguments{
\item{x}{character string giving the nucleotide sequence of interest.}

\item{codStart}{character string with type of start codon: "ATG", "GTG",
"TTG" or "CTG". Default value is "ATG".}
}
\value{
matrix with start and stop positions, length and sequence of
    identified variants of open reading frames.
}
\description{
Identify all possible variants of open reading frames
    in a nucleotide sequence of interest.
}
\examples{
x <- "AAAATGGCTGCGTAATGCAAAATGGCTGCGAATGCAAAATGGCTGCGAATGCCGGCACGTTGCTACGT"
orf <- findORFs(x = x, codStart = "ATG")
}
\author{
Vasily V. Grinev
}
