% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kallisto.R
\name{kallistoIndex}
\alias{kallistoIndex}
\title{Builds a transcriptome index (.idx) from a FASTA-formatted transcriptome.}
\usage{
kallistoIndex(refTranscriptome, indexName = "./transcripts.idx")
}
\arguments{
\item{refTranscriptome}{A character string of the transcriptome file (in
fasta format) to be indexed.}

\item{indexName}{A character string of the output index file to be created.
Default is "./transcripts.idx").}
}
\value{
An index file (.idx) which can be used to quantify expression data using
  \code{\link{kallistoQuant}} pseudoalignment.
}
\description{
\code{kallistoIndex} builds an index of the transcriptome. Once built, the transcriptome index (in the form of a De Bruijn graph) may be used for pseudoalignment of the target RNA-Seq data in order to quantify expression (see \code{kallistoQuant}). 
A transcriptome De Bruijn graph (T-DBG) is constructed from k-mers in the transcriptome, where each k-mer in the transcriptome is a node in the T-DBG. Each transcript is represented by a different colour, and whenever a node appears in a transcript, it is given the colour of that transcript. MORE. 
Kallisto must be installed locally. Kallisto can be downloaded from \url{https://pachterlab.github.io/kallisto/download.html}.
}
\examples{
kallistoIndex("./transcripts.fasta")
kallistoIndex("./transcripts.fa.gz", "./index.idx")
 
}
\seealso{
\url{http://pachterlab.github.io/kallisto/} for more about Kallisto.
}

