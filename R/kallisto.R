#' Builds a transcriptome index from a FASTA formatted file of target sequences
#' 
#' DESCRIPTION GOES HERES. Kallisto must be installed locally in order for this function to work. Kallisto can be downloaded from \url{https://pachterlab.github.io/kallisto/download.html}
#' @param refTranscriptome A character string of the transcriptome file (in
#'   fasta format) to be indexed
#' @param indexName A character string of the output index file to be created.
#'   Default is "./transcripts.idx")
#' @return An index file which can be used to quantify expression data using
#'   \code{\link{kallistoQuant}} pseudoalignment
#' @examples
#' kallistoIndex("./transcripts.fasta")
#' kallistoIndex("./transcripts.fa.gz", "./index.idx")
#'  
#' @seealso \url{http://pachterlab.github.io/kallisto/} for more about Kallisto 
#' @export
kallistoIndex <- function(refTranscriptome, indexName="./transcripts.idx"){
    
    cmd <- paste("kallisto index -i", indexName, refTranscriptome, sep = " ")
    system(cmd)
}

# NOTE TO DO:
# create index using transcriptomes supplied by the package
# list all transcriptomes and indexes supplied by the package


#extdata.dir <- system.file("extdata", package="RNASeqAnalysis")
#fastaFiles <- dir(extdata.dir, pattern=refTranscriptome, full.names=TRUE)
#read.table(fastaFiles)


#' @export
kallistoQuant <- function(){
    print("TO DO")
}

