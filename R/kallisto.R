#' Builds a transcriptome index from a FASTA formatted file of target sequences
#' 
#' DESCRIPTION GOES HERES. Kallisto must be installed locally. Kallisto can be downloaded from \url{https://pachterlab.github.io/kallisto/download.html}
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
    print(cmd)
    system(cmd)
}

# NOTE TO DO:
# create index using transcriptomes supplied by the package


#' @export
kallistoQuant <- function(file1, file2, refIndex, pairedEnd = TRUE, bootstrap = 0, fragmentLength = 200, fragmentSD = 10){
    cmd <- paste("kallisto quant -i", refIndex, "-b", bootstrap, "-o ", sep = " ")
    if (pairedEnd == FALSE){
            outfile <- outfile <- gsub(".fastq.gz", "", file1)
            cmd <- paste(cmd, outfile, "--single -l", fragmentLength, "-s", fragmentSD, file1, sep = " ")
    } else{
            outfile <- gsub("_1.fastq.gz", "", file1)
            cmd <- paste(cmd, outfile, file1, file2, sep = " ")
    } 
    
    
    print(cmd)
    system(cmd)
    abundances <- read.table(file = paste(outfile, "abundance.tsv", sep="/"), header = TRUE, sep = "\t")
    return(abundances)
}

# NOTE TO DO:
# do something with bias correction

#' A number of transcriptome references have already been built, for use with \code{kallistoQuant}
#' 
#' DESCRIPTION GOES HERE
#' 
#' @param None No parameters
#' @return A list of available transcriptomes (.idx files)
#' @examples 
#' availableReferences()
#'  
#' @export
availableReferences <- function(){
        extdata.dir <- system.file("extdata", package="RNASeqAnalysis")
        extdata.dir <- list.files(extdata.dir)
        extdata.dir <- grep("idx", extdata.dir, value = TRUE)
}

