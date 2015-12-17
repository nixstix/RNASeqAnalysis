#' Builds a transcriptome index from a FASTA formatted file of target sequences
#' 
#' DESCRIPTION GOES HERES. Kallisto must be installed locally. Kallisto can be downloaded from \url{https://pachterlab.github.io/kallisto/download.html}.
#' @param refTranscriptome A character string of the transcriptome file (in
#'   fasta format) to be indexed.
#' @param indexName A character string of the output index file to be created.
#'   Default is "./transcripts.idx").
#' @return An index file which can be used to quantify expression data using
#'   \code{\link{kallistoQuant}} pseudoalignment.
#' @examples
#' kallistoIndex("./transcripts.fasta")
#' kallistoIndex("./transcripts.fa.gz", "./index.idx")
#'  
#' @seealso \url{http://pachterlab.github.io/kallisto/} for more about Kallisto.
#' @export
kallistoIndex <- function(refTranscriptome, indexName="./transcripts.idx"){
    cmd <- paste("kallisto index -i", indexName, refTranscriptome, sep = " ")
    print(cmd)
    system(cmd)
}

#' @export
kallistoIndex2 <- function(refTranscriptome, packageTranscriptome = FALSE, indexName="./transcripts.idx"){
        
        if (packageTranscriptome == TRUE){
                print("belongs to package")
                # how do we access this reference?
        } else {
                cmd <- paste("kallisto index -i", indexName, refTranscriptome, sep = " ")
                print(cmd)
                system(cmd)        
        }
}

# NOTE TO DO:
# create index using transcriptomes supplied by the package

#Quantifies transcript expression
#'
#'DESCRIPTION GOES HERE
#'
#' @param file1 A character string of the name of the (fastq.gz) file to be processed.
#' @param file2 A character string of the name of the second (fastq.gz) file to be processed. This is used in the case that paired-end data is available. 
#' @param refIndex A character string of the name of the index file (usually, .idx) against which the fastq files are pseudoaligned. This file can be generated from a transcriptome (fasta) file using \code{\link{kallistoIndex}}.
#' @param pairedEnd is a logical. If \code{true} (default), a paired end protocol is chosen (for this, the file1 and file2 parameters must be specified). If \code{false}, a single end protocol will be run, and only file1 will be processed.
#' @param bootstrap An integer specifying the number of times MORE.
#' @param fragmentLength An integer. MORE.
#' @param fragmentSD A numeric. MORE.
#' @return A data frame of the estimated abundances of each transcript specified in the input file(s). The data frame is also saved to a folder, which is given the title of the files, exlcuding the extensions. The folder contains these abundances (in text format and compressed format), as well as information about the run.
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

