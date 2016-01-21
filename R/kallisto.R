#' A number of transcriptome indices have already been built, for use with \code{kallistoQuant}
#' 
#' A number of transcriptome index files (.idx files) have been pre-built and provided with the package. These files may be used as input to the \code{kallsitoQuant} pseudoalignment function (or the user may provide his/her own transcriptome file). The \code{availableReferences} function provides a list of these ready-available transcriptomes.   
#' 
#' @param None No parameters
#' @return A list of available transcriptome indices (.idx files)
#' @examples 
#' availableReferences()
#'  
#' @export
availableReferences <- function(){
        extdata.dir <- system.file("extdata", package="RNASeqAnalysis")
        extdata.dir <- list.files(extdata.dir)
        extdata.dir <- grep(".idx", extdata.dir, value = TRUE)
        return(as.list(extdata.dir))
        
}

#' Builds a transcriptome index (.idx) from a FASTA-formatted transcriptome.
#' 
#' \code{kallistoIndex} builds an index of  Once built, the transcriptome index may be used for pseudoalignment of the target RNA-Seq data in order to quantify expression (see \code{kallistoQuant}). Kallisto must be installed locally. Kallisto can be downloaded from \url{https://pachterlab.github.io/kallisto/download.html}.
#' @param refTranscriptome A character string of the transcriptome file (in
#'   fasta format) to be indexed.
#' @param indexName A character string of the output index file to be created.
#'   Default is "./transcripts.idx").
#' @return An index file (.idx) which can be used to quantify expression data using
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

#Quantifies transcript expression
#'
#'DESCRIPTION GOES HERE
#'
#' @param file1 A character string of the name of the (fastq.gz) file to be processed.
#' @param file2 A character string of the name of the second (fastq.gz) file to be processed. This is used in the case that paired-end data is available. 
#' @param refIndex A character string of the name of the index file (usually, .idx) against which the fastq files are pseudoaligned. This file can be generated from a fasta-formatted transcriptome file using \code{\link{kallistoIndex}}.
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
#' @param refIndexFromPackage A logical, false by default. MORE.
#' @return A data frame of the estimated abundances of each transcript specified in the input file(s). The data frame is also saved to a folder, which is given the title of the files, exlcuding the extensions. The folder contains these abundances (in text format and compressed format), as well as information about the run.
#' @export
kallistoQuant2 <- function(file1, file2, refIndex, refIndexFromPackage = FALSE, pairedEnd = TRUE, bootstrap = 0, fragmentLength = 200, fragmentSD = 10){
        
        if (refIndexFromPackage == TRUE){
                        refIndex <- system.file("extdata", "kallisto.idx", package="RNASeqAnalysis")
                }
        
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

