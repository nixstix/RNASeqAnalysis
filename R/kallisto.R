#' A number of transcriptome indices have already been built, for use with \code{kallistoQuant}
#' 
#' A number of transcriptome index files (.idx files) have been pre-built and provided with the package. These files may be used as input to the \code{kallsitoQuant} pseudoalignment function (or the user may provide his/her own transcriptome file). The \code{availableReferences} function provides a list of these ready-available transcriptomes.   
#' 
#' @param none No parameters.
#' @return list A list of available transcriptome indices (.idx files).
#' @examples 
#' availableReferences()
#'  
#' @export
availableReferences <- function(){
        extdata.dir <- system.file("extdata", "Kallisto-index", package="RNASeqAnalysis")
        extdata.dir <- list.files(extdata.dir)
        extdata.dir <- grep(".idx", extdata.dir, value = TRUE)
        return(as.list(extdata.dir))
        
}

#' Builds a transcriptome index (.idx) from a FASTA-formatted transcriptome.
#' 
#' \code{kallistoIndex} builds an index of the transcriptome. Once built, the transcriptome index (in the form of a De Bruijn graph) may be used for pseudoalignment of the target RNA-Seq data in order to quantify expression (see \code{kallistoQuant}). 
#' Kallisto must be installed locally. Kallisto can be downloaded from \url{https://pachterlab.github.io/kallisto/download.html}.
#' @param refTranscriptome A character string of the transcriptome file (in
#'   fasta format) to be indexed.
#' @param indexName A character string of the output index file to be created.
#'   Default is "./transcripts.idx").
#' @return An index file (.idx) which can be used to quantify expression data using
#'   \code{\link{kallistoQuant}} pseudoalignment.
#' @examples
#' kallistoIndex("./transcripts.fasta.gz")
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
#' Kallisto quantifies transcript abundance from RNA-Seq data. 
#' 
#' The data is broken down into k-mers and each k-mer is pseudoaligned to k-mers in the index.
#' Because kallisto doesn't rely on full alignment, it is much quicker than other methods, without losing accuracy.
#' @param file1 A character string of the name of the RNA-Seq data file (fastq.gz) to be processed.
#' @param file2 A character string of the RNA-Seq data file (fastq.gz) to be processed - in the case there is paired-end data. 
#' @param refIndex A character string of the name of the index file (usually, .idx) against which the fastq files are pseudoaligned. This file can be generated from a fasta-formatted transcriptome file using \code{\link{kallistoIndex}}.
#' @param pairedEnd is a logical. If \code{true} (default), a paired end protocol is chosen (for this, the file1 and file2 parameters must be specified). If \code{false}, a single end protocol will be run, and only file1 will be processed.
#' @param bootstrap An integer specifying the number of times to bootstrap. The output will provide a measure of variance.
#' @param fragmentLength An integer. Estimated fragment length. Only required for single-end data (for paired-end data, Kallisto is able to calculate the fragment length). Default is 200 bp.
#' @param fragmentSD A numeric. The standard deviation of the fragment length. Only required for single-end data (for paired-end data, Kallisto is able to calculate the standard deviation). Default is 10 bp. 
#' @param refIndexFromPackage A logical, false by default. If the user would like to use an index provided by the package (these can be viewed using \code{availableReferences()}, he/she should specifiy "true" here. Otherwise, the function will assume the index is provided by the user.
#' #' @param mc.cores The number of cores to use when parallelizing. Default is 1 (i.e. no parallelisation)
#' @return A data frame of the estimated abundances of each transcript specified in the input file(s). The data frame is also saved to a folder, which is given the title of the files, exlcuding the extensions. The folder contains these abundances (in text format and compressed format), as well as information about the run.
#' @export
kallistoQuant <- function(dataFile, preFilt = FALSE, refIndex, refIndexFromPackage = FALSE, bootstrap = 0, fragmentLength = 200, fragmentSD = 10, biasCor = FALSE, mc.cores = 1){
        
        dataFileSE <- dataFileSE(dataFile)
        dataFilePE <- dataFilePE(dataFile)
        
        try({
                
                # if user selects index supplied by package   
                if (refIndexFromPackage == TRUE){
                        refIndex <- system.file("extdata", refIndex, package="RNASeqAnalysis")
                }
                
                # if bootstrap is 0, make sure threads is 1   
                if (bootstrap == 0){
                        threads <- 1
                        # otherwise, threads == bootstrap        
                } else {
                        threads <- bootstrap 
                }
                
                # if user selects bias correction        
                if (biasCor == TRUE){
                        cmd <- paste("kallisto quant --plaintext --bias -i", refIndex, "-b", bootstrap, "-t", threads, "-o ", sep = " ")        
                } else {
                        cmd <- paste("kallisto quant --plaintext -i", refIndex, "-b", bootstrap, "-t", threads, "-o ", sep = " ")
                }
                
                if(nrow(dataFileSE) > 0) {
                        runSE <- mcmapply(kallistoQuantRunSE, cmd, preFilt = preFilt, dataFileSE$FILE, dataFileSE$FILTEREDFILE, fragmentLength, fragmentSD, mc.cores = mc.cores, mc.preschedule = FALSE) 
                } else {
                        runSE <- NULL
                }
 
                if(nrow(dataFilePE) > 0) {
                        runPE <- mcmapply(kallistoQuantRunPE, cmd, preFilt = preFilt, prefilt = dataFilePE$FILE.x, postfilt = dataFilePE$FILTEREDFILE.x, prefilt2 = dataFilePE$FILE.y, postfilt2 = dataFilePE$FILTEREDFILE.y, mc.cores = mc.cores, mc.preschedule = FALSE) 
                } else {
                        runPE <- NULL
                }        
        })
}

# private function
kallistoQuantRunSE <- function(cmd = cmd, preFilt = preFilt, prefilt = dataFileSE$FILE, postfilt = dataFileSE$FILTEREDFILE, fragmentLength = fragmentLength, fragmentSD = fragmentSD){
        
        if (preFilt == FALSE){
                file <- postfilt
        } else {
                file <- prefilt
        }
        
        outfile <- gsub(".fastq.gz", "", file)
        cmd <- paste(cmd, outfile, "--single -l", fragmentLength, "-s", fragmentSD, file, sep = " ")
        
        print(cmd)
        system(cmd)
}

# private function
kallistoQuantRunPE <- function(cmd = cmd, preFilt = preFilt, prefilt = dataFilePE$FILE.x, postfilt = dataFilePE$FILTEREDFILE.x, prefilt2 = dataFilePE$FILE.y, postfilt2 = dataFilePE$FILTEREDFILE.y){
        
        if (preFilt == FALSE){
                file <- postfilt
                file2 <- postfilt2
        } else {
                file <- prefilt
                file2 <- prefilt2
        }
        
        outfile <- gsub("_1.fastq.gz", "", file)
        
        cmd <- paste(cmd, outfile, file, file2, sep = " ")
        print(cmd)
        system(cmd)

}


