#' Performs quality checks, then filters reads for quality
#' 
#' The function checks each fastq file specified in the "data file" for quality, and writes findings to a report. Following this, the function trims poor-quality bases and unknown bases from the ends of the sequences. Any reads which are too short, or contain any unknown bases (N), are removed from the file.
#' 
#' @param dataFile An R data frame with the data to be processed. The R object must contain the following headings: File, PE, Sample, Replicate, FilteredFile. MORE INFO ON THE DATA FILE
#' @param pairedEnd A logical. If false (default), a single-end protocol will be run. If true, a paired-end protocol will be run.
#' @param minLength An integer which specifies the minimum length for a read. Reads shorter than this length will be discarded. Default is 30 nucleotides.
#' @param Phred An integer which specifies Phred (ascii) quality score. Any two consecutive nucleotides with a quality score lower than this threshold will be discarded. Default score is 30.
#' @param blockSize An integer which specifies the number of reads to be read at a time when processing. Default is 1e8. 
#' @param readBlockSize An integer which specifies the number of bytes (characters) to be read at one time. Smaller \code{readBlockSize} reduces memory requirements, but is less efficient. Default is 1e5.
#' @return An object of class \code{QualityFilterResults}. Contains quality checks run before the filter, as well as summary statistics of the filtering. Also outputs directories with quality results, and filtered fastq.gz reads. 
#' @seealso \url{https://en.wikipedia.org/wiki/Phred_quality_score} for more about quality scores. 
#' @seealso \code{ShortRead} for more information about quality reports, \code{blockSize} (n) and \code{readerBlockSize}.
#' @details The function should be run in the working directory, where all fastq files are found.
#' @details \code{runQAandFilter} iterates over each file specified in the "datafile". It runs a quality assessment from the \code{ShortRead} package. The \code{ShortRead} package (\url{https://bioconductor.org/packages/release/bioc/html/ShortRead.html}) contains more information about this step.
#' @details At the next step, the function filters and trims the reads for quality. This is done by iterating over chunks of reads in the fastq files at a time. The size of the chunks are decided by the "blockSize" and "readerBlockSize" parameters. More information about how this is done is available in the \code{ShortRead} package.  
#' @details * it removes any trailing or leadining N's from each sequence,
#' @details * it removes any reads wich still contain N's,
#'  
#' @details * it trims the trailing end when it finds a minimum of 2 poor-quality bases in a window of 5. The threshold for poor quality is determined by the parameter "Phred", where the Phred score is logarithmically related to the probability of errors at each base,
#' @details * it removes any reads shorter than a minimum length (this is specified by the "minLength" parameter).
#' @details The function produces a new set of fastq files which have been filtered. The user must specify in the "FILTEREDFILE" column of the data file the output file. The user may specify the same output file for multiple input files - this will append new output to existing files, thereby allowing de-multiplexing of samples which have been run on different lanes.   
#' Quality reports are output to the working directory. A new R object (\code{QualityFilterResults}) is also created, which contains the pre-filtered quality reports, as well as a summary of how many reads have been trimmed or removed. 
#' @export
#' @import ShortRead
runQAandFilter <- function(dataFile, pairedEnd = FALSE, minlength = 30, Phred = 25, blockSize = 1e8, readerBlockSize = 1e5){
        
        # extract single-end names found in datafile to list
        dataFile$PE <- gsub(" ", "", dataFile$PE)
        dataFile <- dataFile[which(dataFile$PE == "SE"), ]
        
        print("QA results will be output to the 'QA' folder")
        
        if (pairedEnd == FALSE){ 
                run <- mapply(SEFilterAndTrim, dataFile$FILE, dataFile$FILTEREDFILE, minlength = minlength, Phred = Phred, blockSize = blockSize, readerBlockSize = readerBlockSize, SIMPLIFY = F)
                
                return(run)               
        }
} 

# private function
SEFilterAndTrim <- function(file, destination, minlength = minlength, Phred = Phred, blockSize = blockSize, readerBlockSize = readerBlockSize){
        
        file <- as.character(file)
        destination <- as.character(destination)
        fileExt <- gsub(".fastq.gz", "", file)
        
        # pre-filter quality check
        QASum_prefilter <- qa(file, type = "fastq")
        
        QADest <- paste("QA/", fileExt, sep = "")
        QASum_prefilterRpt <- report(x = QASum_prefilter, dest = QADest, type = "html")
        
        # open input stream
        stream1 <- FastqStreamer(file, n = blockSize, readerBlockSize = readerBlockSize)
        on.exit(close(stream1))
        
        destination <- as.character(destination)
        
        # define variables
        N_reads_in <- N_filt_reads <- Quality_filt_reads <- minLenFqa <- N_reads_out <- N_trim_N <- N_trim_Q <- 0L
        
        repeat {
                # input chunk
                fq1 <- yield(stream1)
                fq1Len<-length(fq1)
                if (fq1Len == 0)
                        break
                
                # count number and length of reads in
                N_reads_in<-N_reads_in+fq1Len
                init_length <-max(width(fq1))
                
                #trim leading/trailing N
                pos_1 <- trimEnds(sread(fq1), "N", relation="==", ranges=T)
                fqa<-narrow(fq1, start(pos_1), end(pos_1))
                N_trim_N<-N_trim_N+length(fqa[which(width(fqa)<init_length)])
                
                # filter reads containing N
                failed_N_1<-nFilter()(fqa)
                toTrash_1<-which(failed_N_1@.Data==F)
                N_filt_reads<-N_filt_reads+length(toTrash_1)
                if(N_filt_reads>0) {
                        fqa<-fqa[-toTrash_1]
                }
                
                # trim as soon as 2 of 5 nucleotides has quality encoding less than phred score tres
                #determine quality ascii character for trimming
                treschar<-names(encoding(quality(fq1))[which(encoding(quality(fq1))==Phred)])
                fqa <- trimTailw(fqa, 2, treschar, 2)
                N_trim_Q<-N_trim_Q+length(fqa[which(width(fqa)<init_length)])
                
                # drop reads that are less than x nt
                minLenFqa <- minLenFqa + length(fqa[width(fqa) < minlength])
                fqa <- fqa[width(fqa) >= minlength]
                
                N_reads_out <- N_reads_out+length(fqa)
                fullLengthReadsOut <- N_reads_out-N_trim_N-N_trim_Q
                
                writeFastq(object = fqa, file = destination, mode = "a",compress = T)
        }
        
        # collect results into one object
        attr(file, "inputFile") <- file
        attr(file, "prefilterQA") <- QASum_prefilter
        
        attr(file, "filter") <-
                data.frame(readsIn = N_reads_in, filterN = N_filt_reads, filterMinLen = minLenFqa,  readsOut = N_reads_out, trimForN = N_trim_N, trimForQual = N_trim_Q, fullLengthReadsOut = fullLengthReadsOut)
        attr(file, "outputFile") <- destination
        class(file) <- "QualityFilterResults"
        return(file)
}