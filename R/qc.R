#' Performs quality checks, then filters reads for quality
#' 
#' DESCRIPTION GOES HERE
#' 
#' @param dataFile A character string that is the name of the file to be processed. The file used must be in a format identical to that of "data.csv".
#' @param pairedEnd A logical. If false (default), a single-end protocol will be run. If true, a paired-end protocol will be run
#' @param minLength An integer which specifies the minimum length for a read. Reads shorter than this length will be discarded. Default is 30 nucleotides
#' @param Phred An integer which specifies Phred (ascii) quality score. Any two consecutive nucleotides with a quality score lower than this threshold will be discarded. Default score is 30
#' @return An object of class QualityFilterResults. Contains quality checks run before and after the filter, as well as summary statistics of the filtering. Also outputs directories with quality results, and filtered fastq.gz reads 
#' @seealso \url{https://en.wikipedia.org/wiki/Phred_quality_score} for more about quality scores. 
#' @details ADD MORE DETAILS
#' @export
#' @import ShortRead
runQAandFilter <- function(dataFile = "data.csv", pairedEnd = FALSE, minlength = 30, Phred = 25){
        
        # open data file
        dataFile <- read.csv(file = dataFile, sep = ",")
        
        # extract single-end names found in datafile to list
        dataFile$PE <- gsub(" ", "", dataFile$PE)
        SEdata <- dataFile[which(dataFile$PE == "SE"), ]
        SElist <- SEdata[,1]
        
        if (pairedEnd == FALSE){ 
                run <- lapply(SElist, SEFilterAndTrim, minlength = minlength, Phred = Phred)
        }
} 

# private function
SEFilterAndTrim <- function(file, minlength = minlength, Phred = Phred){
        
        fileExt <- paste(file, ".fastq.gz", sep = "")
        
        # pre-filter quality check
        QASum_prefilter <- qa(fileExt, type = "fastq")
        QADest <- paste(file, "/PreFiltQC", sep = "")
        QASum_prefilterRpt <- report(x = QASum_prefilter, dest = QADest, type = "html")
        
        
        # open input stream
        stream1 <- FastqStreamer(fileExt, readerBlockSize = 60000, n = 60000)
        on.exit(close(stream1))
        
        destination <- gsub(pattern = ".fastq.gz", replacement = "-filt.fastq.gz", x = fileExt)
        
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
                
                writeFastq(object = fqa, file = destination, mode = "a",compress = T)
        }
        
        # post-filter quality check
        QASum_postfilter <- qa(destination, type = "fastq")
        QADest <- paste(file, "/PostFiltQC", sep = "")
        QASum_postfilterRpt <- report(x = QASum_postfilter, dest = QADest, type = "html")
        
        
        # collect results into one object
        attr(destination, "name") <- paste(file, "-filt", sep="")
        attr(destination, "prefilterQA") <- QASum_prefilter
        attr(destination, "filter") <-
                data.frame(readsIn = N_reads_in, filterN = N_filt_reads, filterMinLen = minLenFqa,  readsOut = N_reads_out, trimForN = N_trim_N, trimForQual = N_trim_Q, fullLengthReadsOut = readsOut-trimForN-trimForQual)
        attr(destination, "postfilterQA") <- QASum_postfilter
        class(destination) <- "QualityFilterResults"
        return(destination)
}

