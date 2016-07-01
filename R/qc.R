#' Performs quality checks, then filters reads for quality
#' 
#' The function trims poor-quality bases and unknown bases from the ends of the sequences. Any reads which are too short, or contain any unknown bases (N), are removed from the file.
#' 
#' @param dataFile An R data frame with the data to be processed. The R object is a standard format, and must contain the following headings: File, PE, Sample, Replicate, FilteredFile. More information about the file is available at \code{datafileTemplate}. 
#' @param minLength An integer which specifies the minimum length for a read. Reads shorter than this length will be discarded. Default is 30 nucleotides.
#' @param Phred An integer which specifies Phred (ascii) quality score. Any two consecutive nucleotides with a quality score lower than this threshold will be discarded. Default score is 30.
#' @param blockSize An integer which specifies the number of reads to be read at a time when processing. Default is 1e8. 
#' @param readBlockSize An integer which specifies the number of bytes (characters) to be read at one time. Smaller \code{readBlockSize} reduces memory requirements, but is less efficient. Default is 1e5.
#' @param mc.cores The number of cores to use when paralleling. Default is 1 (i.e. no parallelisation)
#' @return An object of class \code{QualityFilterResults}.  
#' @seealso \url{https://en.wikipedia.org/wiki/Phred_quality_score} for more about quality scores. 
#' @seealso \code{ShortRead} for more information about quality reports, \code{blockSize} (n) and \code{readerBlockSize}.
#' @details The function should be run in the working directory, where all fastq files are found.
#' @details \code{runQAandFilter} iterates over each file specified in the "datafile", and filters and trims the reads for quality. This is done by iterating over chunks of reads in the fastq files at a time. The size of the chunks are decided by the "blockSize" and "readerBlockSize" parameters. More information about how this is done is available in the \code{ShortRead} package.  
#' @details * it removes any trailing or leadining N's from each sequence,
#' @details * it removes any reads wich still contain N's,
#' @details * it trims the trailing end when it finds a minimum of 2 poor-quality bases in a window of 5. The threshold for poor quality is determined by the parameter "Phred", where the Phred score is logarithmically related to the probability of errors at each base,
#' @details * it removes any reads shorter than a minimum length (this is specified by the "minLength" parameter).
#' @details The function produces a new set of fastq files which have been filtered. The user must specify in the "FILTEREDFILE" column of the data file the output file. The user may specify the same output file for multiple input files - this will append new output to existing files, thereby allowing de-multiplexing of samples which have been run on different lanes.   
#' A new R object (\code{QualityFilterResults}) is created, which contains pointers to the input and output fastq files, as well as a summary of how many reads have been trimmed or removed. 
#' @export
#' @import ShortRead
filterBadSeqs <- function(dataFile, minlength = 30, Phred = 25, blockSize = 1e8, readerBlockSize = 1e5, mc.cores = 1){
        
        dataFileSE <- dataFileSE(dataFile)
        dataFilePE <- dataFilePE(dataFile)
        
        if(nrow(dataFileSE) > 0) {
                runSE <- mcmapply(filterAndTrimSE, dataFileSE$FILE, dataFileSE$FILTEREDFILE, dataFileSE$PE, minlength = minlength, Phred = Phred, blockSize = blockSize, readerBlockSize = readerBlockSize, SIMPLIFY = F, mc.cores = mc.cores, mc.preschedule = F) 
        } else {
                runSE <- NULL
        }
        
        if(nrow(dataFilePE) > 0) {
                runPE <- mcmapply(filterAndTrimPE, dataFilePE$FILE.x, dataFilePE$FILE.y, dataFilePE$FILTEREDFILE.x, dataFilePE$FILTEREDFILE.y, dataFilePE$PE.x, minlength = minlength, Phred = Phred, blockSize = blockSize, readerBlockSize = readerBlockSize, SIMPLIFY = F, mc.cores = mc.cores, mc.preschedule = F) }
        else{
                runPE <- NULL
        }
        
        
        run <- c(runSE, runPE)
        return(run)               
} 

# private function
filterAndTrimSE <- function(file, destination, PE, minlength = minlength, Phred = Phred, blockSize = blockSize, readerBlockSize = readerBlockSize, mc.cores = mc.cores, mc.preschedule = F){
        cat("Filtering: ", file, "\n", sep = " ")
        # open input stream
        stream1 <- FastqStreamer(file, n = blockSize, readerBlockSize = readerBlockSize)
        on.exit(close(stream1))
        
        # define variables
        N_reads_in <- N_filt_reads <- Q_filt_reads <- minLenFqa <- N_reads_out <- N_trim <- 0L
        
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
                
                # drop reads containing N
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
                Q_filt_reads <- Q_filt_reads + N_reads_in - N_filt_reads - length(fqa)
                
                # drop reads that are less than x nt
                minLenFqa <- minLenFqa + length(fqa[width(fqa) < minlength])
                fqa <- fqa[width(fqa) >= minlength]
                
                N_reads_out <- N_reads_out+length(fqa)
                
                # calculate how many reads have been trimmed
                N_trim <- N_trim + length(fqa[which(width(fqa) < init_length)])
                fullLengthReadsOut <- N_reads_out - N_trim
                writeFastq(object = fqa, file = destination, mode = "a",compress = T)
        }
        
        # collect results into one object
        df <- data.frame(filtFile = destination, readsIn = N_reads_in, filterN = N_filt_reads, filterQ = Q_filt_reads, filterMinLen = minLenFqa,  readsOut = N_reads_out, trim = N_trim, fullLengthReadsOut = fullLengthReadsOut, unmatchedPair = NA)
        rownames(df) <- file
        return(df)
        
}

# private function
filterAndTrimPE <- function(file, file2, destination, destination2, PE, minlength = minlength, Phred = Phred, blockSize = blockSize, readerBlockSize = readerBlockSize, mc.cores = mc.cores, mc.preschedule = F){
        
        cat("Filtering: ", file, file2, "\n", sep = " ")
        
        # halve read size to allow 2 chunks to open at the same time
        blockSize <- blockSize/2
        readerBlockSize <- blockSize/2
        
        # open input stream
        stream1 <- FastqStreamer(file, n = blockSize, readerBlockSize = readerBlockSize)
        on.exit(close(stream1))
        stream2 <- FastqStreamer(file2, n = blockSize, readerBlockSize = readerBlockSize)
        on.exit(close(stream2))
        
        # define variables
        N_reads_in <- N_filt_reads <- Q_filt_reads <- minLenFqa <- N_reads_out <- N_trim <- prRemoval <- 0L
        N_reads_in2 <- N_filt_reads2 <- Q_filt_reads2 <- minLenFqb <- N_reads_out2 <- N_trim2 <- prRemoval2 <- 0L
        
        repeat {
                # input chunk
                fq1 <- yield(stream1)
                fq1Len<-length(fq1)
                if (fq1Len == 0)
                        break
                fq2 <- yield(stream2)
                fq2Len<-length(fq2)
                if (fq2Len == 0)
                        break
                
                # count number and length of reads in
                N_reads_in<-N_reads_in+fq1Len
                init_length <-max(width(fq1))
                N_reads_in2<-N_reads_in2+fq2Len
                init_length2 <-max(width(fq2))
                
                #trim leading/trailing N
                pos_1 <- trimEnds(sread(fq1), "N", relation="==", ranges=T)
                fqa<-narrow(fq1, start(pos_1), end(pos_1))
                pos_2<- trimEnds(sread(fq2), "N", relation="==", ranges=T)
                fqb<-narrow(fq2, start(pos_2), end(pos_2))
                
                # drop reads containing N
                failed_N_1<-nFilter()(fqa)
                toTrash_1<-which(failed_N_1@.Data==F)
                N_filt_reads<-N_filt_reads+length(toTrash_1)
                if(N_filt_reads>0) {
                        fqa<-fqa[-toTrash_1]
                }
                failed_N_2<-nFilter()(fqb)
                toTrash_2<-which(failed_N_2@.Data==F)
                N_filt_reads2<-N_filt_reads2+length(toTrash_2)
                if(N_filt_reads2>0) {
                        fqb<-fqb[-toTrash_2]
                }
                
                # trim as soon as 2 of 5 nucleotides has quality encoding less than phred score tres
                #determine quality ascii character for trimming
                treschar<-names(encoding(quality(fq1))[which(encoding(quality(fq1))==Phred)])
                fqa <- trimTailw(fqa, 2, treschar, 2)
                fqb <- trimTailw(fqb, 2, treschar, 2)
                
                # CHECK THESE STATS COLLECTION LINES
                Q_filt_reads <- Q_filt_reads + N_reads_in - N_filt_reads - length(fqa)
                Q_filt_reads2 <- Q_filt_reads2 + N_reads_in2 - N_filt_reads2 - length(fqb)
                
                # drop reads that are less than x nt
                minLenFqa <- minLenFqa + length(fqa[width(fqa) < minlength])
                fqa <- fqa[width(fqa) >= minlength]
                minLenFqb <- minLenFqb + length(fqb[width(fqb) < minlength])
                fqb <- fqb[width(fqb) >= minlength]
                
                # ensure pair integrity is maintained
                #remove in paired reads  
                id_toKeep <- gsub(pattern = "/2",replacement = "/1",x = as.character(fqb@id))
                prRemoval <- prRemoval + length(fqa) - length(fqa[which(as.character(fqa@id) %in% id_toKeep)])#convert reads names to remove paired reads in other file (change trailing "/2" by "/1")
                fqa <- fqa[which(as.character(fqa@id) %in% id_toKeep)]
                id_toKeep <- gsub(pattern = "/1",replacement = "/2",x = as.character(fqa@id)) #convert reads names to remove paired reads in other file (change trailing "/1" to "/2")
                prRemoval2 <- prRemoval2 + length(fqb) - length(fqb[which(as.character(fqb@id) %in% id_toKeep)]) 
                fqb <- fqb[which(as.character(fqb@id) %in% id_toKeep)]
                
                # calculate number of reads output
                N_reads_out <- N_reads_out+length(fqa)
                N_reads_out2 <- N_reads_out2+length(fqb)
                
                # calculate how many reads have been trimmed
                N_trim <- N_trim + length(fqa[which(width(fqa) < init_length)])
                fullLengthReadsOut <- N_reads_out - N_trim
                N_trim2 <- N_trim2 + length(fqb[which(width(fqb) < init_length2)])
                fullLengthReadsOut2 <- N_reads_out2 - N_trim2
                
                writeFastq(object = fqa, file = destination, mode = "a",compress = T)
                writeFastq(object = fqb, file = destination2, mode = "a",compress = T)
        }
        
        # collect results into one object
        df1 <- data.frame(filtFile = destination, readsIn = N_reads_in, filterN = N_filt_reads, filterQ = Q_filt_reads, filterMinLen = minLenFqa,  readsOut = N_reads_out, trim = N_trim, fullLengthReadsOut = fullLengthReadsOut, unmatchedPair = prRemoval)        
        df2 <- data.frame(filtFile = destination2, readsIn = N_reads_in2, filterN = N_filt_reads2, filterQ = Q_filt_reads2, filterMinLen = minLenFqb,  readsOut = N_reads_out2, trim = N_trim2, fullLengthReadsOut = fullLengthReadsOut2, unmatchedPair = prRemoval2)
        
        df <- rbind(df1,df2)
        rownames(df) <- c(file, file2)
        return(df)
}

#private function
dataFileSE <- function(dataFile){
        dataFile[dataFile$PE == "SE", ]
} 

#private function        
dataFilePE <- function(dataFile){
        dataFile1 <- dataFile[ which( dataFile$PE == "PE" & grepl("_1.fastq.gz", dataFile$FILE)) , ]
        dataFile1$ID <- gsub("_1.fastq.gz", "", dataFile1$FILE)
        dataFile2 <- dataFile[ which( dataFile$PE == "PE" & grepl("_2.fastq.gz", dataFile$FILE)) , ]
        dataFile2$ID <- gsub("_2.fastq.gz", "", dataFile2$FILE)
        dataFilePE <- merge(dataFile1, dataFile2, by = "ID")
}
