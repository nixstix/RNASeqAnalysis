#' Filters single-end reads for quality
#' 
#' DESCRIPTION GOES HERE
#' @param file A character string that is the name of the file to be processed. The input file must be a zipped fastq file (fastq.gz), but this extension should not be included in the character string.
#' @param minLength An integer which specifies the minimum length for a read. Reads shorter than this length will be discarded. Default is 30 nucleotides.
#' @param Phred An integer which specifies Phred (ascii) quality score. Any two consecutive nucleotides with a quality score lower than this threshold will be discarded. Default score is 30.
#' @return A data frame containing statistics on the QC - how many reads were processed, how many reads were removed, the length of the reads. This information is appended to a "filteringStats" document in the working directory. In addition, the function generates a filtered fastq.gz file ("file-filt.fastq.gz").
#' @examples 
#' SEFilterAndTrim(file = "reads")
#' SEFilterAndTrim(file = "reads", minlength = 35, Phred = 30)
#' SEFilterAndTrim(file = "reads", Phred = 30)
#' @seealso \url{https://en.wikipedia.org/wiki/Phred_quality_score} for more about quality scores. \code{\link{PEFilterAndTrim}} processes paired-end data in a similar manner.
#' @export
#' @import ShortRead
SEFilterAndTrim <- function(file, minlength = 30, Phred = 25){
        
        # open input stream
        file <- paste(file, "fastq.gz", sep = ".")
        stream1 <- openStream(file = file)
        
        # define variables
        N_reads_in<-0
        N_filt_reads<-0
        Quality_filt_reads<-0
        Length_filt_reads<-0
        Intact_reads<-0
        Total_trim_reads<-0
        N_reads_out<-0
        on.exit(close(stream1))
        
        destination1 <- gsub(pattern = ".fastq.gz", replacement = "-filt.fastq.gz", x = file)
        
        repeat {
                # input chunk
                
                fq1 <- yield(stream1)
                fq1Len<-length(fq1)
                if (fq1Len == 0)
                        break
                N_reads_in<-N_reads_in+fq1Len
                init_length<-max(width(fq1))
                
                
                #trim leading/trailing N
                pos_1 <- trimEnds(sread(fq1), "N", relation="==", ranges=T)
                fqa<-narrow(fq1, start(pos_1), end(pos_1))
                ##     Ntrim_reads<-Ntrim_reads+length(fqa[which(width(fqa)<init_length)])
                ##     Ntrim_reads_bl_tresh<-Ntrim_reads_bl_tresh+length(fqa[which(width(fqa)<minlength_tres)])
                
                # filter reads containing N
                failed_N_1<-nFilter()(fqa)
                toTrash_1<-which(failed_N_1@.Data==F)
                N_filt_reads<-N_filt_reads+length(toTrash_1)
                if(N_filt_reads>0) {
                        fqa<-fqa[-toTrash_1]
                }
                
                # trim as soon as 2 of 7 nucleotides has quality encoding less than phred score tres
                #determine quality ascii character for trimming
                treschar<-names(encoding(quality(fq1))[which(encoding(quality(fq1))==Phred)])
                #based on file 1 quality
                fqa <- trimTailw(fqa, 2, treschar, 2)
                Quality_filt_reads<-Quality_filt_reads+N_reads_in-N_filt_reads-length(fqa)
                
                # drop reads that are less than x nt
                Length_filt_reads<-Length_filt_reads+length(fqa[which(width(fqa)<minlength)])
                Intact_reads<-Intact_reads+length(fqa[which(width(fqa)==init_length)])
                Total_trim_reads<-N_reads_in-Intact_reads-N_filt_reads-Length_filt_reads-Quality_filt_reads
                fqa <- fqa[width(fqa) >= minlength]
                N_reads_out<-N_reads_out+length(fqa)
                ##     N_reads_out==N_reads_in-Nint_reads-Quality_filt_reads-Length_filt_reads
                #     Total_trim_reads==N_reads_out-Intact_reads
                #     
                ## write to destinations
                writeFastq(object = fqa, file = destination1, mode = "a",compress = T)
                
        }
        #create summary stats report
        
        if(file.exists("filteringStatsSE")){
                dfStats <- read.table("filteringStatsSE", sep = "\t", header = TRUE)
                df <- data.frame("File"=file, 
                                 "No_reads_in"=N_reads_in, 
                                 "No_filt_reads"=N_filt_reads, 
                                 "Qual_filt_reads"=Quality_filt_reads, 
                                 "Len_filt_reads"=Length_filt_reads, 
                                 "No_intact_reads"=Intact_reads, 
                                 "Trimmed_reads"=Total_trim_reads, 
                                 "No_reads_out"=N_reads_out)
                dfStats <- rbind(dfStats, df)
                
        } else {
                
                dfStats <- data.frame(
                        File=character(),
                        No_reads_in=numeric(),
                        No_filt_reads=numeric(),
                        Qual_filt_reads=numeric(),
                        Len_filt_reads=numeric(),
                        No_intact_reads=numeric(),
                        Trimmed_reads=numeric(),
                        No_reads_out=numeric(),
                        stringsAsFactors=FALSE
                )
                dfStats[1,] <- c(file, N_reads_in, N_filt_reads, Quality_filt_reads, Length_filt_reads, Intact_reads, Total_trim_reads, N_reads_out)
                
        }
        write.table(dfStats, file = "filteringStatsSE", quote = FALSE, sep = "\t", col.names = TRUE)
        dfStats <- read.table("filteringStatsSE", sep = "\t")
        return(dfStats)
}


#*****************************************************************************

#' Filters single-end reads for quality
#' 
#' DESCRIPTION GOES HERE
#' @param file A character string that is the name of the file to be processed. Files should be provided in pairs, with the extensions "_1.fastq.gz" and "_2.fastq.gz" - however, these extensions should not be included in the character string. 
#' @param minLength An integer which specifies the minimum length for a read. Reads shorter than this length will be discarded. Default is 30 nucleotides.
#' @param Phred An integer which specifies Phred (ascii) quality score. Any two consecutive nucleotides with a quality score lower than this threshold will be discarded. Default score is 25.
#' @return A data frame containing statistics on the QC - how many reads were processed, how many reads were removed, the length of the reads. This information is appended to a "filteringStats" document in the working directory. In addition, the function generates a filtered fastq.gz file ("file_1-filt.fastq.gz", "file_2-filt.fastq.gz").
#' @examples 
#' SEFilterAndTrim(file = "reads")
#' SEFilterAndTrim(file = "reads", minlength = 25, Phred = 20)
#' SEFilterAndTrim(file = "reads", Phred = 40)
#' @seealso \url{https://en.wikipedia.org/wiki/Phred_quality_score} for more about quality scores. \code{\link{PEFilterAndTrim}} processes paired-end data in a similar manner.
#' @export
#' @import ShortRead
PEFilterAndTrim <- function(file, minlength = 30, Phred = 25)  
{
        
        ## open input stream
        file1 <- paste(file, "_1.fastq.gz", sep = "")
        stream1 <- openStream(file = file1)
        file2 <- paste(file, "_2.fastq.gz", sep = "")
        stream1 <- openStream(file = file2)
        
        # define variables
        N_reads_in_1<-0
        N_reads_in_2<-0
        
        repeat {
                ## input chunk
                fq1 <- yield(stream1)
                fq1Len<-length(fq1)
                if (fq1Len == 0)
                        break
                N_reads_in_1<-N_reads_in_1+fq1Len
                init_length_1<-max(width(fq1))
                
                fq2 <- yield(stream1)
                fq2Len<-length(fq2)
                if (fq2Len == 0)
                        break
                N_reads_in_2<-N_reads_in_2+fq2Len
                init_length_2<-max(width(fq2))
                
                ##trim leading/trailing N
                pos_1 <- trimEnds(sread(fq1), "N", relation="==", ranges=T)
                fqa<-narrow(fq1, start(pos_1), end(pos_1))
                pos_2<- trimEnds(sread(fq2), "N", relation="==", ranges=T)
                fqb<-narrow(fq2, start(pos_2), end(pos_2))
                ## filter reads containing N
                failed_N_1<-nFilter()(fqa)
                toTrash_1<-which(failed_N_1@.Data==F)
                failed_N_2<-nFilter()(fqb)
                toTrash_2<-which(failed_N_2@.Data==F)
                toTrash<-unique(union(toTrash_1,toTrash_2))
                if(length(toTrash)>0) {
                        fqa<-fqa[-toTrash]
                        fqb<-fqb[-toTrash]
                }
                
                ## trim as soon as 2 of 7 nucleotides has quality encoding less than phred score tres
                ##determine quality ascii character for trimming
                treschar<-names(encoding(quality(fq1))[which(encoding(quality(fq1))==Phred)])
                #based on file 1 quality
                fqa <- trimTailw(fqa, 2, treschar, 2)
                #remove in paired reads  
                id_toKeep<-gsub(pattern = "/1",replacement = "/2",x = as.character(fqa@id)) #convert reads names to remove paired reads in other file (change trailing "/1" by "/2")
                fqb<-fqb[which(as.character(fqb@id) %in% id_toKeep)]
                #based on file 2 quality
                fqb <- trimTailw(fqb, 2, treschar, 2)
                #remove in paired reads  
                id_toKeep<-gsub(pattern = "/2",replacement = "/1",x = as.character(fqb@id)) #convert reads names to remove paired reads in other file (change trailing "/2" by "/1")
                fqa<-fqa[which(as.character(fqa@id) %in% id_toKeep)]
                ## drop reads that are less than 30nt
                fqa <- fqa[width(fqa) >= minlength]
                #remove in paired reads  
                id_toKeep<-gsub(pattern = "/1",replacement = "/2",x = as.character(fqa@id)) #convert reads names to remove paired reads in other file (change trailing "/1" by "/2")
                fqb<-fqb[which(as.character(fqb@id) %in% id_toKeep)]
                
                fqb <- fqb[width(fqb) >= minlength]
                #remove in paired reads  
                id_toKeep<-gsub(pattern = "/2",replacement = "/1",x = as.character(fqb@id)) #convert reads names to remove paired reads in other file (change trailing "/2" by "/1")
                fqa<-fqa[which(as.character(fqa@id) %in% id_toKeep)]
                
                ## write to destination
                destination_1 <- paste(file, "-filt_1.fastq.gz", sep = "")
                writeFastq(object = fqa, file = destination_1, mode = "w",compress = T)
                destination_2 <- paste(file, "-filt_2.fastq.gz", sep = "")
                writeFastq(object = fqb, file = destination_2, mode = "w",compress = T)
        }
        
        # collect stats
        
        if(file.exists("filteringStatsPE")){
                dfStats <- read.table("filteringStatsPE", sep = "\t", header = TRUE)
                df <- data.frame("File"=file, 
                                 "No_reads_in"=N_reads_in_1 )
                dfStats <- rbind(dfStats, df)
                
        } else {
                
                dfStats <- data.frame(
                        File=character(),
                        No_reads_in=numeric(),
                        #No_filt_reads=numeric(),
                        #Qual_filt_reads=numeric(),
                        #Len_filt_reads=numeric(),
                        #No_intact_reads=numeric(),
                        #Trimmed_reads=numeric(),
                        #No_reads_out=numeric(),
                        stringsAsFactors=FALSE
                )
                dfStats[1,] <- c(file1, N_reads_in_1)
                
        }
        write.table(dfStats, file = "filteringStatsPE", quote = FALSE, sep = "\t", col.names = TRUE)
        dfStats <- read.table("filteringStatsPE", sep = "\t")
        return(dfStats)
}

# work on stats to be returned

# private function
openStream <- function(file){
        stream <- FastqStreamer(file, readerBlockSize = 2e7,n=2e7)
        on.exit(close(stream))
        return(stream)
}
