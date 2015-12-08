#' Filters single-end reads for quality
#' @export
#' @import ShortRead
SEFilterAndTrim <- function(file1, minlength = 30, Phred = 25){
        
        # open input stream
        file1 <- paste(file1, "fastq.gz", sep = ".")
        stream1 <- FastqStreamer(file1, readerBlockSize = 2e7,n=2e7)
        
        # define variables
        N_reads_in<-0
        N_filt_reads<-0
        Quality_filt_reads<-0
        Length_filt_reads<-0
        Intact_reads<-0
        Total_trim_reads<-0
        N_reads_out<-0
        on.exit(close(stream1))
        
        destination1 <- gsub(pattern = ".fastq.gz", replacement = "-filt.fastq.gz", x = file1)
        
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
                treschar<-names(encoding(quality(fq1))[which(encoding(quality(fq1))==Phred_tres)])
                #based on file 1 quality
                fqa <- trimTailw(fqa, 2, treschar, 2)
                Quality_filt_reads<-Quality_filt_reads+N_reads_in-N_filt_reads-length(fqa)
                
                # drop reads that are less than x nt
                Length_filt_reads<-Length_filt_reads+length(fqa[which(width(fqa)<minlength_tres)])
                Intact_reads<-Intact_reads+length(fqa[which(width(fqa)==init_length)])
                Total_trim_reads<-N_reads_in-Intact_reads-N_filt_reads-Length_filt_reads-Quality_filt_reads
                fqa <- fqa[width(fqa) >= minlength_tres]
                N_reads_out<-N_reads_out+length(fqa)
                ##     N_reads_out==N_reads_in-Nint_reads-Quality_filt_reads-Length_filt_reads
                #     Total_trim_reads==N_reads_out-Intact_reads
                #     
                ## append to destinations
                writeFastq(object = fqa, file = destination1, mode = "a",compress = T)
                
        }
        #create summary stats report
        
        if(file.exists("filteringStats")){
                dfStats <- read.table("filteringStats", sep = "\t", header = TRUE)
                df <- data.frame("File"=file1, 
                                 "No_reads_in"=N_reads_in, 
                                 "No_filt_reads"=N_filt_reads, 
                                 "Qual_filt_reads"=Quality_filt_reads, 
                                 "Len_filt_reads"=Length_filt_reads, 
                                 "No_intact_reads"=Intact_reads, 
                                 "Trimmed_reads"=Total_trim_reads, 
                                 "No_reads_out"=N_reads_out)
                dfStats <- rbind(dfStats, df)
                #write.table(dfStats, file = "filteringStats", quote = FALSE, sep = "\t")
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
                dfStats[1,] <- c(file1, N_reads_in, N_filt_reads, Quality_filt_reads, Length_filt_reads, Intact_reads, Total_trim_reads, N_reads_out)
                
        }
        write.table(dfStats, file = "filteringStats", quote = FALSE, sep = "\t", col.names = TRUE)
        dfStats <- read.table("filteringStats", sep = "\t")
        return(dfStats)
}



