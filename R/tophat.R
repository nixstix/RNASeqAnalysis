
#' Collects stats from Tophat and Cufflinks runs 
#' 
#' DESCRIPTION GOES HERE. MENTION DATA FILE
#' @param WD A character string that is the path to the working directory, where the data file (containing summary information of the files that were processed by Tophat and Cufflinks) and all output fromTophat and Cufflinks can be found. The default is the current working directory.
#' @return A data frame containing statistics of all files that were processed by Tophat and Cufflinks - the mimimum and maximum lengths of the reads, the number (and percentages) of reads aligned, the number of genes and isoforms identified, and the number of genes and isoforms for which the FPKM was successfully calculated.    
#' @examples 
#' tophatStats()
#' tophatStats("path/to/working/directory")
#' @seealso \code{\link{tophat}} runs alignment and expression calculations on sequence data.
#' @export
tophatStats <- function(WD = "./"){
        
        # change to working directory, open data file                
        currentWD <- setwd("./")
        setwd(WD)
        dataFile <- read.csv(file = "data.csv", sep = ",")
        
        # collect stats for each line of data file
        dfStat <- apply(dataFile, 1, tophatStatsRun)
        
        # return to current directory
        setwd(currentWD)
        return(dfStat)
}

tophatStatsRun <- function(dataFile){
        
        #open Tophat output
        fileName <- as.character(dataFile[1])
        prepReads <- read.table(paste(fileName, "/th-out/prep_reads.info", sep=""), sep = "=")
        
        # if single end data,
        # collect mimimum read length, maximum read length, number and percentage of reads mapped
        if (dataFile[2] == "SE" && file.exists(fileName)){
                min_read_len <- as.numeric(prepReads[1, 2])
                max_read_len <- as.numeric(prepReads[2, 2])
                no_reads_mapped <- as.numeric(prepReads[4,2])
                percent_reads_mapped <- as.numeric(prepReads[4,2]/prepReads[3,2]*100)
                percent_reads_mapped <- round(percent_reads_mapped,digits=1)
                
                # if paired end data,
                # collect mimimum read length (lowest value of left or right), 
                # maximum read length (highest value of left or right), 
                # number and percentage of reads mapped (average of left and right)      
        } else if (dataFile[2] == "PE" && file.exists(fileName)){
                min_read_len <- as.numeric(min(prepReads[1, 2], prepReads[5, 2]))
                max_read_len <- as.numeric(min(prepReads[2, 2], prepReads[6, 2]))
                no_reads_mapped <- as.numeric(mean(prepReads[4, 2], prepReads[8, 2]))
                percent_reads_mapped <- as.numeric(no_reads_mapped/mean(prepReads[3,2], prepReads[7,2])*100)
                percent_reads_mapped <- round(percent_reads_mapped,digits=1)
        }
        
        #collect Cufflinks stats (same for paired end and single end data)
        # number of genes counted, plus how many genes where FPKM status is "OK"
        # number of isoforms counted, plus how many isoforms where FPKM status is "OK"
        geneCount <- read.table(paste(fileName, "/cl-out/genes.fpkm_tracking", sep=""), sep = "\t", header=TRUE)
        no_genes <- as.numeric(length(geneCount[,1]))
        no_genes_ok <- as.numeric(length(no_genes[geneCount$FPKM_status == "OK"]))
        
        isoformCount <- read.table(paste(fileName, "/cl-out/isoforms.fpkm_tracking", sep=""), sep = "\t", header=TRUE)
        no_isoform <- as.numeric(length(isoformCount[,1]))
        no_isoform_ok <- as.numeric(length(no_isoform[isoformCount$FPKM_status == "OK"]))
        
        # if stats file already exists, append data (saves overwriting if user mistakenly runs function twice with same data)
        if(file.exists("tophatStats")){
                dfStats <- read.table("tophatStats", sep = "\t", header = TRUE)
                dfStatsCurrent <- data.frame("File"=fileName, 
                                             "Min_read_len"=min_read_len, 
                                             "Max_read_len"=max_read_len, 
                                             "No_reads_mapped"=no_reads_mapped,
                                             "Percent_reads_mapped"=percent_reads_mapped, 
                                             "No_genes"=no_genes, 
                                             "Gene_FPKM_stat_OK"=no_genes_ok, 
                                             "No_isoforms"=no_isoform,
                                             "Isoform_FPKM_stat_OK"=no_isoform_ok)
                dfStats <- rbind(dfStats, dfStatsCurrent)
                
                # if stats data files doesn't exist, create new file, append data
        } else {
                dfStats <- data.frame(
                        File=character(),
                        Min_read_len=numeric(),
                        Max_read_len=numeric(),
                        No_reads_mapped=numeric(),
                        Percent_reads_mapped=numeric(),
                        No_genes=numeric(),
                        Gene_FPKM_stat_OK=numeric(),
                        No_isoforms=numeric(),
                        Isoform_FPKM_stat_OK=numeric(),
                        stringsAsFactors = FALSE
                )
                dfStats[1,] <- c(File = fileName, Min_read_len=min_read_len, Max_read_len=max_read_len, No_reads_mapped=no_reads_mapped, Percent_reads_mapped=percent_reads_mapped, No_genes=no_genes, Gene_FPKM_stat_OK=no_genes_ok, No_isoforms=no_isoform, Isoform_FPKM_stat_ok=no_isoform)
        }
        
        #write data to external table
        write.table(dfStats, file = "tophatStats", quote = FALSE, sep = "\t", col.names = TRUE)
        dfStats <- read.table("tophatStats", sep = "\t", header = TRUE)
        return(dfStats)
}

#' TopHat
#'
#' @export
tophat <- function(WD = "./"){
        
        # change to working directory
        currentWD <- setwd("./")
        setwd(WD)
        
        #open necessary files, define variables
        dataFile <- "data.csv"
        df <- read.csv(file = dataFile, sep = ",")
        
        parameterFile <- read.table(file = "parameters.txt", sep = "=", row.names = "V1")
        mateInnerDist <- parameterFile[c("MateInnerDist"), 1]        
        threads <- parameterFile[c("Threads"), 1]
        transcriptome <- parameterFile[c("TranscriptomeIndex"), 1]
        Genome <- parameterFile[c("Genome"), 1]
        
        transcriptomeIndex <- buildIndex(transcriptome, Genome)
        
        apply(df, 1, tophatRun, mateInnerDist = mateInnerDist, threads = threads, Genome = Genome,  transcriptomeIndex = transcriptomeIndex)
        
        #return to original directory
        setwd(currentWD)
}

tophatRun <- function(df, mateInnerDist, threads, Genome, transcriptomeIndex){
        
        file1 <- paste(df[1], "1.fastq.gz", sep = "_")
        file2 <- paste(df[1], "2.fastq.gz", sep = "_")
        file3 <- paste(df[1], ".fastq.gz", sep = "")
        
        pathToFile <- as.character(df[1])
        dir.create(path = pathToFile, showWarnings = FALSE)
        setwd(pathToFile)
        
        cat("\n ******************** \n")
        cat(paste("Running Tophat on", pathToFile, "\n", sep=" "))
        cat("******************** \n")
        
        
        
        if (df[2] == "PE" && file.exists(paste("../", file1, sep = "")) && file.exists(paste("../", file2, sep = ""))){
                
                cat("Paired-end alignment \n")
                file1 <- paste("../", file1, sep="")   
                file2 <- paste("../", file2, sep="")
                thCmd <- paste("tophat --no-coverage-search --no-novel-juncs -r", mateInnerDist , "-T -p", threads, "-o th-out --transcriptome-index", transcriptomeIndex, Genome, file1, file2, sep=" ")
                system(thCmd)
                
        } else if (df[2] == "SE" && file.exists(paste("../", file3, sep = ""))){
                cat("Single read alignment \n")
                file3 <- paste("../", file3, sep="")
                thCmd <- paste("tophat --no-coverage-search --no-novel-juncs -T -p", threads, "-o th-out --transcriptome-index", transcriptomeIndex, Genome, file3, sep=" ")
                system(thCmd)
                
        }
        else {
                cat("Incorrect data entered - check your data file \n")
                
        }
        
        # run Cufflinks
        
        if(file.exists("th-out/accepted_hits.bam")){
                cat("\n Cufflinks quantification \n \n")
        clCmd <- paste("cufflinks -o cl-out -p", threads, "-G", transcriptomeIndex, "th-out/accepted_hits.bam", sep = " ")
        system(clCmd)
        } else {
                cat("\n Tophat output not found - cannot perform Cufflinks run \n \n")
        }
        
        #### add th-out to assemblies.txt
        
        setwd("../")
}        

buildIndex <- function(transcriptome, Genome){
        if(grepl(".gtf", transcriptome) == TRUE){
                cat("Building transcriptome Index \n")
                buildIndexCmd <- paste("tophat2 -G", transcriptome, "--transcriptome-index transcriptome", Genome, sep = " ")
                system(buildIndexCmd)
                TIpath <- strsplit(as.character(transcriptome), "/")
                TI <- TI[[1]][length(TI[[1]])]
                TI <- paste("transcriptome/", TI, sep = "")
                transcriptomeIndex <- gsub(pattern = "gtf", replacement = "gff", TI)
        } else if(grepl(".gff", transcriptome) == TRUE){
                transcriptomeIndex <- transcriptome 
        }
        return(transcriptomeIndex)
}
