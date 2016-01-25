
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

#' private function
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

