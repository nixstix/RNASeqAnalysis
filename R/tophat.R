 rm(list = ls())

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
        TranscriptomeIndex <- parameterFile[c("TranscriptomeIndex"), 1]
        Genome <- parameterFile[c("Genome"), 1]
        
        # run tophat
        tophatRun(df = df, mateInnerDist = mateInnerDist, threads = threads, Genome = Genome,  TranscriptomeIndex = TranscriptomeIndex)
   
        
        #return to original directory
        setwd(currentWD)
}

tophatRun <- function(df, mateInnerDist, threads, Genome, TranscriptomeIndex){
        
                        
        
        file1 <- paste(df$FILE, "1_filtered_20000.fastq.gz", sep = "_")
        file2 <- paste(df$FILE, "2_filtered_20000.fastq.gz", sep = "_")
        file3 <- paste(df$FILE, "_filtered_20000.fastq.gz", sep = "")
        
        pathToFile <- as.character(df$FILE)
        dir.create(path = pathToFile, showWarnings = FALSE)
        setwd(pathToFile)
        
        print(paste("Running Tophat on", pathToFile, sep=" "))
        if (df$PE == "PE" && file.exists(paste("../", file1, sep="")) && file.exists(paste("../", file2, sep=""))){
                
                
                print("Paired-end alignment")
                file1 <- paste("../", file1, sep="")   
                file2 <- paste("../", file2, sep="")
                thCmd <- paste("tophat2 --no-coverage-search --no-novel-juncs -r", mateInnerDist , "-T -p", threads, "-o th-out --transcriptome-index", TranscriptomeIndex, Genome, file1, file2, sep=" ")
                
        } else if (df$PE == "SE" && file.exists(file3)){
                print("Single read alignment")
                file3 <- paste("../", file3, sep="")
                thCmd <- paste("tophat2 --no-coverage-search --no-novel-juncs, -T -p", threads, "-o th-out --transcriptome-index", TranscriptomeIndex, Genome, file3, sep=" ")
                
        } else{
                
                stop("Data is either incorrect or files cannot be located - please check your input.", call. = FALSE)
        }
        
        system(thCmd)
        print("Cufflinks quantification")
        TranscriptomeIndex <- paste(TranscriptomeIndex, "gff", sep = ".")
        clCmd <- paste("cufflinks -o cl-out -p", threads, "-G", TranscriptomeIndex, "th-out/accepted_hits.bam", sep = " ")
        system(clCmd)
        
        #### add th-out to assemblies.txt
        
        setwd("../")
}

tophat("/Users/thrupp/Desktop/Data&Experiments/RPackage/WorkingDir/")
