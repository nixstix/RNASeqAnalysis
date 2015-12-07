#' TopHat
#'
#' @export
topHat <- function(){

    # change to working directory
    currentWD <- setwd("./")
    setwd("../../WorkingDir/")

    ### REMOVE - ADD AS ARGUMENT
    dataFile <- "data.csv"
    df <- read.csv(file = dataFile, sep = ",")
    assign("df", df, envir = .GlobalEnv)


    # create new function here , then lapply(df)
    file1 <- paste(df$FILE, "1_filtered_20000.fastq.gz", sep = "_")
    file2 <- paste(df$FILE, "2_filtered_20000.fastq.gz", sep = "_")
    file3 <- paste(df$FILE, "_filtered_20000.fastq.gz", sep = "")

    pathToFile <- as.character(df$FILE)
    dir.create(path = pathToFile, showWarnings = FALSE)
    setwd(pathToFile)

    print(paste("Running Tophat on", pathToFile, sep=" "))
    if (df$PE == "PE" && file.exists(paste("../", file1, sep="")) && file.exists(paste("../", file2, sep=""))){
        print("Paired-end alignment")
        #system("tophat2 --no-coverage-search --no-novel-juncs -r 360 -T -p 2 -o th-out --transcriptome-index ~/Documents/Personal/Computing/Thesis/Homo_sapiens_NCBI_GRCh38/TopHatTranscriptome/genes ~/Documents/Personal/Computing/Thesis/Homo_sapiens_NCBI_GRCh38/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome ../ERR315325_1_filtered_20000.fastq.gz ../ERR315325_2_filtered_20000.fastq.gz")
    } else if (df$PE == "SE" && file.exists(file3)){
        print("Single read alignment")
    } else{
        print("Data is either incorrect or files cannot be located - please check your input.")
    }

    setwd("../")

    #change back to current working directory
    setwd(currentWD)
}



