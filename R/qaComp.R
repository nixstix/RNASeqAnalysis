#' Compares quality before and after filtering fastq files.
#' Porvides a graphic of the change in quality of fastq files, before and after running the \code{\link{filterBadSeqs}} function.
#' @param dataFile An R data frame with the data to be processed. The R object is a standard format, and must contain the following headings: File, PE, Sample, Replicate, FilteredFile. More information about the file is available at \code{datafileTemplate}. The \code{\link{checkDataFile}} function should be run on the file before use.
#' @param prefiltData A FastqQA object resulting from \code{\link{runQA}}. The object contains QA data on fastq files before the filtering procedure (\code{\link{filterBadSeqs}}).
#' @param postfiltData A FastqQA object resulting from \code{\link{runQA}}. The object contains QA data on fastq files after the filtering procedure (\code{\link{filterBadSeqs}}). 
#' @return a datafrome containing densities and quality scores for each fastq file (both before and after filtering). It also outputs a graph which shows comparisons of quality before and after filtering. 
#' @import lattice 
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline barplot legend
#' @importFrom stats end median start
#' @export
compareQA <- function(dataFile, prefiltData, postfiltData){
        
        dfPre <- prefiltData[["readQualityScore"]]
        dfPost <- postfiltData[["readQualityScore"]]
        
        lsQA <- lapply(1:nrow(dataFile), getResults, dataFile, dfPre, dfPost)
        dfQA <- do.call(rbind, lsQA)
        
        # prepare PDF for writing
        print("Results are graphically represented in the plotQAafterFiltering.pdf document")
        pdf(file = "./plotQAafterFiltering.pdf", paper = "a4r", bg = "white")
        plt <- xyplot(density ~ quality | factor, dfQA, auto.key = T, type = "l", lwd = 2, groups = run, layout = c(4,round(nrow(dataFile)/3)))
        print(plt)
        dev.off()
        
        return(dfQA)
        
}

# private function
getResults <- function(x, dataFile, data1 = dfPre, data2 = dfPost){
        
        preMatch <- dataFile[x, "FILE"]
        postMatch <- dataFile[x,"FILTEREDFILE"]
        
        data1 <- data1[data1$lane == preMatch, ]
        data1$run <- factor("Pre-filtering")
        data2 <- data2[data2$lane == postMatch, ]
        data2$run <- factor("Post-filtering")
        dfAll <- rbind(data1, data2)
        dfAll$factor <- gsub("fastq.gz", "", postMatch) 
        return(dfAll)
        
}
