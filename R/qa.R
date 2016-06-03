#' Performs quality checks
#' 
#' The function checks each fastq file specified in the "data file" for quality, and writes findings to a report. 
#' 
#' @param dataFile An R data frame with the data to be processed. The R object is a standard format, and must contain the following headings: File, PE, Sample, Replicate, FilteredFile. More information about the file is available at \code{datafileTemplate}. 
#' @param preFilter A logical - if true (default), the function will select and analyse files which have not yet been processed for quality. If false, the function will select and analyse those files which have been processed for quality, i.e. the "filtered file" in the data file.
#' @return Outputs directories with quality results in the form of raw data and HTLM format
#' @details The function should be run in the working directory, where all fastq files are found.
#' @details \code{runQA} iterates over each file specified in the "datafile". It runs a quality assessment from the \code{ShortRead} package. The \code{ShortRead} package (\url{https://bioconductor.org/packages/release/bioc/html/ShortRead.html}) contains more information about this step. The quality assessment may be performed before and after the filtering step, by setting the "pre-filter" parameter to true or to false, respectively. All quality assessment data is output to the "QA" directory. 
#' Quality reports are output to the working directory (under the QA directory). R objects of the raw data used to generate the reports are also saved to this directory. 
#' @export
#' @import ShortRead
runQA <- function(dataFile, preFilter = TRUE){
        print("QA results will be output to the 'QA' folder")
        
        if (preFilter == TRUE){ 
                # pre-filter quality check
                QASum_prefilter <- qa(dirPath = dataFile$FILE, type = "fastq")
                QASum_prefilterRpt <- report(x = QASum_prefilter, dest = "QA/prefilter", type = "html")
                save(QASum_prefilter, file = "./QA/prefilter/QASum_prefilter.RData")
        } 
        
        else if(preFilter == FALSE){
                # pre-filter quality check
                QASum_postfilter <- qa(dirPath = unique(dataFile$FILTEREDFILE), type = "fastq")
                QASum_postfilterRpt <- report(x = QASum_postfilter, dest = "QA/postfilter/", type = "html")
                save(QASum_postfilter, file = "./QA/postfilter/QASum_postfilter.RData")        
        }
}