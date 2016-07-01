#' Performs quality checks
#' 
#' The function checks each fastq file specified in the "data file" for quality, and writes findings to a report. 
#' 
#' @param dataFile An R data frame with the data to be processed. The R object is a standard format, and must contain the following headings: File, PE, Sample, Replicate, FilteredFile. More information about the file is available at \code{datafileTemplate}. 
#' @param preFilter A logical - if true (default), the function will select and analyse files which have not yet been processed for quality. If false, the function will select and analyse those files which have been processed for quality, i.e. the "filtered file" in the data file.
#' @param mc.cores A number specifying the number of cores to be used during parallel processing. 
#' @return FastqQA object. Outputs quality results in the form of raw data (an R FastqQA object) and HTML format (saved to "QA" directory).
#' @details The function should be run in the working directory, where all fastq files are found.
#' @details \code{runQA} iterates over each file specified in the "datafile". It runs a quality assessment from the \code{ShortRead} package. The \code{ShortRead} package (\url{https://bioconductor.org/packages/release/bioc/html/ShortRead.html}) contains more information about this step. The quality assessment may be performed before and after the filtering step, by setting the "pre-filter" parameter to true or to false, respectively. All quality assessment data is output to the "QA" directory. 
#' Quality reports are output to the working directory (under the QA directory). R objects of the raw data used to generate the reports are also saved to this directory. 
#' @export
#' @import ShortRead
runQA <- function(dataFile, preFilter = TRUE, mc.cores = NULL){
        
        if(!missing(mc.cores)){
                BPPARAM = MulticoreParam(workers = mc.cores)
        }
        else{
                BPPARAM = registered()[1]
        }
        
        print("QA results will be output to the 'QA' folder")
        
        if (preFilter == TRUE){ 
                # pre-filter quality check
                print("QA on pre-filtered files")
                QASum_filter <- qa(dirPath = dataFile$FILE, type = "fastq", BPPARAM = BPPARAM)
                print("QA on pre-filtered files completed")
                QASum_prefilterRpt <- report(x = QASum_filter, dest = "QA/prefilter", type = "html")
                print("QA report and data are now available in the 'QA' folder")
                save(QASum_filter, file = "./QA/prefilter/QASum_prefilter.RData")
        } 
        
        else if(preFilter == FALSE){
                # post-filter quality check
                print("QA on post-filtered files")
                QASum_filter <- qa(dirPath = unique(dataFile$FILTEREDFILE), type = "fastq", BPPARAM = BPPARAM)
                print("QA on post-filtered files completed")
                QASum_postfilterRpt <- report(x = QASum_filter, dest = "QA/postfilter/", type = "html")
                print("QA report and data are now available in the 'QA' folder")
                save(QASum_filter, file = "./QA/postfilter/QASum_postfilter.RData")        
        }
        return(QASum_filter)
}