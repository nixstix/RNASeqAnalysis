#' Performs quality checks
#' 
#' The function checks each fastq file specified in the "data file" for quality, and writes findings to a report. 
#' 
#' @param dataFile An R data frame with the data to be processed. The R object is a standard format, and must contain the following headings: File, PE, Sample, Replicate, FilteredFile. More information about the file is available at \code{\link{datafileTemplate}}. 
#' @param preFilter A logical - if true (default), the function will select and analyse files which have not yet been processed for quality. If false, the function will select and analyse those files which have been processed for quality, i.e. the "filtered file" in the data file.
#' @return FastqQA object. Outputs quality results in the form of raw data (an R FastqQA object) and HTML format (saved to "QA" directory).
#' @details The function should be run in the working directory, where all fastq files are found.
#' @details \code{\link{runQA}} iterates over each file specified in the "datafile". It runs a quality assessment from the \code{\link{ShortRead}} package. The \code{\link{ShortRead}} package (\url{https://bioconductor.org/packages/release/bioc/html/ShortRead.html}) contains more information about this step. The quality assessment may be performed before and after the filtering step, by setting the "pre-filter" parameter to true or to false, respectively. All quality assessment data is output to the "QA" directory. 
#' Quality reports are output to the working directory (under the QA directory). R objects of the raw data used to generate the reports are also saved to this directory. 
#' @export
#' @import ShortRead
runQA <- function(dataFile, preFilter = TRUE){
        

##### PARALLELISATION NOT WORKING, DISABLED FOR NOW. IF RE-INTRODUCED, DON'T FORGET TO ADD MC.CORES PARAMETER TO RUNQA FUNCTION, AND ADD BPPARAM TO QA FUNCTION        
#         if(is.null(mc.cores)){
#                 BPPARAM = registered()[1]
#         }
#         else{
#                 BPPARAM = MulticoreParam(workers = mc.cores)
#         }
        
        print("QA results will be output to the 'QA' folder")
        
        if (preFilter == TRUE){ 
                # pre-filter quality check
                print("QA on pre-filtered files")
                QASum_filter <- qa(dirPath = dataFile$FILE, type = "fastq", BPPARAM=bpparam())
                print("QA on pre-filtered files completed")
                
                QASum_prefilterRpt <- report(x = QASum_filter, dest = "QA/prefilter", type = "html")
                print("QA report and data are now available in the 'QA' folder")
                QASum_prefilter <- QASum_filter
                save(QASum_prefilter, file = "./QA/prefilter/QASum_prefilter.RData")
        } 
        
        else if(preFilter == FALSE){
                # post-filter quality check
                print("QA on post-filtered files")
                QASum_filter <- qa(dirPath = unique(dataFile$FILTEREDFILE), type = "fastq", BPPARAM = bpparam())
                print("QA on post-filtered files completed")
                QASum_postfilterRpt <- report(x = QASum_filter, dest = "QA/postfilter/", type = "html")
                print("QA report and data are now available in the 'QA' folder")
                QASum_postfilter <- QASum_filter
                save(QASum_postfilter, file = "./QA/postfilter/QASum_postfilter.RData")        
        }
        return(QASum_filter)
}