#' @name datafileTemplate
#' @title (Empty) template data file which summarises the data to be analysed
#' @docType data
#' @author Nicky Thrupp \email{nixthrupp@gmail.com}
#' @keywords datasets
#' @description This is an empty dataset, the purpose of which is to act as a template for users. In the data file, the name of the file, the name of the sample and the technical replicate, and whether the sample is paired-end or single-end should be included. In addition, if the user is planning to filter the input file (see \code{runQAandFilter}), a column is provided for the user to add the name of the output fastq file.
#' @description This data file acts as input to many of the package's functions. The format of the data file should not be altered, save to add observations (rows).
#' @description The data set \code{datafileExample} is an example of a data file containing observations.  
#' @details \itemize{
#' \item File. The name of the fastq (fastq.gz) file to be processed
#' \item PE. "PE" should be entered into this column if the sample is from a paired-end run, or "SE" if the sample is from a single-end run.
#' \item Sample. The name of the sample (e.g. "control", "lung").
#' \item Replicate. The name of the replicate (e.g. "a", "2").
#' \item FilteredFile. The name of the fastq (fastq.gz) file that the results from \code{runQAandFilter} function are written to. 
#' }
#' @usage datafileTemplate 
#' @format a dataframe with 0 observations of 5 vairables
NULL


#' @name datafileExample
#' @title Example data file which summarises the data to be analysed
#' @docType data
#' @author Nicky Thrupp \email{nixthrupp@gmail.com}
#' @keywords datasets
#' @description This is a dataset, the purpose of which is to demonstrate the types of information that should be stored in such a data file. In the data file, the name of the file, the name of the sample and the technical replicate, and whether the sample is paired-end or single-end should be included. In addition, if the user is planning to filter the input file (see \code{runQAandFilter}), a column is provided for the user to add the name of the output fastq file.
#' @description This data file acts as input to many of the package's functions. The format of the data file should not be altered, save to add observations (rows).
#' @description The data set \code{datafileExample} is an example of a data file containing observations.  
#' @details \itemize{
#' \item File. The name of the fastq (fastq.gz) file to be processed
#' \item PE. "PE" should be entered into this column if the sample is from a paired-end run, or "SE" if the sample is from a single-end run.
#' \item Sample. The name of the sample (e.g. "control", "lung").
#' \item Replicate. The name of the replicate (e.g. "a", "2").
#' \item FilteredFile. The name of the fastq (fastq.gz) file that the results from \code{runQAandFilter} function are written to. 
#' }
#' @usage datafileTemplate 
#' @format a dataframe with 6 observations of 5 vairables
NULL