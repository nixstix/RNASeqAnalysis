#' Checks that the format of the "data file" (which is used as input in many functions) is correct.
#' @description The function takes as input a data file (see \code{datafileTemplate}), and checks whether it is in a format acceptable to be used as input to other functions.
#' @param x. A data frame containing: the name of the file to be processed, whether it is single or paried-end data, the sample and replicate ID, and (optional) the name of an output file which results from the \code{runQAandFilter} function will be written to.
#' @return A data frame which has been modified - the column names and classes have been verified (and corrected if necessary), all white space has been removed, observations relating to files which do not exist in the directory have been removed. The function also prints out any duplicated files.
#' @seealso \code{datafileTemplate}
#' @export 
checkDataFile <- function(x){
        
        # check min number of columns exists
        if (ncol(x) < 5){
                stop("You need at least 5 columns in the data file. See datafileTemplate for more information", call. = TRUE)
        }
        
        # check column headings
        x <- checkColHeadings(x)
        
        # remove white space
        cat(" Removing white space ", "\n\n")
        x$FILE <- gsub(" ", "", x$FILE)
        x$PE <- gsub(" ", "", x$PE)
        x$SAMPLE <- gsub(" ", "", x$SAMPLE)
        x$REPLICATE <- gsub(" ", "", x$REPLICATE)
        x$FILTEREDFILE <- gsub(" ", "", x$FILTEREDFILE)
        
        # convert all values to factors
        cat(" Converting columns 2, 3, and 4 to factors ; converting columns 1 and 5 to character strings ", "\n\n")
        x$PE <- as.factor(x$PE)
        x$SAMPLE <- as.factor(x$SAMPLE)
        x$REPLICATE <- as.factor(x$REPLICATE)
        x$FILE <- as.character(x$FILE)
        x$FILTEREDFILE <- as.character(x$FILTEREDFILE)
        
        # check row names are unique
        cat(" Checking for any duplicates (the user should remove these duplicates to prevent difficulties with downstream processing):", "\n")
        print(x$FILE[duplicated(x$FILE)])
        cat("\n")
        
        # check for SE/PE integrity
        cat(" Checking for any data which is not labeled as PE or as SE (this data will not be processed):", "\n")
        print(x[!grepl("^(S|P)E$", x$PE),])
        cat("\n")
       
        # check all files mentioned exist in the directory
        cat(" Checking for any files missing from the directory (this data will be REMOVED from the data set):", "\n")
        print(x[!file.exists(as.character(x$FILE)), ])
        x <- x[file.exists(as.character(x$FILE)), ]
        
        return(x)
}

#private function
checkColHeadings <- function(x){
        
        cat("\n","Checking column headers:", "\n")
        
        if (!colnames(x)[1] == "FILE"){
                colnames(x)[1] <- "FILE"
                cat("Column 1 header has been changed", "\n")
        }
        if (!colnames(x)[2] == "PE"){
                colnames(x)[2] <- "PE"
                cat("Column 2 header has been changed", "\n")
        }
        if (!colnames(x)[3] == "SAMPLE"){
                colnames(x)[3] <- "SAMPLE"
                cat("Column 3 header has been changed", "\n")
        }
        if (!colnames(x)[4] == "REPLICATE"){
                colnames(x)[4] <- "REPLICATE"
                cat("Column 4 header has been changed", "\n")
        }
        if (!colnames(x)[5] == "FILTEREDFILE"){
                colnames(x)[5] <- "FILTEREDFILE"
                cat("Column 5 header has been changed", "\n")
        }
        cat(" All columns are correctly named", "\n\n")
        
        return(x)
}
