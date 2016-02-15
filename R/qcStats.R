#' Collects summary statistics on filtering
#' 
#' DESCRIPTION GOES HERE
#' 
#' @param x An object of class QualityFilterResults (e.g. returned from \code{runQAandFilter}). 
#' @return A dataframe. The dataframe contains the number of reads which have been trimmed (on leading or trailing tails) for poor quality or unknown bases, the number of reads which have been removed entirely due to poort quality, possession of unknown bases, or length too short. Also produces graphs of absolute number of reads output, and relative numbers of reads which have been trimmed or removed (in PDF format).
#' @seealso \url{runQAandFilter} 
#' @export
summariseFilteringStats <- function(x){
        
        # summarise stats into table
        stats <- do.call(rbind, lapply(x, attr, "filter"))
        rownames(stats) <- do.call(rbind, lapply(x, attr, "name"))
        
        # prepare PDF for writing
        pdf(file = "./Barplots_filtering_stats.pdf",paper = "a4r",bg = "white")
        
        # bar plot of absolute number of reads output
        mean_reads<-mean(stats$readsOut)
        median_reads<-median(stats$readsOut)
        statsMat <- as.matrix(t(stats))
        barplot(height = statsMat[4,], 
                border=NA,
                axes = TRUE, cex.names = 0.65, 
                xlab = "Sample", ylab = "Nr. reads", las = 2,
                main = "Number of reads output per sample", col = "darkred")
        abline(h = median_reads, col="blue",lwd =2, lty = 5)
        abline(h = mean_reads, col="black",lwd =2, lty = 2)
        legend(legend = c("median", "mean"), lty = c(6,2,1), lwd=2, 
               col = c("blue","black"),"bottomright")
        
        # stacked bar plot of percentages of reads in, reads out, reads trimmed
        statsPercent <- stats[ , c(1:3, 5:7)]
        statsPercent <- apply(statsPercent,1,function(x) 
                x[2:6]/rep(x[1],length(x[2:6])))
        statsPercent <- statsPercent[c(5,4,3,2,1), ]
        
        barplot(height = statsPercent, 
                border = NA, axes = TRUE, las = 2, cex.names = 0.65, cex.axis = 0.8,
                xlab = "Sample", ylab = "Proportion of reads",
                main = "Filtering by sample", col = c("darkgreen", "green", "lightgreen", "red", "darkred"))
        legend(legend = rev(c("Intact reads","trimmed for trailing quality","trimmed for trailing N", "removed for length", "removed for N")),pch = 15, col = rev(c("darkgreen", "green","lightgreen", "red", "darkred")), "bottomright")
        
        print("A barplot of results has been saved to Barplots_filtering_stats.pdf")
        
        dev.off()
        
        # return statistics
        return(stats)
}
