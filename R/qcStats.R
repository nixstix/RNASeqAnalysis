#' Collects summary statistics on filtering
#' This function collects summary stats from the \code{filterBadSeqs} function, and prints out data to a PDF. Two graphs are output. In the first, the data shows, for each fastq read file, how many reads went in, how many went out, and the mean and median sizes of the output fastq files. The second shows the proportion of reads trimmed or removed for each fastq file.
#' 
#' @param x A dataframe returned from \code{filterBadSeqs}). 
#' @return A dataframe. The dataframe contains the number of reads which have been trimmed (on leading or trailing tails) for poor quality or unknown bases, the number of reads which have been removed entirely due to poor quality, possession of unknown bases, or length too short. Also produces graphs of absolute number of reads output, and relative numbers of reads which have been trimmed or removed (in PDF format).
#' @seealso \url{filterBadSeqs} 
#' @export
summariseFilteringStats <- function(x){
        
        # summarise stats into table
        attr(x,"names") <- NULL
        stats <- do.call(rbind, x)
        
        # prepare PDF for writing
        pdf(file = "./Barplots_filtering_stats.pdf", paper = "a4r", bg = "white")
        
        # bar plot of absolute number of reads output
        mean_reads<-mean(stats$readsOut)
        median_reads<-median(stats$readsOut)
        
        stats$diff <- stats$readsIn - stats$readsOut
        nrReadsOut <- t(stats[, c("readsOut", "diff")])
        
        barplot(height = nrReadsOut[c("readsOut", "diff"),],
                axes = TRUE, cex.names = 0.65, 
                xlab = "Sample", ylab = "Nr. reads", las = 2, names.arg = gsub(".fastq.gz", "",  rownames(stats)),
                ylim = c(0, 1.1 * max(stats[,"readsIn"])),
                main = "Number of reads output per sample", col = c("orange", "yellow"))
        abline(h = median_reads, col="blue",lwd =2, lty = 5)
        abline(h = mean_reads, col="black",lwd =2, lty = 2)
        legend(legend = c("Reads in", "Reads out", "median output", "mean output"), lty = c(1,1, 5,2) , lwd=c(5,5,2,2), 
               col = c("yellow", "orange", "blue","black"),"bottomright")
        
        # stacked bar plot of percentages of reads in, reads out, reads trimmed
        
        statsPercent <- stats[, c("readsIn", "filterN", "filterQ", "filterMinLen", "trim", "fullLengthReadsOut", "unmatchedPair")]
        statsPercent[is.na(statsPercent)] <- 0
        statsPercent <- apply(statsPercent,1,function(x) 
                x[2:7] / x[1] )
        statsPercent <- statsPercent[c(5,4,3,2,1,6), ]
        
        barplot(height = statsPercent, 
                border = NA, axes = TRUE, las = 2, cex.names = 0.65, cex.axis = 0.8,
                xlab = "Sample", ylab = "Proportion of reads", names.arg = gsub(".fastq.gz", "",  rownames(stats)),
                main = "Filtering by sample", col = c("darkgreen", "green", "red", "darkred", "indianred",  "lightpink"))
        legend(legend = rev(c("Intact reads","trimmed for trailing quality or trailing N", "removed for length", "removed for quality", "removed for N", "removed for pair mismatch")),pch = 15, col = rev(c("darkgreen", "green", "red", "darkred", "indianred", "lightpink")), "bottomright")
        
        print("A barplot of results has been saved to Barplots_filtering_stats.pdf")
        
        dev.off()
        
        # return statistics
        return(stats)
}