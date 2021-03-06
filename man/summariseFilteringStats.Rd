% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qcStats.R
\name{summariseFilteringStats}
\alias{summariseFilteringStats}
\title{Collects summary statistics on filtering
This function collects summary stats from the \code{filterBadSeqs} function, and prints out data to a PDF. Two graphs are output. In the first, the data shows, for each fastq read file, how many reads went in, how many went out, and the mean and median sizes of the output fastq files. The second shows the proportion of reads trimmed or removed for each fastq file.}
\usage{
summariseFilteringStats(x)
}
\arguments{
\item{x}{A dataframe returned from \code{filterBadSeqs}).}
}
\value{
A dataframe. The dataframe contains the number of reads which have been trimmed (on leading or trailing tails) for poor quality or unknown bases, the number of reads which have been removed entirely due to poor quality, possession of unknown bases, or length too short. Also produces graphs of absolute number of reads output, and relative numbers of reads which have been trimmed or removed (in PDF format).
}
\description{
Collects summary statistics on filtering
This function collects summary stats from the \code{filterBadSeqs} function, and prints out data to a PDF. Two graphs are output. In the first, the data shows, for each fastq read file, how many reads went in, how many went out, and the mean and median sizes of the output fastq files. The second shows the proportion of reads trimmed or removed for each fastq file.
}
\seealso{
\url{filterBadSeqs}
}

