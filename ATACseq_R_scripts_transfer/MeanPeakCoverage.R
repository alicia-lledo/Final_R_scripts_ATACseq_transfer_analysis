getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/pipeline_output_all_samples_QC/correct_blood_bwa/")
getwd ()


library(ggplot2)

# Read in peak coverage data
peak.coverage.data <- read.delim("blood_bwa_coverage_peaks_250_width.800bp_window.narrowPeak", header=FALSE)
colnames(peak.coverage.data) <- c("Chromosome", "Start", "Stop", "Peak", "Base", "Coverage")

# For each position in a 400bp window around the centre of the peak, calculate the average read
# coverage across all peaks
# "by" applies a function, in this case calculating the mean read count, to subsets of a data
# frame - so the peak.coverage.data dataframe is split by row using the values in the Base column and
# the mean coverage calculated for each Base position
average.coverage.by.position <- by(peak.coverage.data, peak.coverage.data$Base, 
                                   function(x) mean(x$Coverage))

barplot(average.coverage.by.position)

average.coverage.by.position <- data.frame(vapply(average.coverage.by.position,unlist,unlist(average.coverage.by.position[[1]])))
colnames(average.coverage.by.position) <- "Average.Coverage"

plot<-ggplot(average.coverage.by.position, aes(x=1:800, y=Average.Coverage)) + geom_bar(stat="identity") +
  theme_bw() + scale_x_continuous(name="Position", breaks=c(1, seq(50, 800, 50))) + 
  coord_cartesian(xlim=c(1, 800))

pdf("monocytes_peak_coverage.pdf", width=20, height=10)
plot
dev.off ()
