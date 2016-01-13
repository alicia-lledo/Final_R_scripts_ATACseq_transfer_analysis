#Reads centered across the TSS

getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/pipeline_output_all_samples_QC/correct_blood_bwa/")
getwd ()

# Use ggplot library to draw graphs
library(ggplot2)

#Read the txt file the first column is the x-axis, while the second one is the y-axis. 
  #It should generate a histogram, peaked on the zero position. 

Reads.centered.tss<-read.delim ("tssBackgroundHist.20.txt", skip = 1, header=FALSE) #skip the first line

plot(x=Reads.centered.tss$V1, y=Reads.centered.tss$V2)


