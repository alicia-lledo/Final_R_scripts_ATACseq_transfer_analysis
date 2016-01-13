getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots")
getwd ()

# Use ggplot library to draw graphs
library(ggplot2)


###NORMALISATION BASED ON CHROMOSOME LENGTH###
######################################################

###NUMBER OF PEAKS PER CHROMOSOME
my.data <- read.delim("number_peaks_per_chromosome.txt", header=FALSE)
head(my.data)
class (my.data$V1)
str(my.data$V1)
class (my.data$V2)
str(my.data$V2)
my.data$V1<-factor (my.data$V1, levels = my.data$V1) #I am telling the levels of the factor being ordered in an specific order as they are in the col V1 of my data
my.data$V1 #level orders of the column have now changed


#Bar graph of chromosomes and number of peaks per chromosome
#fill= to tell the criteria to colour the bars, in this case different colours for different chromosomes
#geom_bar (stat="identity") it means that the statistics showing is only the value of the data, if ="bin" it will calculate the frequencies of each type of data
#aes () for telling which data in each axis
#xlab and ylab for labels names
#theme (legend.title) to specify about legend title
#text.x = element_text(angle = 90, hjust = 1) to write the x axis labels in vertical
plot_peaks_per_chromosome <- ggplot(data=my.data, aes(x=my.data$V1,y=my.data$V2, fill=my.data$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of peaks") + ggtitle("Peaks per chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

plot_peaks_per_chromosome

pdf("Peaks per chromosome.pdf", onefile = TRUE)
plot_peaks_per_chromosome
dev.off() #close the pdf


###NUMBER OF READS PER ENTIRE CHROMOSOME
my.data2 <- read.delim("reads_per_entire_chromosome.txt", header=FALSE)
head(my.data2)
class (my.data2$V1)
str(my.data2$V1)
class (my.data2$V2)
str(my.data2$V2)
my.data2$V1<-factor (my.data2$V1, levels = my.data2$V1) #I am telling the levels of the factor being ordered in an specific order as they are in the col V1 of my data
my.data2$V1 #level orders of the column have now changed


#Bar graph of chromosomes and number of peaks per chromosome
#fill= to tell the criteria to colour the bars, in this case different colours for different chromosomes
#geom_bar (stat="identity") it means that the statistics showing is only the value of the data, if ="bin" it will calculate the frequencies of each type of data
#aes () for telling which data in each axis
#xlab and ylab for labels names
#theme (legend.title) to specify about legend title
#text.x = element_text(angle = 90, hjust = 1) to write the x axis labels in vertical
plot_reads_per_entire_chromosome <- ggplot(data=my.data2, aes(x=my.data2$V1,y=my.data2$V2, fill=my.data2$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of reads") + ggtitle("Reads per entire chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

plot_reads_per_entire_chromosome

pdf("Reads mapping to each entire chromosome.pdf", onefile = TRUE)
plot_reads_per_entire_chromosome
dev.off() #close the pdf



###NUMBER OF READS MAPPING TO ALL PEAKS IN EACH CHROMOSOME
my.data3 <- read.delim("sum_number_of_reads_mapped_to_all_peaks_per_chromosome.txt", header=FALSE)
head(my.data3)
dim(my.data3)
class (my.data2$V1)
str(my.data3$V1)
class (my.data3$V2)
str(my.data3$V2)
my.data3$V1<-factor (my.data3$V1, levels = my.data3$V1) #Function factor () to tell the levels of the factor being ordered in an specific order as they are in the col V1 of my data
my.data3$V1 #level orders of the column have now changed


#Bar graph of chromosomes and number of peaks per chromosome
#fill= to tell the criteria to colour the bars, in this case different colours for different chromosomes
#geom_bar (stat="identity") it means that the statistics showing is only the value of the data, if ="bin" it will calculate the frequencies of each type of data
#aes () for telling which data in each axis
#xlab and ylab for labels names
#theme (legend.title) to specify about legend title
#text.x = element_text(angle = 90, hjust = 1) to write the x axis labels in vertical
plot_number_reads_mapping_to_peaks_per_chromosome <- ggplot(data=my.data3, aes(x=my.data3$V1,y=my.data3$V2, fill=my.data3$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of reads") + ggtitle("Reads mapping to peaks in each chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

plot_number_reads_mapping_to_peaks_per_chromosome

pdf("Reads mapping to peaks in each chromosome.pdf", onefile = TRUE)
plot_number_reads_mapping_to_peaks_per_chromosome
dev.off() #close the pdf



