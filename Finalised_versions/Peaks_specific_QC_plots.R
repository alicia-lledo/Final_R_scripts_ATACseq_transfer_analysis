getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/pipeline_output_all_samples_QC/correct_blood_bwa/")
getwd ()

# Use ggplot library to draw graphs
library(ggplot2)

####READING TABLES###

###NUMBER OF PEAKS PER CHROMOSOME
peak.per.chr <- read.delim("my_data_number_of_merged_peaks_per_chromosome.txt", header=FALSE)
peak.per.chr$V1<-factor (peak.per.chr$V1, levels = peak.per.chr$V1) #I am telling the levels of the factor being ordered in an specific order as they are in the col V1 of my data
peak.per.chr$V1 #level orders of the column have now changed



###NUMBER OF READS PER ENTIRE CHROMOSOME
reads.per.chr <- read.delim("correctreads_per_entire_chromosome.txt", header=FALSE) #it doesn't change based on the peak calling
reads.per.chr$V1<-factor (reads.per.chr$V1, levels = reads.per.chr$V1) #I am telling the levels of the factor being ordered in an specific order as they are in the col V1 of my data
reads.per.chr$V1 #level orders of the column have now changed




###NUMBER OF READS MAPPING TO ALL PEAKS IN EACH CHROMOSOME
reads.mapping.peaks.per.chr <- read.delim("sum_number_of_reads_mapped_to_peaks_per_chromosome.txt", header=FALSE)
reads.mapping.peaks.per.chr$V1<-factor (reads.mapping.peaks.per.chr$V1, levels = reads.mapping.peaks.per.chr$V1) #Function factor () to tell the levels of the factor being ordered in an specific order as they are in the col V1 of my data
reads.mapping.peaks.per.chr$V1 #level orders of the column have now changed



###NUMBER OF READS MAPPING TO NO PEAKS (OFF-TARGET) PER CHROMOSOME
reads.mapping.off.peaks.per.chr <- read.delim("sum_number_of_correctreads_off-peak_per_chromosome.txt", header=FALSE)
reads.mapping.off.peaks.per.chr$V1<-factor (reads.mapping.off.peaks.per.chr$V1, levels = reads.mapping.off.peaks.per.chr$V1) #Function factor () to tell the levels of the factor being ordered in an specific order as they are in the col V1 of my data
reads.mapping.off.peaks.per.chr$V1 #level orders of the column have now changed



###CHROMOSOME size raw and ungapped

chr.size <- read.delim ("hg19.chrom.sizes.UCSC.txt", header=FALSE)
levels(chr.size$V1)[23] <-"chr23" #to change the name of a particular level of a factor
chr.size$V2<- chr.size$V2/1000000 #chr size in Mb

chr.size.ungapped <- read.delim ("ungapped_chromosome_size.txt", header=FALSE)
levels(chr.size.ungapped$V1)[23] <-"chr23" #to change the name of a particular level of a factor
chr.size.ungapped$V2<- chr.size.ungapped$V2/1000000 #chr size in Mb


###READ COUNT FOR EACH PEAK (WHOLE GENOME)

peak_name_and_read_count <- read.delim ("peak_name_and_read_count.bed", header=FALSE)
nature_original_peak_name_and_read_count <- read.delim ("nature_peak_name_and_read_count.bed", header=FALSE)



median(peak_name_and_read_count$V2)
median(nature_original_peak_name_and_read_count$V2)
median(nature_fastq_peak_name_and_read_count$V2)

###Q5. DISTRIBUTION OF READS IN ALL PEAKS PER CHROMOSOME###
#############################################################################################################

#Generate a boxplot plot for the counts in the tables showing counts of reads per peak

peak_name_and_read_count$V3 <- "JKnight lab"
nature_original_peak_name_and_read_count$V3 <- "Greenleaf lab"

merged_data <- rbind (peak_name_and_read_count, nature_original_peak_name_and_read_count[which(nature_original_peak_name_and_read_count$V2<1000),])

merged.plot<-ggplot(merged_data, aes(x=merged_data$V3, y=merged_data$V2, fill=merged_data$V3)) + geom_boxplot()+
  theme_bw()+theme (axis.title.x = element_blank())+ ylab("Read counts in peaks")+ theme(plot.title = element_text(size=20,face="bold"),
                                                                                                             axis.text.x = element_text(angle = 0, face="bold", size=18),
      axis.text.y = element_text(face="bold",size=18), axis.title.y = element_text(face="bold",size=18))+scale_fill_manual(values=c("#999999", "#E69F00"), 
                  name="Data source",labels=c("JKnight lab","Greenleaf lab"))
merged.plot
theme(legend.key.size = unit(2.5, "cm"))

median (peak_name_and_read_count$V2)
median (nature_original_peak_name_and_read_count$V2)

#To generate a density plot for the same data
max(peak_name_and_read_count$V2)
max(nature_original_peak_name_and_read_count$V2)
length(which(nature_original_peak_name_and_read_count$V2>300))

#density plot
ggplot(merged_mydata_nature_original, aes(x=merged_mydata_nature_original$V2, colour=merged_mydata_nature_original$V3))+
  geom_density()+xlim(0, 300)+ xlab("Number of PE reads") + ylab("Density")+theme(legend.title=element_blank())+
  ggtitle ("Density of number of reads per peak")

#frequency plot
ggplot(merged_mydata_nature_original, aes(x=merged_mydata_nature_original$V2, fill=merged_mydata_nature_original$V3)) +
  geom_histogram(binwidth=5, alpha=.5, position="identity")+xlim (0, 300)+xlab("Number of PE reads") + ylab("Frequency")+
  ggtitle ("Frequency of number of reads per peak")+theme(legend.title=element_blank())