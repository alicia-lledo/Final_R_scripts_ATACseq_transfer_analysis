getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/QC_files")
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

sum (reads.mapping.peaks.per.chr$V2)



###NUMBER OF READS MAPPING TO NO PEAKS (OFF-TARGET) PER CHROMOSOME
reads.mapping.off.peaks.per.chr <- read.delim("sum_number_of_correctreads_off-peak_per_chromosome.txt", header=FALSE)
reads.mapping.off.peaks.per.chr$V1<-factor (reads.mapping.off.peaks.per.chr$V1, levels = reads.mapping.off.peaks.per.chr$V1) #Function factor () to tell the levels of the factor being ordered in an specific order as they are in the col V1 of my data
reads.mapping.off.peaks.per.chr$V1 #level orders of the column have now changed

sum (reads.mapping.off.peaks.per.chr$V2)



###CHROMOSOME size raw and ungapped

chr.size <- read.delim ("hg19.chrom.sizes.UCSC.txt", header=FALSE)
levels(chr.size$V1)[23] <-"chr23" #to change the name of a particular level of a factor
chr.size$V2<- chr.size$V2/1000000 #chr size in Mb

chr.size.ungapped <- read.delim ("ungapped_chromosome_size.txt", header=FALSE)
levels(chr.size.ungapped$V1)[23] <-"chr23" #to change the name of a particular level of a factor
chr.size.ungapped$V2<- chr.size.ungapped$V2/1000000 #chr size in Mb


###READ COUNT FOR EACH PEAK (WHOLE GENOME)

peak_name_and_read_count <- read.delim ("peak_name_and_read_count.txt", header=FALSE)
nature_original_peak_name_and_read_count <- read.delim ("nature_original_peak_name_and_read_count.txt", header=FALSE)



####NUMBER OF READS MAPPING OFF PEAKS (same length) IN EACH CHROMOSOME
reads.mapping.offpeaks.same.length.per.chr <- read.delim("sum_number_of_reads_mapping_equal_offpeak_regions_bwa.txt", header=FALSE)
reads.mapping.offpeaks.same.length.per.chr$V1<-factor (reads.mapping.offpeaks.same.length.per.chr$V1, levels = reads.mapping.offpeaks.same.length.per.chr$V1) #Function factor () to tell the levels of the factor being ordered in an specific order as they are in the col V1 of my data
reads.mapping.offpeaks.same.length.per.chr$V1 #level orders of the column have now changed



###PLOTS

###Q1.NORMALISATION BASED ON NUMBER OF PEAKS FROM EACH CHROMOSOME PER MAPPED READS TO THAT CHROMOSOME###
#####################################################################################################


table.num_peaks_per_mapped_reads_per_chrom<- peak.per.chr
table.num_peaks_per_mapped_reads_per_chrom$V2<- (peak.per.chr$V2/reads.per.chr$V2)*1000
class(table.num_peaks_per_mapped_reads_per_chrom$V2)
table.num_peaks_per_mapped_reads_per_chrom[24,2]=0 #replace missing value in chrY
plot<- ggplot(data=table.num_peaks_per_mapped_reads_per_chrom, aes(x=table.num_peaks_per_mapped_reads_per_chrom$V1,y=table.num_peaks_per_mapped_reads_per_chrom$V2, fill=table.num_peaks_per_mapped_reads_per_chrom$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of peaks per number reads (x1,000)") + ggtitle("Number of peaks per number of PE reads in each chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
plot
pdf("Q1_Number of peaks per nunmber of PE reads in each chromosome.pdf", onefile = TRUE)
plot
dev.off() #close the pdf





###Q2.NUMBER OF READS MAPPING TO PEAKS PER NUMBER OF TOTAL READS MAPPING TO THAT CHROMOSOME###
#############################################################################################################


table.num_reads_in_peaks_per_total_reads_in_chr<- reads.mapping.peaks.per.chr
table.num_reads_in_peaks_per_total_reads_in_chr$V2<- (reads.mapping.peaks.per.chr$V2/reads.per.chr$V2)*100
class(table.num_reads_in_peaks_per_total_reads_in_chr$V2)
table.num_reads_in_peaks_per_total_reads_in_chr[24,2]=0
plot<- ggplot(data=table.num_reads_in_peaks_per_total_reads_in_chr, aes(x=table.num_reads_in_peaks_per_total_reads_in_chr$V1,y=table.num_reads_in_peaks_per_total_reads_in_chr$V2, fill=table.num_reads_in_peaks_per_total_reads_in_chr$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of PE reads in peaks per total number reads (%)") + ggtitle("Number of PE reads in peaks per number of reads in chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
plot
pdf("Q2_Number of reads in peaks per number of PE reads in Chromosome.pdf", onefile = TRUE)
plot
dev.off() #close the pdf

#A variation of that graph can be done to resemble Fig6 of the ChIPQC report showing the total number of reads as 100% and different colours for the on target and off-target reads

modified.table.num_reads_in_peaks_per_total_reads_in_chr <- table.num_reads_in_peaks_per_total_reads_in_chr
modified.table.num_reads_in_peaks_per_total_reads_in_chr$V3 <- 100-modified.table.num_reads_in_peaks_per_total_reads_in_chr$V2
  #Find the parametres to plot both cols in one bar




###Q3. READS MAPPING TO PEAKS PER NUMBER OF READS MAPPING TO NO PEAKS REGION PER CHROMOSOME (dirty)###
############################################################################################################

table.reads_on_target_vs_reads_off_target_per_chr<- reads.mapping.peaks.per.chr
table.reads_on_target_vs_reads_off_target_per_chr$V2<- (reads.mapping.peaks.per.chr$V2/reads.mapping.off.peaks.per.chr$V2)*100
class(table.reads_on_target_vs_reads_off_target_per_chr$V2)
plot<- ggplot(data=table.reads_on_target_vs_reads_off_target_per_chr, aes(x=table.reads_on_target_vs_reads_off_target_per_chr$V1,y=table.reads_on_target_vs_reads_off_target_per_chr$V2, fill=table.reads_on_target_vs_reads_off_target_per_chr$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of PE reads on peaks per number of reads off-peak (%)") + ggtitle("Number of PE reads in peaks versus off-peak per chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
plot
pdf("Q3_Number of PE reads in peaks versus off-peak per chromosome.pdf", onefile = TRUE)
plot
dev.off() #close the pdf



###Q4. NUMBER OF READS IN PEAKS PER NUMBER OF READS OFF-PEAKS ( non target regions of the same length that peaks) PER CHROMOSOME###
##############################################################################################################################

table.reads.in.peaks.per.reads.off.peak.same.length<- reads.mapping.peaks.per.chr
table.reads.in.peaks.per.reads.off.peak.same.length$V2<- (reads.mapping.peaks.per.chr$V2/reads.mapping.offpeaks.same.length.per.chr$V2)
class(table.reads.in.peaks.per.reads.off.peak.same.length$V2)
table.reads.in.peaks.per.reads.off.peak.same.length[24,2]=0 #replace missing value in chrY
plot<- ggplot(data=table.reads.in.peaks.per.reads.off.peak.same.length, aes(x=table.reads.in.peaks.per.reads.off.peak.same.length$V1,
                                                                            y=table.reads.in.peaks.per.reads.off.peak.same.length$V2, 
                                                                            fill=table.reads.in.peaks.per.reads.off.peak.same.length$V1)) + 
  geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("PE reads fold-enrichment in peak versus off-peak") + ggtitle("Fold enrichment of PE in peaks versus off-peak regions bwa")+ 
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))
plot









##################################################################
########################CORRELATIONS##############################
#################################################################


###1.CORRELATION BETWEEN CHROMOSOME SIZE AND Q1 (PEAKS VS READS MAPPED PER CHROMOSOME)

#Create a data fram with 3 columns: 
  #col1 chr name
  #col2: number of peaks per number of reads from Q1
  #col3: chromosome size
correlation.Q1.chr.size <- data.frame(chromosome_size=chr.size$V2[1:23], Q1=table.num_peaks_per_mapped_reads_per_chrom$V2[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation.Q1.chr.size2 <- data.frame(chromosome=chr.size$V1[1:23], chromosome_size=chr.size$V2[1:23], Q1=table.num_peaks_per_mapped_reads_per_chrom$V2[1:23])
correlation.Q1.chr.size2$chromosome <-factor (correlation.Q1.chr.size2$chromosome, levels=correlation.Q1.chr.size2$chromosome)
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.Q1.chr.size, use="everything", method="pearson")
ggplot(correlation.Q1.chr.size2, aes(x=correlation.Q1.chr.size2$chromosome_size, y=correlation.Q1.chr.size2$Q1, color=correlation.Q1.chr.size2$chromosome)) +
  xlab("Chromosome size (Mb)") + 
  ylab("Q1.Number of peaks per number of total PE reads (%)") + ggtitle("Correlation between chromosome size and Q1") +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))    # Don't add shaded confidence region

#Repeat for the ungapped chr size

correlation.Q1.chr.size <- data.frame(chromosome_size=chr.size.ungapped$V2[1:23], Q1=table.num_peaks_per_mapped_reads_per_chrom$V2[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation.Q1.chr.size2 <- data.frame(chromosome=chr.size.ungapped$V1[1:23], chromosome_size=chr.size.ungapped$V2[1:23], Q1=table.num_peaks_per_mapped_reads_per_chrom$V2[1:23])
correlation.Q1.chr.size2$chromosome <-factor (correlation.Q1.chr.size2$chromosome, levels=correlation.Q1.chr.size2$chromosome)
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.Q1.chr.size, use="everything", method="pearson")
ggplot(correlation.Q1.chr.size2, aes(x=correlation.Q1.chr.size2$chromosome_size, y=correlation.Q1.chr.size2$Q1, color=correlation.Q1.chr.size2$chromosome)) +
  xlab("Chromosome size (Mb)") + 
  ylab("Q1.Number of peaks per number of total PE reads (%)") + ggtitle("Correlation between chromosome size ungapped and Q1") +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))    # Don't add shaded confidence region



###2.CORRELATION BETWEEN CHROMOSOME SIZE AND NUMBER OF PEAKS/READS

######PEAKS
#Create a data fram with 3 columns: 
#col1 chr name
#col2: number of peaks 
#col3: chromosome size
correlation.chr.size.peaks <- data.frame(chromosome_size=chr.size$V2[1:23], num.peaks=peak.per.chr$V2[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation.chr.size.peaks2 <- data.frame(chromosome=chr.size$V1[1:23],chromosome_size=chr.size$V2[1:23], num.peaks=peak.per.chr$V2[1:23])
correlation.chr.size.peaks2$chromosome <-factor (correlation.chr.size.peaks2$chromosome, levels=correlation.chr.size.peaks2$chromosome) #to order the chr descendent
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.chr.size.peaks, use="everything", method="pearson")
ggplot(correlation.chr.size.peaks2, aes(x=correlation.chr.size.peaks2$chromosome_size, y=correlation.chr.size.peaks2$num.peaks, color=correlation.chr.size.peaks2$chromosome)) +
  xlab("Chromosome size (Mb)") + 
  ylab("Number of peaks") + ggtitle("Correlation between chromosome size and number of peaks") +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))    # Don't add shaded confidence region


######READS
#Create a data fram with 3 columns: 
#col1 chr name
#col2: number of reads
#col3: chromosome size
correlation.chr.size.reads <- data.frame(chromosome_size=chr.size$V2[1:23], num.reads=reads.per.chr$V2[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation.chr.size.reads2 <- data.frame(chromosome=chr.size$V1[1:23],chromosome_size=chr.size$V2[1:23], num.reads=reads.per.chr$V2[1:23])
correlation.chr.size.reads2$chromosome <-factor (correlation.chr.size.reads2$chromosome, levels=correlation.chr.size.reads2$chromosome) #to order the chr descendent
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.chr.size.reads, use="everything", method="pearson")
ggplot(correlation.chr.size.reads2, aes(x=correlation.chr.size.reads2$chromosome_size, y=correlation.chr.size.reads2$num.reads, color=correlation.chr.size.reads2$chromosome)) +
  xlab("Chromosome size (Mb)") + 
  ylab("Number of PE reads") + ggtitle("Correlation between chromosome size and number of PE reads") +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))    # Don't add shaded confidence region



###3.CORRELATION BETWEEN Q1 AND Q2

#Create a data fram with 3 columns: 
#col1 chr name
#col2:Q1
#col3: Q3
correlation.Q1.and.Q2 <- data.frame(Q1=table.num_peaks_per_mapped_reads_per_chrom$V2[1:23], Q2=table.num_reads_in_peaks_per_total_reads_in_chr$V2[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation.Q1.and.Q2_2 <- data.frame(chromosome=chr.size$V1[1:23],Q1=table.num_peaks_per_mapped_reads_per_chrom$V2[1:23], Q2=table.num_reads_in_peaks_per_total_reads_in_chr$V2[1:23])
correlation.Q1.and.Q2_2$chromosome <-factor (correlation.Q1.and.Q2_2$chromosome, levels=correlation.Q1.and.Q2_2$chromosome)
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.Q1.and.Q2, use="everything", method="pearson")
ggplot(correlation.Q1.and.Q2_2, aes(x=correlation.Q1.and.Q2_2$Q1, y=correlation.Q1.and.Q2_2$Q2, color=correlation.Q1.and.Q2_2$chromosome)) +
  xlab("Q1") + 
  ylab("Q2") + ggtitle("Correlation between Q1 and Q2") +
  geom_point(shape=1) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))    # Don't add shaded confidence region


