getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/bowtie2_analysis")
getwd ()

# Use ggplot library to draw graphs
library(ggplot2)

####READING TABLES###

###NUMBER OF PEAKS PER CHROMOSOME
peak.per.chr <- read.delim("peaks_per_chromosome_bowtie2.txt", header=FALSE)
peak.per.chr$V1<-factor (peak.per.chr$V1, levels = peak.per.chr$V1) #I am telling the levels of the factor being ordered in an specific order as they are in the col V1 of my data
peak.per.chr$V1 #level orders of the column have now changed

###CHROMOSOME SIZE UNGAPPED
chr.size.ungapped <- read.delim ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/QC_files/ungapped_chromosome_size.txt", header=FALSE)
#levels(chr.size.ungapped$V1)[23] <-"chr23" #to change the name of a particular level of a factor
chr.size.ungapped$V2<- chr.size.ungapped$V2/1000000 #chr size in Mb


###NUMBER OF READS PER ENTIRE CHROMOSOME
reads.per.chr <- read.delim("sum_correctreads_per_entire_chromosome.txt", header=FALSE) #it doesn't change based on the peak calling
reads.per.chr$V1<-factor (reads.per.chr$V1, levels = reads.per.chr$V1) #I am telling the levels of the factor being ordered in an specific order as they are in the col V1 of my data
reads.per.chr$V1 #level orders of the column have now changed


####NUMBER OF READS MAPPING TO ALL PEAKS IN EACH CHROMOSOME
reads.mapping.peaks.per.chr <- read.delim("sum_number_of_reads_mapped_to_peaks_per_chromosome_bowtie2.txt", header=FALSE)
reads.mapping.peaks.per.chr$V1<-factor (reads.mapping.peaks.per.chr$V1, levels = reads.mapping.peaks.per.chr$V1) #Function factor () to tell the levels of the factor being ordered in an specific order as they are in the col V1 of my data
reads.mapping.peaks.per.chr$V1 #level orders of the column have now changed


####NUMBER OF READS MAPPING OFF PEAKS (same length) IN EACH CHROMOSOME
reads.mapping.offpeaks.same.length.per.chr <- read.delim("sum_number_of_reads_mapping_equal_offpeak_regions_bowtie2.txt", header=FALSE)
reads.mapping.offpeaks.same.length.per.chr$V1<-factor (reads.mapping.offpeaks.same.length.per.chr$V1, levels = reads.mapping.offpeaks.same.length.per.chr$V1) #Function factor () to tell the levels of the factor being ordered in an specific order as they are in the col V1 of my data
reads.mapping.offpeaks.same.length.per.chr$V1 #level orders of the column have now changed



###PLOTS####

###PLOT OF DISTRIBUTION OF NUMBER OF PEAKS PER CHROMOSOME

plot1<- ggplot(data=peak.per.chr, aes(x=peak.per.chr$V1,y=peak.per.chr$V2, fill=peak.per.chr$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of peaks") + ggtitle("Number of peaks per chromosome")+ 
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))
plot1

pdf("Number of peaks per chromosome (bowtie2).pdf", onefile = TRUE)
plot1
dev.off() #close the pdf
  
#CORRELATION OF NUMBER OF PEAKS PER CHROMOSOME AND UNGAPPED CHROMOSOME SIZE

correlation.peaksperchr.and.chr.size <- data.frame(chromosome_size=chr.size.ungapped$V2, peaks.per.chr=peak.per.chr$V2)

#Create the same data frame adding the chromosome name for ggplot2
correlation.peaksperchr.and.chr.size2 <- data.frame(chromosome=chr.size.ungapped$V1, chromosome_size=chr.size.ungapped$V2, peaks.per.chr=peak.per.chr$V2)
correlation.peaksperchr.and.chr.size2$chromosome <-factor (correlation.peaksperchr.and.chr.size2$chromosome, levels=correlation.peaksperchr.and.chr.size2$chromosome)
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.peaksperchr.and.chr.size, use="everything", method="pearson")
cor.test(correlation.peaksperchr.and.chr.size$chromosome_size, correlation.peaksperchr.and.chr.size$peaks.per.chr)
plot2<-ggplot(correlation.peaksperchr.and.chr.size2, aes(x=correlation.peaksperchr.and.chr.size2$chromosome_size, y=correlation.peaksperchr.and.chr.size2$peaks.per.chr, color=correlation.peaksperchr.and.chr.size2$chromosome)) +
  xlab("Chromosome size (Mb)") + 
  ylab("Number of peaks per chromosome") + ggtitle("Correlation between chromosome size\n and number of peaks per chr") +
  geom_point(size=3) +    # Use hollow circles geom_point(shape=1) empty for filled
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))+ # Don't add shaded confidence region
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))
plot2

pdf("Correlation between chromosome size and number of peaks per chr (bowtie2).pdf", onefile = TRUE)
plot2
dev.off() #close the pdf


###Q1.NORMALISATION BASED ON NUMBER OF PEAKS FROM EACH CHROMOSOME PER MAPPED READS TO THAT CHROMOSOME###
#####################################################################################################


table.num_peaks_per_mapped_reads_per_chrom<- peak.per.chr
table.num_peaks_per_mapped_reads_per_chrom$V2<- (peak.per.chr$V2/reads.per.chr$V2)*100
class(table.num_peaks_per_mapped_reads_per_chrom$V2)
table.num_peaks_per_mapped_reads_per_chrom[24,2]=0 #replace missing value in chrY
plot3<- ggplot(data=table.num_peaks_per_mapped_reads_per_chrom, aes(x=table.num_peaks_per_mapped_reads_per_chrom$V1,y=table.num_peaks_per_mapped_reads_per_chrom$V2, fill=table.num_peaks_per_mapped_reads_per_chrom$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of peaks per number of PE reads (%)") + ggtitle("Number of peaks per number of\n PE reads (%) in each chromosome (bowtie2)")+ 
  theme(plot.title = element_text(size=16,face="bold"), legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.title.x = element_text(face="bold",size=15),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=14))
plot3

pdf("Number of peaks per number of PE reads in each chromosome (Q1 bowtie2).pdf", onefile = TRUE)
plot3
dev.off() #close the pdf


  #CORRELATION CHROMOSOME SIZE AND Q1
#for the ungapped chr size

correlation.Q1.chr.size <- data.frame(chromosome_size=chr.size.ungapped$V2[1:23], Q1=table.num_peaks_per_mapped_reads_per_chrom$V2[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation.Q1.chr.size2 <- data.frame(chromosome=chr.size.ungapped$V1[1:23], chromosome_size=chr.size.ungapped$V2[1:23], Q1=table.num_peaks_per_mapped_reads_per_chrom$V2[1:23])
correlation.Q1.chr.size2$chromosome <-factor (correlation.Q1.chr.size2$chromosome, levels=correlation.Q1.chr.size2$chromosome)
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.Q1.chr.size, use="everything", method="pearson")
cor.test(correlation.Q1.chr.size$chromosome_size, correlation.Q1.chr.size$Q1)
plot4<-ggplot(correlation.Q1.chr.size2, aes(x=correlation.Q1.chr.size2$chromosome_size, y=correlation.Q1.chr.size2$Q1, color=correlation.Q1.chr.size2$chromosome)) +
  xlab("Chromosome size (Mb)") + 
  ylab("Peaks per PE reads (%)") + ggtitle("Correlation between chromosome size and Q1") +
  geom_point() +    # Use filled circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))+    # Don't add shaded confidence region
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), 
        axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18)) #no legend title


plot4

pdf("Correlation between chromosome size and Q1 (bowtie2).pdf", onefile = TRUE)
plot4
dev.off() #close the pdf



###Q2.NUMBER OF READS MAPPING TO PEAKS PER NUMBER OF TOTAL READS MAPPING TO THAT CHROMOSOME###
#############################################################################################################


table.num_reads_in_peaks_per_total_reads_in_chr<- reads.mapping.peaks.per.chr
table.num_reads_in_peaks_per_total_reads_in_chr$V2<- (reads.mapping.peaks.per.chr$V2/reads.per.chr$V2)*100
class(table.num_reads_in_peaks_per_total_reads_in_chr$V2)
table.num_reads_in_peaks_per_total_reads_in_chr[24,2]=0
plot5<- ggplot(data=table.num_reads_in_peaks_per_total_reads_in_chr, aes(x=table.num_reads_in_peaks_per_total_reads_in_chr$V1,y=table.num_reads_in_peaks_per_total_reads_in_chr$V2, fill=table.num_reads_in_peaks_per_total_reads_in_chr$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of PE reads in peaks\n per total number PE reads (%)") + 
  ggtitle("Number of PE reads in peaks per number of\n PE reads in each chromosome")+ 
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))
plot5

pdf("Number of PE reads in peaks per number of PE reads in each chromosome (Q2 bowtie2).pdf", onefile = TRUE)
plot5
dev.off() #close the pdf



#A variation of that graph can be done to resemble Fig6 of the ChIPQC report showing the total number of reads as 100% and different colours for the on target and off-target reads

modified.table.num_reads_in_peaks_per_total_reads_in_chr <- table.num_reads_in_peaks_per_total_reads_in_chr
modified.table.num_reads_in_peaks_per_total_reads_in_chr$V3 <- 100-modified.table.num_reads_in_peaks_per_total_reads_in_chr$V2
#Find the parametres to plot both cols in one bar





###Q4. NUMBER OF READS IN PEAKS PER NUMBER OF READS OFF-PEAKS ( non target regions of the same length that peaks) PER CHROMOSOME###
##############################################################################################################################

table.reads.in.peaks.per.reads.off.peak.same.length<- reads.mapping.peaks.per.chr
table.reads.in.peaks.per.reads.off.peak.same.length$V2<- (reads.mapping.peaks.per.chr$V2/reads.mapping.offpeaks.same.length.per.chr$V2)
class(table.reads.in.peaks.per.reads.off.peak.same.length$V2)
table.reads.in.peaks.per.reads.off.peak.same.length[24,2]=0 #replace missing value in chrY
plot6<- ggplot(data=table.reads.in.peaks.per.reads.off.peak.same.length, aes(x=table.reads.in.peaks.per.reads.off.peak.same.length$V1,
                                                                            y=table.reads.in.peaks.per.reads.off.peak.same.length$V2, 
                                                                            fill=table.reads.in.peaks.per.reads.off.peak.same.length$V1)) + 
geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("PE reads fold-enrichment in peak versus off-peak") + ggtitle("Fold enrichment of PE in peaks\n versus off-peak regions bowtie2")+ 
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))
plot6

pdf("Fold enrichment of PE in peaks versus off-peak regions (Q4 bowtie2).pdf", onefile = TRUE)
plot6
dev.off() #close the pdf







###CORRELATION BETWEEN CHROMOSOME SIZE AND NUMBER OF PEAKS/READS

######PEAKS
#Create a data fram with 3 columns: 
#col1 chr name
#col2: number of peaks 
#col3: chromosome size
correlation.chr.size.peaks <- data.frame(chromosome_size=chr.size.ungapped$V2 [1:23], num.peaks=peak.per.chr$V2[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation.chr.size.peaks2 <- data.frame(chromosome=chr.size.ungapped$V1[1:23],chromosome_size=chr.size.ungapped$V2[1:23], num.peaks=peak.per.chr$V2[1:23])
correlation.chr.size.peaks2$chromosome <-factor (correlation.chr.size.peaks2$chromosome, levels=correlation.chr.size.peaks2$chromosome) #to order the chr descendent
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.chr.size.peaks, use="everything", method="pearson")
cor.test (correlation.chr.size.peaks$chromosome_size, correlation.chr.size.peaks$num.peaks)
plot7<-ggplot(correlation.chr.size.peaks2, aes(x=correlation.chr.size.peaks2$chromosome_size, y=correlation.chr.size.peaks2$num.peaks, color=correlation.chr.size.peaks2$chromosome)) +
  xlab("Chromosome size (Mb)") + 
  ylab("Number of peaks") + ggtitle("Correlation between chromosome size\n and number of peaks") +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))+
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), 
        axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))# Don't add shaded confidence region

plot7

pdf("Correlation between chromosome size and number of peaks (bowtie2).pdf", onefile = TRUE)
plot7
dev.off() #close the pdf


######READS
#Create a data fram with 3 columns: 
#col1 chr name
#col2: number of reads
#col3: chromosome size
correlation.chr.size.reads <- data.frame(chromosome_size=chr.size.ungapped$V2[1:23], num.reads=reads.per.chr$V2[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation.chr.size.reads2 <- data.frame(chromosome=chr.size.ungapped$V1[1:23],chromosome_size=chr.size.ungapped$V2[1:23], num.reads=reads.per.chr$V2[1:23])
correlation.chr.size.reads2$chromosome <-factor (correlation.chr.size.reads2$chromosome, levels=correlation.chr.size.reads2$chromosome) #to order the chr descendent
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.chr.size.reads, use="everything", method="pearson")
cor.test(correlation.chr.size.reads$chromosome_size, correlation.chr.size.reads$num.reads)
plot8<-ggplot(correlation.chr.size.reads2, aes(x=correlation.chr.size.reads2$chromosome_size, y=correlation.chr.size.reads2$num.reads, color=correlation.chr.size.reads2$chromosome)) +
  xlab("Chromosome size (Mb)") + 
  ylab("Number of PE reads") + ggtitle("Correlation between chromosome size and\n number of PE reads") +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))+
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), 
        axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))# Don't add shaded confidence region


plot8

pdf("Correlation between chromosome size and number of reads (bowtie2).pdf", onefile = TRUE)
plot8
dev.off() #close the pdf



###CORRELATION BETWEEN Q1 AND Q2: 

#Create a data fram with 3 columns: 
#col1 chr name
#col2:Q1
#col3: Q3
correlation.Q1.and.Q2 <- data.frame(Q1=table.num_peaks_per_mapped_reads_per_chrom$V2[1:23], Q2=table.num_reads_in_peaks_per_total_reads_in_chr$V2[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation.Q1.and.Q2_2 <- data.frame(chromosome=chr.size.ungapped$V1[1:23],Q1=table.num_peaks_per_mapped_reads_per_chrom$V2[1:23], Q2=table.num_reads_in_peaks_per_total_reads_in_chr$V2[1:23])
correlation.Q1.and.Q2_2$chromosome <-factor (correlation.Q1.and.Q2_2$chromosome, levels=correlation.Q1.and.Q2_2$chromosome)
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.Q1.and.Q2, use="everything", method="pearson")
cor.test(correlation.Q1.and.Q2$Q1, correlation.Q1.and.Q2$Q2)
plot9<-ggplot(correlation.Q1.and.Q2_2, aes(x=correlation.Q1.and.Q2_2$Q1, y=correlation.Q1.and.Q2_2$Q2, color=correlation.Q1.and.Q2_2$chromosome)) +
  xlab("Q1") + 
  ylab("Q2") + ggtitle("Correlation between Q1 and Q2") +
  geom_point() + 
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))+  # Don't add shaded confidence region
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), 
        axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))# Don't add shaded confidence region


plot9

pdf("Correlation between Q1 and Q2 (bowtie2).pdf", onefile = TRUE)
plot9
dev.off() #close the pdf



###CORRELATION BETWEEN CHROMOSOME SIZE AND NUMBER OF GENES PER CHROMOSOME: 

#First:
#To pull only genes in chr1 to chr22,X and Y from the file /well/jknight/ATACseq/ATACseq_001/Analysis/QC_plots/refseq_genes_hg19.txt with refseq gene list
#To coundt the number of genes (overlapping) contained in each chromosome

gene.file <- read.delim("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/refseq_genes_hg19.txt") # 54439 16
head(gene.file)

# remove genes in blacklist regions
to.remove <- grep("_.*", gene.file$chrom)
gene.file <- gene.file[-(to.remove), ]
gene.file <- droplevels(gene.file)#droplevels to drop unused levels from a factor or, more commonly, from factors in a data frame.
head(gene.file)
# count the number of genes on each chromosome
gene.count <- data.frame(table(gene.file$chrom)) #table counts frequency of events at a particular factor combination, in this case chr name
head(gene.count)
gene.count$Sort<-gene.count$Var1 #create a new column
  #Sort (or order) a vector or factor (partially) into ascending or descending order
gene.count$Sort<-gsub ("chr", "",gene.count$Sort) #search for a pattern of text and replace with the desired text or nothing
gene.count$Sort<-gsub ("X", "23",gene.count$Sort) #chr X and Y need to be replaced separately as characters with number 23 and 24
gene.count$Sort<-gsub ("Y", "24",gene.count$Sort)
gene.count$Sort<-as.numeric(gene.count$Sort)#convert the sort col to numeric to be able to apply the order function
gene.count<-gene.count[order(gene.count$Sort),]#replace gene.counts by a data frame where the rows will be ordered by the value of the Sort column and all the col will be printed


#Create a data fram with 3 columns: 
#col1 chr name
#col2: number of genes
#col3: chromosome size
correlation.chr.size.number.genes <- data.frame(chromosome_size=chr.size.ungapped$V2[1:23], num.genes=gene.count$Freq[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation.chr.size.number.genes2 <- data.frame(chromosome=chr.size.ungapped$V1[1:23],chromosome_size=chr.size.ungapped$V2[1:23], num.genes=gene.count$Freq[1:23])
correlation.chr.size.number.genes2$chromosome <-factor (correlation.chr.size.number.genes2$chromosome, levels=correlation.chr.size.number.genes2$chromosome) #to order the chr descendent
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.chr.size.number.genes, use="everything", method="pearson")
cor.test(correlation.chr.size.number.genes$chromosome_size, correlation.chr.size.number.genes$num.genes)
plot10<-ggplot(correlation.chr.size.number.genes2, aes(x=correlation.chr.size.number.genes2$chromosome_size, y=correlation.chr.size.number.genes$num.genes, color=correlation.chr.size.number.genes2$chromosome)) +
  xlab("Chromosome size (Mb)") + 
  ylab("Number of genes") + ggtitle("Correlation between chromosome size and\n number of genes") +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))+
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), 
        axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))# Don't add shaded confidence region
plot10

pdf("Correlation between chromosme size and number of genes per chromosome (bowtie2).pdf", onefile = TRUE)
plot10
dev.off() #close the pdf




###CORRELATION BETWEEN NUMBER OF GENES AND PEAKS MAPPED TO CHROMOSOME 

correlation.num.genes.peaks.per.chr<- data.frame(num.peaks=peak.per.chr$V2[1:23], num.genes=gene.count$Freq[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation.num.genes.peaks.per.chr2<- data.frame(chromosome=chr.size.ungapped$V1[1:23], num.peaks=peak.per.chr$V2[1:23], num.genes=gene.count$Freq[1:23])
correlation.num.genes.peaks.per.chr2$chromosome <-factor (correlation.num.genes.peaks.per.chr2$chromosome, levels=correlation.num.genes.peaks.per.chr2$chromosome) #to order the chr descendent
#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation.num.genes.peaks.per.chr, use="everything", method="pearson")
cor.test(correlation.num.genes.peaks.per.chr$num.genes, correlation.num.genes.peaks.per.chr$num.peaks)
plot11<-ggplot(correlation.chr.size.number.genes2, aes(x=correlation.chr.size.number.genes2$num.genes, y=correlation.num.genes.peaks.per.chr2$num.peaks, color=correlation.chr.size.number.genes2$chromosome)) +
  xlab("Number of genes") + 
  ylab("Number of peaks") + ggtitle("Correlation between number of genes and\n number of peaks per chromosome") +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))+
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), 
        axis.title.x = element_text(face="bold",size=18),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))#Don't add shaded confidence region
plot11

pdf("Correlation between number of genes and number of peaks per chromosome (bowtie2).pdf", onefile = TRUE)
plot10
dev.off() #close the pdf


###CORRELATION BETWEEN NUMBER OF PEAKS AND NUMBER OF READS

correlation.reads.peaks <- data.frame (num.reads=reads.per.chr$V2[1:23], num.peaks=peak.per.chr$V2[1:23])
correlation.reads.peaks2 <- data.frame (chromosome=chr.size.ungapped$V1[1:23], num.read=reads.per.chr$V2[1:23], num.peaks=peak.per.chr$V2[1:23])
cor(correlation.reads.peaks, use="everything", method="pearson")
cor.test(correlation.reads.peaks$num.reads,correlation.reads.peaks$num.peaks)


###CORRELATION Q1 (NUM OF PEAKS PER READS) VS GENE DENSITY (NUM GENES PER CHR LENGTH)

gene.density <- data.frame (chromosome=chr.size.ungapped$V1[1:23], chr.length=chr.size.ungapped$V2[1:23], genes=gene.count$Freq[1:23])
gene.density$chromosome <-factor (gene.density$chromosome, levels=gene.density$chromosome)
gene.density$V4 <- gene.density$genes/gene.density$chr.length

plot12<- ggplot(data=gene.density, aes(x=gene.density$chromosome,y=gene.density$V4, fill=gene.density$chromosome)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of genes/chromosome size (Mb)") + ggtitle("Gene density")+ 
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))
plot12

correlation.gene.density.Q1 <- data.frame (Q1=table.num_peaks_per_mapped_reads_per_chrom$V2[1:23], density=gene.density$V4)
cor(correlation.gene.density.Q1, use="everything", method="pearson" )
cor.test (correlation.gene.density.Q1$Q1, correlation.gene.density.Q1$density)



###CORRELATION GENE DENSITY AND CHR SIZE

correlation.gene.density.chr.size <- data.frame (chr.size=gene.density$chr.length[1:23], density=gene.density$V4)
cor(correlation.gene.density.chr.size, use="everything", method="pearson" )
cor.test (correlation.gene.density.chr.size$chr.size, correlation.gene.density.chr.size$density)

