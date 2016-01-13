getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots")
getwd ()

# Use ggplot library to draw graphs
library(ggplot2)


###NORMALISATION BASED ON GENE NUMBER PER CHROMOSOME###
######################################################

#First part of the script:
#First:
#To pull only genes in chr1 to chr22,X and Y from the file /well/jknight/ATACseq/ATACseq_001/Analysis/QC_plots/refseq_genes_hg19.txt
#To coundt the number of genes (overlapping) contained in each chromosome

gene.file <- read.delim("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/refseq_genes_hg19.txt") # 54439 16
head(gene.file)

# remove genes in blacklist regions
to.remove <- grep("_.*", gene.file$chrom)
gene.file <- gene.file[-(to.remove), ]
gene.file <- droplevels(gene.file)#droplevels to drop unused levels from a factor or, more commonly, from factors in a data frame.
head(gene.file)
# count the number of genes on each chromosome
gene.count <- data.frame(table(gene.file$chrom))
head(gene.count)
gene.count$Sort<-gene.count$Var1 #create a new column
gene.count$Sort<-gsub ("chr", "",gene.count$Sort) #search for a pattern of text and replace with the desired text or nothing
gene.count$Sort<-gsub ("X", "23",gene.count$Sort) #chr X and Y need to be replaced separately as characters with number 23 and 24
gene.count$Sort<-gsub ("Y", "24",gene.count$Sort)
gene.count$Sort<-as.numeric(gene.count$Sort)#convert the sort col to numerci to be able to apply the order function
gene.count<-gene.count[order(gene.count$Sort),]#replace gene.counts by a data frame where the rows will be ordered by the value of the Sort column and all the col will be printed

#Second part:
# This converts your table of genes into a GRanges object - so each peak is 
# represented as a range from its start to stop position.  
library(GenomicRanges)

gene.file$Width <- gene.file$txEnd - gene.file$txStart
genes <- with(gene.file, GRanges(seqnames=Rle(paste0(chrom)), 
                                 ranges=IRanges(start=txStart, width=Width)))

# Some of these may be overlapping e.g. on chromosome 1:
# coverage counts the number of ranges that cover each position - so counts how 
# many genes include each base
base.coverage <- coverage(genes) #it counst how many ranges can be mapped to each bp in all the chromosomes. Ideally if ranges weren't overlapping each bp should only be contained in one range or in none
hist(base.coverage@listData$chr1@values)

# Reduce first orders the ranges from x to right, then merges the overlapping or 
# adjacent ones.  Therefore each base that in within a gene is only counted once
genes <- reduce(genes) #similar to intersect merge tool to merge all the intervals overlapping to make a single one

# To check that no position is counted more than once:
base.coverage <- coverage(genes)
hist(base.coverage@listData$chr1@values)

# To calculate the total number of bases that are in genes for each chromosome
Chr <- rep(genes@seqnames@values, as.numeric(genes@seqnames@lengths)) #a loop to pull the chr names
gene.table <- data.frame(cbind(Chr, genes@ranges@width))#loop to retrieve all the gene widths for each chrom and add all the bp 
colnames(gene.table) <- c("Chr", "Width")
gene_bp<-aggregate(gene.table$Width, by=list(Category=gene.table$Chr), FUN=sum) #for each crhomosome it will add up all the widths and assign them to the category which is the chr number


#Calculate gene density per chromosome
#Read the h19 chromosome size file and sort it by chromosome number
hg19.chromosome.size<- read.delim ("hg19.chrom.sizes.UCSC.txt", header=FALSE)
hg19.chromosome.size$Sort <- hg19.chromosome.size$V1#factor levels not following chromosome ordered numbers
hg19.chromosome.size$Sort<-gsub ("chr", "",hg19.chromosome.size$Sort) #search for a pattern of text and replace with the desired text or nothing
hg19.chromosome.size$Sort<-gsub ("X", "23",hg19.chromosome.size$Sort) #chr X and Y need to be replaced separately as characters with number 23 and 24
hg19.chromosome.size$Sort<-gsub ("Y", "24",hg19.chromosome.size$Sort)
hg19.chromosome.size$Sort<-as.numeric(hg19.chromosome.size$Sort)
class(hg19.chromosome.size$Sort)
hg19.chromosome.size<- hg19.chromosome.size[order(hg19.chromosome.size$Sort),]
gene.density.per.chr <- gene.count$Freq/(hg19.chromosome.size$V2/1000000) #gene density per Mb
chr.label<-c(1:24)
chr.label<- as.character(chr.label)
table.gene.density<-data.frame(cbind(gene.density.per.chr, chr.label))
table.gene.density$gene.density.per.chr <- as.numeric(table.gene.density$gene.density.per.chr)
table.gene.density$chr.label<-factor (table.gene.density$chr.label, levels = table.gene.density$chr.label)
plot.gene.density <- ggplot(data=table.gene.density, aes(x=table.gene.density$chr.label,y=table.gene.density$gene.density.per.chr, fill=table.gene.density$chr.label))+
  geom_bar(stat="identity")
plot.gene.density







###NUMBER OF PEAKS PER CHROMOSOME
my.data <- read.delim("number_peaks_per_chromosome.txt", header=FALSE)
head(my.data)
class (my.data$V1)
str(my.data$V1)
class (my.data$V2)
str(my.data$V2)
my.data$V1<-factor (my.data$V1, levels = my.data$V1) #I am telling the levels of the factor being ordered in an specific order as they are in the col V1 of my data
my.data$V1 #level orders of the column have now changed

#To normalise to total number of overlapping genes per chromosome
my.data$V2
my.data$V2 <- my.data$V2/gene.count$Freq
my.data$V2


#Bar graph of chromosomes and number of peaks per chromosome
#fill= to tell the criteria to colour the bars, in this case different colours for different chromosomes
#geom_bar (stat="identity") it means that the statistics showing is only the value of the data, if ="bin" it will calculate the frequencies of each type of data
#aes () for telling which data in each axis
#xlab and ylab for labels names
#theme (legend.title) to specify about legend title
#text.x = element_text(angle = 90, hjust = 1) to write the x axis labels in vertical
plot_peaks_per_chromosome <- ggplot(data=my.data, aes(x=my.data$V1,y=my.data$V2, fill=my.data$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of peaks/total number of genes") + ggtitle("Peaks per chromosome normalised to total number of genes per chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

plot_peaks_per_chromosome

pdf("Number of genes normalised peaks per chromosome.pdf", onefile = TRUE)
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

#To normalise to total number of overlapping genes per chromosome
my.data2$V2
my.data2$V2 <- my.data2$V2/gene.count$Freq
my.data2$V2

#Bar graph of chromosomes and number of peaks per chromosome
#fill= to tell the criteria to colour the bars, in this case different colours for different chromosomes
#geom_bar (stat="identity") it means that the statistics showing is only the value of the data, if ="bin" it will calculate the frequencies of each type of data
#aes () for telling which data in each axis
#xlab and ylab for labels names
#theme (legend.title) to specify about legend title
#text.x = element_text(angle = 90, hjust = 1) to write the x axis labels in vertical
plot_reads_per_entire_chromosome <- ggplot(data=my.data2, aes(x=my.data2$V1,y=my.data2$V2, fill=my.data2$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of reads/total number of genes") + ggtitle("Reads per entire chromosome normalised to total number of genes per chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

plot_reads_per_entire_chromosome

pdf("Normalised reads per chromosome entire chromosome.pdf", onefile = TRUE)
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

#To normalise to total number of overlapping genes per chromosome
my.data3$V2
my.data3$V2 <- my.data3$V2/gene.count$Freq
my.data3$V2

#Bar graph of chromosomes and number of peaks per chromosome
#fill= to tell the criteria to colour the bars, in this case different colours for different chromosomes
#geom_bar (stat="identity") it means that the statistics showing is only the value of the data, if ="bin" it will calculate the frequencies of each type of data
#aes () for telling which data in each axis
#xlab and ylab for labels names
#theme (legend.title) to specify about legend title
#text.x = element_text(angle = 90, hjust = 1) to write the x axis labels in vertical
plot_number_reads_mapping_to_peaks_per_chromosome <- ggplot(data=my.data3, aes(x=my.data3$V1,y=my.data3$V2, fill=my.data3$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of reads/total number of genes") + ggtitle("Reads mapping to peaks in each chromosome normalised to total number of genes per chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

plot_number_reads_mapping_to_peaks_per_chromosome

pdf("Number of genes normalised reads mapping to peaks in each chromosome.pdf", onefile = TRUE)
plot_number_reads_mapping_to_peaks_per_chromosome
dev.off() #close the pdf







###NORMALISATION BASED ON NUMBER OF BP (NON OVERLAPPING) MAPPING  TO GENES PER CHROMOSOME###
############################################################################################



###NUMBER OF PEAKS PER CHROMOSOME
my.data4 <- read.delim("number_peaks_per_chromosome.txt", header=FALSE)
head(my.data4)
class (my.data4$V1)
str(my.data4$V1)
class (my.data4$V2)
str(my.data4$V2)
my.data4$V1<-factor (my.data4$V1, levels = my.data4$V1) #I am telling the levels of the factor being ordered in an specific order as they are in the col V1 of my data
my.data4$V1 #level orders of the column have now changed

#To normalise to total number of overlapping genes per chromosome
my.data4$V2
my.data4$V2 <- my.data4$V2/gene_bp$x
my.data4$V2


#Bar graph of chromosomes and number of peaks per chromosome
#fill= to tell the criteria to colour the bars, in this case different colours for different chromosomes
#geom_bar (stat="identity") it means that the statistics showing is only the value of the data, if ="bin" it will calculate the frequencies of each type of data
#aes () for telling which data in each axis
#xlab and ylab for labels names
#theme (legend.title) to specify about legend title
#text.x = element_text(angle = 90, hjust = 1) to write the x axis labels in vertical
plot_peaks_per_chromosome <- ggplot(data=my.data4, aes(x=my.data4$V1,y=my.data4$V2, fill=my.data4$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of peaks/bp mapping to genes ") + ggtitle("Peaks per chromosome normalised to total number of bp mapping to genes per chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

plot_peaks_per_chromosome

pdf("BpInGene normalised peaks per chromosome.pdf", onefile = TRUE)
plot_peaks_per_chromosome
dev.off() #close the pdf


###NUMBER OF READS PER ENTIRE CHROMOSOME
my.data5 <- read.delim("reads_per_entire_chromosome.txt", header=FALSE)
head(my.data5)
class (my.data5$V1)
str(my.data5$V1)
class (my.data5$V2)
str(my.data5$V2)
my.data5$V1<-factor (my.data5$V1, levels = my.data5$V1) #I am telling the levels of the factor being ordered in an specific order as they are in the col V1 of my data
my.data5$V1 #level orders of the column have now changed

#To normalise to total number of overlapping genes per chromosome
my.data5$V2
my.data5$V2 <- my.data5$V2/gene_bp$x
my.data5$V2

#Bar graph of chromosomes and number of peaks per chromosome
#fill= to tell the criteria to colour the bars, in this case different colours for different chromosomes
#geom_bar (stat="identity") it means that the statistics showing is only the value of the data, if ="bin" it will calculate the frequencies of each type of data
#aes () for telling which data in each axis
#xlab and ylab for labels names
#theme (legend.title) to specify about legend title
#text.x = element_text(angle = 90, hjust = 1) to write the x axis labels in vertical
plot_reads_per_entire_chromosome <- ggplot(data=my.data5, aes(x=my.data5$V1,y=my.data5$V2, fill=my.data5$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of reads/bp mapping to genes") + ggtitle("Reads per entire chromosome normalised to total number of bp mapping to genes per chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

plot_reads_per_entire_chromosome

pdf("BpInGene normalised reads per chromosome entire chromosome.pdf", onefile = TRUE)
plot_reads_per_entire_chromosome
dev.off() #close the pdf



###NUMBER OF READS MAPPING TO ALL PEAKS IN EACH CHROMOSOME
my.data6 <- read.delim("sum_number_of_reads_mapped_to_all_peaks_per_chromosome.txt", header=FALSE)
head(my.data6)
dim(my.data6)
class (my.data6$V1)
str(my.data6$V1)
class (my.data6$V2)
str(my.data6$V2)
my.data6$V1<-factor (my.data6$V1, levels = my.data6$V1) #Function factor () to tell the levels of the factor being ordered in an specific order as they are in the col V1 of my data
my.data3$V1 #level orders of the column have now changed

#To normalise to total number of overlapping genes per chromosome
my.data6$V2
my.data6$V2 <- my.data3$V2/gene_bp$x
my.data6$V2

#Bar graph of chromosomes and number of peaks per chromosome
#fill= to tell the criteria to colour the bars, in this case different colours for different chromosomes
#geom_bar (stat="identity") it means that the statistics showing is only the value of the data, if ="bin" it will calculate the frequencies of each type of data
#aes () for telling which data in each axis
#xlab and ylab for labels names
#theme (legend.title) to specify about legend title
#text.x = element_text(angle = 90, hjust = 1) to write the x axis labels in vertical
plot_number_reads_mapping_to_peaks_per_chromosome <- ggplot(data=my.data6, aes(x=my.data6$V1,y=my.data6$V2, fill=my.data6$V1)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of reads/bp mapping to genes") + ggtitle("Reads mapping to peaks in each chromosome normalised to total number of bp mapping to genes per chromosome")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

plot_number_reads_mapping_to_peaks_per_chromosome

pdf("BpInGene normalised reads mapping to peaks in each chromosome.pdf", onefile = TRUE)
plot_number_reads_mapping_to_peaks_per_chromosome
dev.off() #close the pdf




###NORMALISATION BASED ON NUMBER OF PEAKS FROM EACH CHROMOSOME PER MAPPED READS TO THAT CHROMOSOME###
####################################################################################################

num_peaks_per_mapped_reads_per_chrom<- my.data$V2/my.data2$V2
table.peaks.per.mapped.reads<- data.frame(cbind(my.data$V1, num_peaks_per_mapped_reads_per_chrom))
table.peaks.per.mapped.reads$num_peaks_per_mapped_reads_per_chrom<- as.numeric(table.peaks.per.mapped.reads$num_peaks_per_mapped_reads_per_chrom)
table.peaks.per.mapped.reads$V1<-factor (table.peaks.per.mapped.reads$V1, levels = table.peaks.per.mapped.reads$V1)
table.peaks.per.mapped.reads$V1
plot_num_peaks_per_mapped_reads_per_chrom<-ggplot(data=table.peaks.per.mapped.reads, aes(x=table.peaks.per.mapped.reads$V1,y=table.peaks.per.mapped.reads$num_peaks_per_mapped_reads_per_chrom, fill=table.peaks.per.mapped.reads$V1)) 
+ geom_bar(stat="identity")+xlab("Chromosome") + ylab("Peaks/total mapping reads") + ggtitle("Number of peaks normalised to number of total reads mapping to each chromosome")
+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))






