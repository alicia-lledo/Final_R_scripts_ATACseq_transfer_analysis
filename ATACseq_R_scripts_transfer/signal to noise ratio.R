getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data analysis/QC_plots/QC_files") 
getwd ()


library(ggplot2)

#reads_mapping_in_TssA <- read.delim("reads_mapping_to_TssA_CD14.txt", header=FALSE)
#reads_not_mapping_in_TssA <- read.delim("reads_not_mapping_to_TssA_CD14.txt", header=FALSE)

reads_per_chrom <- read.delim("correctreads_per_entire_chromosome.txt", header=FALSE)
per_chr_reads_mapping_in_TssA <- read.delim("per_chr_reads_mapping_to_TssA_CD14.txt", header=FALSE)
per_chr_reads_not_mapping_in_TssA <- read.delim("per_chr_reads_not_mapping_to_TssA_CD14.txt", header=FALSE)

#To express results as number of reads
reads_per_chrom$V2[1:23] <- (reads_per_chrom$V2[1:23])*2
per_chr_reads_mapping_in_TssA$V2[1:23] <- (per_chr_reads_mapping_in_TssA$V2[1:23])*2                        
per_chr_reads_not_mapping_in_TssA$V2[1:23] <- (per_chr_reads_not_mapping_in_TssA$V2[1:23])*2

#To creat one column with chr name column in appropriate order
ordered.chromosomes <- factor(x=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"), levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))

###CHROMOSOME size

chr.size <- read.delim ("hg19.chrom.sizes.UCSC.txt", header=FALSE)
levels(chr.size$V1)[23] <-"chr23" #to change the name of a particular level of a factor
chr.size$V2<- chr.size$V2/1000000 #chr size in Mb



#####################################################################################
##################SIGNAL TO NOISE RATIOS############################################
####################################################################################

#1. Number of reads mapping to TssA divided by number of reads mapped to no TssA
TssA_vs_not_TssA <- data.frame(chromosome= ordered.chromosomes, Ratio= per_chr_reads_mapping_in_TssA$V2[1:23]/per_chr_reads_not_mapping_in_TssA$V2[1:23])
TssA_vs_not_TssA$Ratio <- TssA_vs_not_TssA$Ratio*100

ggplot(data=TssA_vs_not_TssA, aes(x=TssA_vs_not_TssA$chromosome,y=TssA_vs_not_TssA$Ratio, fill=TssA_vs_not_TssA$chromosome)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("# Reads in TssA/#reads not in TssA (x100)") + ggtitle("Reads in TssA per reads not in TssA")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

  
    ###CORRELATION WITH CHROMOSOME SIZE

#Create a data fram with 3 columns: 
#col1 chr name
#col2: chr size
#col3: TssA_vs_not_TssA ratio
correlation <- data.frame(chromosome_size=chr.size$V2[1:23], TssA_vs_not_TssA_ratio=TssA_vs_not_TssA$Ratio[1:23])

#Create the same data frame adding the chromosome name for ggplot2
correlation2 <- data.frame(chromosome=ordered.chromosomes,chromosome_size=chr.size$V2[1:23], TssA_vs_not_TssA_ratio=TssA_vs_not_TssA$Ratio[1:23])

#Use the cor(x,y) function to test for linear pearson correlation between two variables

cor(correlation, use="everything", method="pearson")
ggplot(correlation2, aes(x=correlation2$chromosome_size, y=correlation2$TssA_vs_not_TssA_ratio, color=correlation2$chromosome)) +
  xlab("Chromosome size (Mb)") + 
  ylab("TssA reads per not TssA reads (x100)") + ggtitle("Correlation between chromosome size and TssA/not TssA") +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))    # Don't add shaded confidence region


      ##It seems the best measure for signal to noise ratio as it looks even across different chromosomes and is not correlated with size




#2.Also divide reads mapped to TssA to total reads mapped to any annotation feature (PE reads 6,893,004 or 13,786,008 single reads)   
annotated_reads <-13786008
TssA_vs_all_annotated<- data.frame(chromosomes=ordered.chromosomes,TssA.vs.annotated.total=per_chr_reads_mapping_in_TssA$V2[1:23]/annotated_reads)
TssA_vs_all_annotated$TssA.vs.annotated.total <-TssA_vs_all_annotated$TssA.vs.annotated.total*100
ggplot(data=TssA_vs_all_annotated, aes(x=TssA_vs_all_annotated$chromosome,y=TssA_vs_all_annotated$TssA.vs.annotated.total, fill=TssA_vs_all_annotated$chromosome)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("#Reads in TssA/total reads annotated (x100)") + ggtitle("Reads in TssA per total annotated reads")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))



#3.Also divide reads to TssA by the total number of reads mapping to the entire chromosome.Total number of correct PE reads are 7,226,241 (14,452,482)
total_reads <- 14452482
TssA_vs_total_reads <- data.frame(chromosomes=ordered.chromosomes, TssA.vs.total=per_chr_reads_mapping_in_TssA$V2[1:23]/total_reads)
TssA_vs_total_reads$TssA.vs.total <- TssA_vs_total_reads$TssA.vs.total*100
ggplot(data=TssA_vs_total_reads, aes(x=TssA_vs_total_reads$chromosome,y=TssA_vs_total_reads$TssA.vs.total, fill=TssA_vs_total_reads$chromosome)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("#Reads in TssA/total reads (x100)") + ggtitle("Reads in TssA per total reads")+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))


