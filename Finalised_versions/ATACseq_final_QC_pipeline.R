getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/pipeline_output_all_samples_QC/correct_blood_bwa/")
getwd ()

# Use ggplot library to draw graphs
library(ggplot2)

####READING AND CREATING TABLES###

###READING COMBINED FILE
combined.table<- read.delim("non-lesional_new_filt.stat_new.txt", header=TRUE)

###CHROMOSOMES
chromosome.names <- combined.table[,1]
chromosome.names <-factor (chromosome.names, levels = chromosome.names)
chromosome.names


##NUMBER OF PEAKS PER CHROMOSOME
peak.per.chr <- combined.table[,2] #subset using only col name,without table
#peak.per.chr$Chromosome<-factor (peak.per.chr$Chromosome, levels = peak.per.chr$Chromosome) #I am telling the levels of the factor being ordered in an specific order as they are in the col V1 of my data
#peak.per.chr$Chromosome #level orders of the column have now changed

###CHROMOSOME SIZE UNGAPPED
chr.size.ungapped <- read.delim ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/QC_files/ungapped_chromosome_size.txt", header=FALSE)
#levels(chr.size.ungapped$V1)[23] <-"chr23" #to change the name of a particular level of a factor
chr.size.ungapped$V2<- chr.size.ungapped$V2/1000000 #chr size in Mb


###NUMBER OF READS PER ENTIRE CHROMOSOME
reads.per.chr <- combined.table[,4]
reads.per.chr

#reads.per.chr$Chromosome<-factor (reads.per.chr$Chromosome, levels = reads.per.chr$Chromosome) #I am telling the levels of the factor being ordered in an specific order as they are in the col V1 of my data
#reads.per.chr$Chromosome #level orders of the column have now changed


####NUMBER OF READS MAPPING TO ALL PEAKS IN EACH CHROMOSOME
reads.mapping.peaks.per.chr <- combined.table[,3]
reads.mapping.peaks.per.chr
#reads.mapping.peaks.per.chr$Chromosome<-factor (reads.mapping.peaks.per.chr$Chromosome, levels = reads.mapping.peaks.per.chr$Chromosome) #Function factor () to tell the levels of the factor being ordered in an specific order as they are in the col V1 of my data
#reads.mapping.peaks.per.chr$Chromosome #level orders of the column have now changed


####NUMBER OF READS MAPPING OFF PEAKS (same length) IN EACH CHROMOSOME
reads.mapping.offpeaks.same.length.per.chr <- combined.table[,5]
reads.mapping.offpeaks.same.length.per.chr
#reads.mapping.offpeaks.same.length.per.chr$V1<-factor (reads.mapping.offpeaks.same.length.per.chr$V1, levels = reads.mapping.offpeaks.same.length.per.chr$V1) #Function factor () to tell the levels of the factor being ordered in an specific order as they are in the col V1 of my data
#reads.mapping.offpeaks.same.length.per.chr$V1 #level orders of the column have now changed


###PEAK WIDTH
peak.width.table <- read.delim("./blood_bwa.peaks_width.hist", header=FALSE)
min(peak.width.table$V1)
max(peak.width.table$V1)
median (peak.width.table$V1)


###INSERT SIZE

# read in count data for insert size
insert.size.table <- read.delim("length_table.txt", header=TRUE)  #output from InsertSizeMetrics

# Convert to density (normalisation) by dividing by total number of counts

insert.size.table$All_Reads.fr_count <- insert.size.table$All_Reads.fr_count/sum(insert.size.table$All_Reads.fr_count)

#Create another column which contains the same normalised densities multiplied *1000
insert.size.table$Normalised.read.density.1000 <- insert.size.table$All_Reads.fr_count*1000






#########################PLOTS########
#####################################



###DENSITY PLOTS (PERIODICITY)

##Fig2a Nature paper
# Plot insert size on x axis, normalised read densityx10^-3 total on y
# i.e. a density plot of insert size
#aes to specify axis parametres (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
#by indicating discrete it is like a categorical vector and therefore it plots every single value
plot1 <- ggplot(insert.size.table, aes(x=insert_size, y=Normalised.read.density.1000)) + theme_bw() + 
  geom_line(stat="identity", colour = "red", alpha=1) + scale_x_continuous(name="Fragment length (bp)", breaks=seq(0,800,by=200))+ 
  scale_y_continuous(name="Normalised insert size density (x10^-3)") + ggtitle("Normalised insert size density uninvolved skin")+
  theme(plot.title = element_text(size=20,face="bold"), 
        legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), 
        axis.title.x = element_text(face="bold",size=18),
        axis.text.y = element_text(hjust = 1, size=15), 
        axis.title.y = element_text(face="bold",size=18))# Don't add shaded confidence region
plot1 #to make it appear in the graph screen



##Fig2a inset Nature paper
#geom_point(stat="identity", colour = "red") if I wanted to plot each data point as a dot and then link them with the line
# Plot insert size on x axis, normalised read density
# y axis on log scale
plot2 <-ggplot(insert.size.table, aes(x=insert_size, y=All_Reads.fr_count)) + theme_bw() +
  geom_line(stat="identity", colour = "red") + scale_y_log10(name="Normalised read density")+ scale_x_continuous(name="Fragment length (bp)", breaks=seq(0,1000,by=200))+
  ggtitle("Normalised insert size density\n  uninvolved skin (y axis in log10 scale)")+
  theme(plot.title = element_text(size=20,face="bold"), 
        legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), 
        axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), 
        axis.title.y = element_text(face="bold",size=18))# Don't add shaded confidence region
plot2




###PEAK WIDTH DISTRIBUTION

min (peak.width.table$V1)
max (peak.width.table$V1)

x_axis_breaks <- seq(from =min(peak.width.table$V1), to = 2000, by = 100)

length(peak.width.table$V1[which(peak.width.table$V1<2000)])



plot3<-hist (peak.width.table$V1[which(peak.width.table$V1<2000)], breaks=x_axis_breaks,col="white", main="Peak width distribution", 
                                          xlab="length (bp)", ylab="Number of peaks", xaxp=c(200,2000,9))
#xlim=c(200,2000) it only indicates start and end not sequence
plot3
peak.width.table$V2<-"JKnight lab"

plot3b<-plot<-ggplot(peak.width.table, aes(x=peak.width.table$V1, y=peak.width.table$V2, fill=peak.width.table$V2)) + geom_boxplot() +
  guides(fill=FALSE)+ theme_bw()+theme (axis.title.x = element_blank())+ ylab("Peak width (bp)")+ theme(plot.title = element_text(size=20,face="bold"), 
                                                                                                        axis.text.y = element_text(face="bold",size=18), 
                                                                                                        axis.title.y = element_text(face="bold",size=18),
                                                                                                        legend.title=element_blank(), axis.text.x = element_text(angle = 0, face="bold", size=18))
plot3b  



###PLOT OF DISTRIBUTION OF NUMBER OF PEAKS PER CHROMOSOME

plot4<- ggplot(data=combined.table, aes(x=combined.table$Chromosome,y=combined.table$Called_peaks, fill=combined.table$Chromosome)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of peaks") + ggtitle("Number of peaks per chromosome")+ 
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))
plot4



#CORRELATION OF NUMBER OF PEAKS PER CHROMOSOME AND UNGAPPED CHROMOSOME SIZE

correlation.peaksperchr.and.chr.size <- data.frame(chromosome_size=chr.size.ungapped$V2, peaks.per.chr=combined.table$Called_peaks)

#Create the same data frame adding the chromosome name for ggplot2
correlation.peaksperchr.and.chr.size2 <- data.frame(chromosome=chromosome.names, chromosome_size=chr.size.ungapped$V2, peaks.per.chr=combined.table$Called_peaks)
#correlation.peaksperchr.and.chr.size2$chromosome <-factor (correlation.peaksperchr.and.chr.size2$chromosome, levels=correlation.peaksperchr.and.chr.size2$chromosome)
#Use the cor(x,y) function to test for linear pearson correlation between two variables
pval<-cor.test(correlation.peaksperchr.and.chr.size$chromosome_size, correlation.peaksperchr.and.chr.size$peaks.per.chr)
plot5<-ggplot(correlation.peaksperchr.and.chr.size2, aes(x=correlation.peaksperchr.and.chr.size2$chromosome_size, 
                                                         y=correlation.peaksperchr.and.chr.size2$peaks.per.chr, color=correlation.peaksperchr.and.chr.size2$chromosome)) +
  xlab("Chromosome size (Mb)") + 
  ylab("Number of peaks per chromosome") + ggtitle("Correlation between chromosome size\n and number of peaks per chr") +
  geom_point(size=3) +    # Use hollow circles geom_point(shape=1) empty for filled
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE, aes (group=1))+ # Don't add shaded confidence region
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))

plot5





###Q1.NORMALISATION BASED ON NUMBER OF PEAKS FROM EACH CHROMOSOME PER MAPPED READS TO THAT CHROMOSOME###
#####################################################################################################


table.num_peaks_per_mapped_reads_per_chrom<- data.frame (chromosome=chromosome.names, peaks.per.reads=((peak.per.chr/reads.per.chr)*100))

#table.num_peaks_per_mapped_reads_per_chrom[24,2]=0 #replace missing value in chrY
plot6<- ggplot(data=table.num_peaks_per_mapped_reads_per_chrom, aes(x=table.num_peaks_per_mapped_reads_per_chrom$chromosome,
                                                                    y=table.num_peaks_per_mapped_reads_per_chrom$peaks.per.reads, 
                                                                    fill=table.num_peaks_per_mapped_reads_per_chrom$chromosome))+
  geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of peaks per number of PE reads (%)") + 
  ggtitle("Number of peaks per number of\n PE reads (%) in each chromosome\n")+ 
  theme(plot.title = element_text(size=16,face="bold"), legend.title=element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.title.x = element_text(face="bold",size=15),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=14))
plot6




###Q2.NUMBER OF READS MAPPING TO PEAKS PER NUMBER OF TOTAL READS MAPPING TO THAT CHROMOSOME###
#############################################################################################################


table.num_reads_in_peaks_per_total_reads_in_chr<-data.frame (chromosome=chromosome.names, reads.in.peaks.per.total.reads=((reads.mapping.peaks.per.chr/reads.per.chr)*100))

#table.num_reads_in_peaks_per_total_reads_in_chr[24,2]=0
plot7<- ggplot(data=table.num_reads_in_peaks_per_total_reads_in_chr, aes(x=table.num_reads_in_peaks_per_total_reads_in_chr$chromosome,
                                                                         y=table.num_reads_in_peaks_per_total_reads_in_chr$reads.in.peaks.per.total.reads, 
                                                                         fill=table.num_reads_in_peaks_per_total_reads_in_chr$chromosome)) + geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Number of PE reads in peaks per total number fragments (%)") + 
  ggtitle("Number of PE reads in peaks per number of\n PE reads in each chromosome")+theme_bw()+
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))
plot7

#var(data)=42.9

###Q4. NUMBER OF READS IN PEAKS PER NUMBER OF READS OFF-PEAKS ( non target regions of the same length that peaks) PER CHROMOSOME###
##############################################################################################################################

table.reads.in.peaks.per.reads.off.peak.same.length<-data.frame (chromosome=chromosome.names, 
                                                                 ratio.peaks.off.peaks=(reads.mapping.peaks.per.chr/reads.mapping.offpeaks.same.length.per.chr))
#var(data)

#table.reads.in.peaks.per.reads.off.peak.same.length<-data.frame (chromosome=chromosome.names,reads.in.peaks.per.reads.off.peaks=(reads.mapping.peaks.per.chr/reads.mapping.offpeaks.same.length.per.chr))

#table.reads.in.peaks.per.reads.off.peak.same.length[24,2]=0 #replace missing value in chrY


plot8<- ggplot(data=table.reads.in.peaks.per.reads.off.peak.same.length, aes(x=table.reads.in.peaks.per.reads.off.peak.same.length$chromosome,
                                                                             y=table.reads.in.peaks.per.reads.off.peak.same.length$ratio.peaks.off.peaks, 
                                                                             fill=table.reads.in.peaks.per.reads.off.peak.same.length$chromosome)) + theme_bw()+geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Fold enrichment of fragments on peaks versus off peaks")+
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))+ scale_y_continuous(breaks=seq(0, 3, 0.5))+expand_limits(y=c(0,3))

plot8
mean (table.reads.in.peaks.per.reads.off.peak.same.length$ratio.peaks.off.peaks)

#Percentage version


readsQ4<- data.frame (reads.in.peaks=combined.table$Fragments_in_peaks, reads.off=combined.table$Fragments_per_off_peaks_peak_length)
readsQ4$total<-rowSums(readsQ4)

table.reads.in.peaks.per.reads.off.peak.same.length<-data.frame (chromosome=chromosome.names, 
                                                                 percentage.reads.in.peaks=((reads.mapping.peaks.per.chr)*100)/readsQ4$total)
plot9<- ggplot(data=table.reads.in.peaks.per.reads.off.peak.same.length, aes(x=table.reads.in.peaks.per.reads.off.peak.same.length$chromosome,
                                                                         y=table.reads.in.peaks.per.reads.off.peak.same.length$percentage.reads.in.peaks, 
                                                                         fill=table.reads.in.peaks.per.reads.off.peak.same.length$chromosome)) + theme_bw()+geom_bar(stat="identity")+
  xlab("Chromosome") + ylab("Percentage of reads in peaks")+theme_bw()+ scale_y_continuous(breaks=seq(0, 100, 20))+
  theme(plot.title = element_text(size=20,face="bold"), legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), axis.title.y = element_text(face="bold",size=18))

plot9

#ggtitle("Fold enrichment of fragments in peaks\n versus off-peak regions")

######



###CORRELATION BETWEEN JKNIGHT AND GREENLEAF READ COUNT IN PEAKS

getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/Peak_correlation_jknight_greenleaf/")
getwd ()

# Use ggplot library to draw graphs
library(ggplot2)

#read tables

correlation.peaks.counts<- read.delim("correlation_read_counts_jknight_nature.bed", header=FALSE)
correlation.peaks.counts$V1<-"Jknight"
correlation.peaks.counts$V3<-"Greenleaf"

##Not normalised to total counts
#correlation.peaks.counts$V2<- log10((correlation.peaks.counts$V2*2)) #to transform to number of reads instead of fragments
#correlation.peaks.counts$V4<- log10((correlation.peaks.counts$V4*2))

##Normalised to totalc counts in each dataset
correlation.peaks.counts$V2<- log10((correlation.peaks.counts$V2*2)/sum(correlation.peaks.counts$V2)) #to transform to number of reads instead of fragments
correlation.peaks.counts$V4<- log10((correlation.peaks.counts$V4*2)/sum(correlation.peaks.counts$V4))

#Colouring
# Correlation plot colours

similarity <- abs((correlation.peaks.counts$V2/correlation.peaks.counts$V4))
my.colours <- colorRampPalette(c("aliceblue", "darkblue"))
colour.vector <- my.colours(100)[as.numeric(cut(similarity, breaks=100))]


#To colour darker those peaks with strongest correlation
#my.colours <- rep("steelblue", length(correlation.peaks.counts$V2))
#for (i in 1: length(my.colours)){
 # if(correlation.peaks.counts$V2[i]/correlation.peaks.counts$V4[i] >= 0.95 &
  #   correlation.peaks.counts$V2[i]/correlation.peaks.counts$V4[i] <= 1.05)
   # my.colours[i] <- "darkblue"
#}

#To remove inf values after having corrected for total number of reads in each dataset

correlation.peaks.counts$V2[which(correlation.peaks.counts$V2=="Inf")]<-NA 
correlation.peaks.counts$V2[which(correlation.peaks.counts$V2=="-Inf")]<-NA 

correlation.peaks.counts$V4[which(correlation.peaks.counts$V4=="Inf")]<-NA 
correlation.peaks.counts$V4[which(correlation.peaks.counts$V4=="-Inf")]<-NA 


plot (x=correlation.peaks.counts$V2, y=correlation.peaks.counts$V4, pch=16, col=colour.vector, 
      xlab="ATAC-seq JKnight (log10 reads)", ylab="ATAC-seq Greenleaf (log10 reads)",xlim=c(-5.5,-3), ylim = c(-5.5,-3))
abline(0, 1)



testing1 <- cor.test(correlation.peaks.counts$V2, correlation.peaks.counts$V4)
testing2 <- cor.test(correlation.peaks.counts$V2, correlation.peaks.counts$V4, exact=FALSE, method=c("spearman"))

testing1
testing2


###PEAK OVERLAP JKNIGHT GREENLEAF VENN DIAGRAM

install.packages("VennDiagram")
library(VennDiagram)

grid.newpage()
vennplot <- draw.pairwise.venn(area1=45940, area2=36761, cross.area=27903, 
                               category=c("JKnight lab", "Greenleaf lab"), scaled=TRUE, 
                               col=c("darkgoldenrod", "steelblue2"), fill=c("darkgoldenrod", "steelblue2"), 
                               alpha=c(0.3, 0.2), label.col=rep("white", 3))
grid.draw(vennplot)




###PEAK OVERLAP JKNIGHT GREENLEAF vs DUKE AND UW DHS VENN DIAGRAM

install.packages("VennDiagram")
library(VennDiagram)
#1=UW, 2=JKnight and 3=Greenleaf
#1=DUKE, 2=JKnight and 3=Greenleaf
grid.newpage()
vennplot.UW <- draw.triple.venn(area1=209204, area2=45940, area3 =36761,  n12 = 43558, n23= 27903,n13= 31270,n123= 25823, 
                               category=c("UW","JKnight lab","Greenleaf lab"), scaled=TRUE, 
                               col=c("lightsalmon","darkgoldenrod", "steelblue2"), fill=c("lightsalmon","darkgoldenrod", "steelblue2"),
                              alpha=c(0.3, 0.2, 0.3), label.col=rep("white", 7))
grid.draw(vennplot.UW)

#In percentage is an overlap of 94.8% for my data and 85% for nature


grid.newpage()
vennplot.Duke <- draw.triple.venn(area1=159457, area2=45940, area3 =36761,  n12 = 29796, n23= 27903,n13= 24370,n123= 21547, 
                                category=c("Duke","JKnight lab","Greenleaf lab"), scaled=TRUE, 
                                col=c("darkolivegreen3","darkgoldenrod", "steelblue2"), fill=c("darkolivegreen3","darkgoldenrod", "steelblue2"),
                                alpha=c(0.3, 0.2, 0.3), label.col=rep("white", 7))
grid.draw(vennplot.Duke)


#In percentage is an overlap of 64.85% for my data and 66.29% for nature



###PEAK ANNOTATION CHROMATIN STATES JKNIGHT AND GREENLEAF DATASETS

peak.annotation<- read.delim("chr_states_peak_annotation_jknight_greenleaf.txt", header=TRUE)
###install.packages("reshape2")


ggplot(peak.annotation, stat="identity", aes(x=reorder(chrom_state,fold.enrichment), y=peak.annotation$fold.enrichment, fill=peak.annotation$Dataset)) +  
  stat_identity(position="dodge", geom="bar")+coord_flip()+scale_fill_manual(values=c("#999999", "#E69F00"),guide = guide_legend(reverse=TRUE))+theme_bw()




#PEAK ENRICHMENT OF TFBS

top_10_TFBS <- read.delim("top_10_TFBS_jknight.txt", header=TRUE)


ggplot(top_10_TFBS, stat="identity", aes(x=top_10_TFBS$TF, y=top_10_TFBS$fold_enrichment))+coord_flip()+
  geom_bar(aes(fill=top_10_TFBS$TF), stat="identity")+xlab("Fold enrichment") + ylab("TF")+theme_bw()

