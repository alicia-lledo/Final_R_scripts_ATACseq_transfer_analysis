getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/pipeline_output_all_samples_QC/")
getwd ()

# Use ggplot library to draw graphs
library(ggplot2)

###PEAK OVERLAP UNINVOLVED AND LESIONAL

install.packages("VennDiagram")
library(VennDiagram)

grid.newpage()
vennplot <- draw.pairwise.venn(area1=37616, area2=36211, cross.area=21721, 
                               category=c("Uninvolved", "Lesional"), scaled=TRUE, 
                               col=c("darkgoldenrod", "steelblue2"), fill=c("darkgoldenrod", "steelblue2"), 
                               alpha=c(0.3, 0.2), label.col=rep("white", 3))
grid.draw(vennplot)




###PEAK OVERLAP UNINVOLVED AND LESIONAL VERSUS NHEKS (ROADMAP)

install.packages("VennDiagram")
library(VennDiagram)
#1=NHEKS, 2=uninvolved and 3=lesional

grid.newpage()
vennplot.UW <- draw.triple.venn(area1=290414, area2=37616, area3 =36211,  n12 = 25182, n23= 21721,n13= 25707,n123= 17020, 
                                category=c("NHEKs","uninvolved","lesional"), scaled=TRUE, 
                                col=c("lightsalmon","darkgoldenrod", "steelblue2"), fill=c("lightsalmon","darkgoldenrod", "steelblue2"),
                                alpha=c(0.3, 0.2, 0.3), label.col=rep("white", 7))
grid.draw(vennplot.UW)

#In percentage is an overlap of 66.9% for uninvolved and 70.9% for lesional





###DISTRIBUTION OF READS IN ALL PEAKS PER CHROMOSOME###


###READ COUNT FOR EACH PEAK (WHOLE GENOME)

uninvolved_peak_name_and_read_count <- read.delim ("uninvolved_bwa/uninvolved_peaks_name_and_read_count", header=FALSE)
lesional_peak_name_and_read_count <- read.delim ("lesional_bwa/lesional_peaks_name_and_read_count", header=FALSE)



#Generate a boxplot plot for the counts in the tables showing counts of reads per peak

uninvolved_peak_name_and_read_count$V3 <- "uninvolved"
lesional_peak_name_and_read_count$V3 <- "lesional"

hist (x=uninvolved_peak_name_and_read_count$V2)
hist(lesional_peak_name_and_read_count$V2)

merged_data <- rbind (uninvolved_peak_name_and_read_count, lesional_peak_name_and_read_count)
#merged_data <- rbind (uninvolved_peak_name_and_read_count, lesional_peak_name_and_read_count[which(nature_original_peak_name_and_read_count$V2<1000),])

merged.plot<-ggplot(merged_data, aes(x=merged_data$V3, y=merged_data$V2, fill=merged_data$V3)) + geom_boxplot()+
  theme_bw()+theme (axis.title.x = element_blank())+ ylab("Fragment counts in peaks")+ theme(plot.title = element_text(size=20,face="bold"),
                                                                                         axis.text.x = element_text(angle = 0, face="bold", size=18),
                                                                                         axis.text.y = element_text(face="bold",size=18), axis.title.y = element_text(face="bold",size=18))+scale_fill_manual(values=c("#999999", "#E69F00"), 
                                                                                                                                                                                                              name="Data source",labels=c("uninvolved","lesional"))

merged.plot

mean (uninvolved_peak_name_and_read_count$V2)
mean (lesional_peak_name_and_read_count$V2)





###PEAK ANNOTATION CHROMATIN STATES UNINVOLVED AND LESIONAL

peak.annotation<- read.delim("uninvolved_bwa/peak_annotation_chr_states_uninvolved_lesional.txt", header=TRUE)
###install.packages("reshape2")


ggplot(peak.annotation, stat="identity", aes(x=reorder(annotation,fold), y=peak.annotation$fold, fill=peak.annotation$track)) +  
  stat_identity(position="dodge", geom="bar")+coord_flip()+scale_fill_manual(values=c("#999999", "#E69F00"),guide = guide_legend(reverse=TRUE))+theme_bw()




###CORRELATION BETWEEN JKNIGHT AND GREENLEAF READ COUNT IN PEAKS

getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/QC_plots/pipeline_output_all_samples_QC/uninvolved_bwa/")
getwd ()

# Use ggplot library to draw graphs
library(ggplot2)

#read tables

correlation.peaks.counts<- read.delim("correlation_read_counts_skin.bed", header=FALSE)
correlation.peaks.counts$V1<-"uninvolved"
correlation.peaks.counts$V3<-"lesional"

##Not normalised to total counts
#correlation.peaks.counts$V2<- log10((correlation.peaks.counts$V2*2)) #to transform to number of reads instead of fragments
#correlation.peaks.counts$V4<- log10((correlation.peaks.counts$V4*2))

##Normalised to totalc counts in each dataset
correlation.peaks.counts$V2<- log10((correlation.peaks.counts$V2*2)/sum(correlation.peaks.counts$V2)) #to transform to number of reads instead of fragments
correlation.peaks.counts$V4<- log10((correlation.peaks.counts$V4*2)/sum(correlation.peaks.counts$V4))



#To remove inf values after having corrected for total number of reads in each dataset

correlation.peaks.counts$V2[which(correlation.peaks.counts$V2=="Inf")]<-NA 
correlation.peaks.counts$V2[which(correlation.peaks.counts$V2=="-Inf")]<-NA 

correlation.peaks.counts$V4[which(correlation.peaks.counts$V4=="Inf")]<-NA 
correlation.peaks.counts$V4[which(correlation.peaks.counts$V4=="-Inf")]<-NA 


#Colouring
# Correlation plot colours

similarity <- abs((correlation.peaks.counts$V2/correlation.peaks.counts$V4))
my.colours <- colorRampPalette(c("aliceblue", "darkblue"))
colour.vector <- my.colours(100)[as.numeric(cut(similarity, breaks=100))]


plot (x=correlation.peaks.counts$V2, y=correlation.peaks.counts$V4, pch=16, col=colour.vector, 
      xlab="ATAC-seq uninvolved (log10 reads)", ylab="ATAC-seq lesional (log10 reads)",xlim=c(-5.5,-3), ylim = c(-5.5,-3))
abline(0, 1)



testing1 <-cor.test(correlation.peaks.counts$V2, correlation.peaks.counts$V4)
testing2 <- cor.test(correlation.peaks.counts$V2, correlation.peaks.counts$V4, exact=TRUE, method=c("kendall"))

testing1
testing2

