getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/Mapping")
getwd ()

# Use ggplot library to draw graphs
library(ggplot2)

# read in count data for insert size
my.data <- read.delim("fragment_length_and_frquency_of_reads_my_data_CD14.txt", header=TRUE)
# Convert to density (normalisation) by dividing by total number of counts
my.data$All_Reads.fr_count <- my.data$All_Reads.fr_count/sum(my.data$All_Reads.fr_count)
#Create another column which contains the same normalised densities multiplied *1000
my.data$Normalised.read.density.1000 <- my.data$All_Reads.fr_count*1000
head (my.data)


# Plot insert size on x axis, normalised read densityx10^-3 total on y
# i.e. a density plot of insert size
#aes to specify axis parametres (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
#by indicating discrete it is like a categorical vector and therefore it plots every single value
fig2a <- ggplot(my.data, aes(x=insert_size, y=Normalised.read.density.1000)) + theme_bw() + 
  geom_line(stat="identity", colour = "red") + scale_x_continuous(name="Fragment length (bp)")+ 
  scale_y_continuous(name="Normalised insert size density (x10^-3)") + ggtitle("Normalised insert size density")+
  theme(plot.title = element_text(size=20,face="bold"), 
                                                                     legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), 
                                                                     axis.title.x = element_text(face="bold",size=20),
                                                                     axis.text.y = element_text(hjust = 1, size=15), 
                                                                     axis.title.y = element_text(face="bold",size=18))# Don't add shaded confidence region
fig2a #to make it appear in the graph screen
pdf("Normalised insert size density (bwa).pdf", onefile = TRUE)
fig2a
dev.off() #close the pdf




#geom_point(stat="identity", colour = "red") if I wanted to plot each data point as a dot and then link them with the line

# Plot insert size on x axis, normalised read density
# y axis on log scale
fig2alog <-ggplot(my.data, aes(x=insert_size, y=All_Reads.fr_count)) + theme_bw() +
  geom_line(stat="identity", colour = "red") + scale_y_log10(name="Normalised read density")+ scale_x_continuous(name="Fragment length (bp)")+
  ggtitle("Normalised insert size density (y axis in log10 scale")+
  theme(plot.title = element_text(size=20,face="bold"), 
        legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15), 
        axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15), 
        axis.title.y = element_text(face="bold",size=18))# Don't add shaded confidence region
fig2alog
pdf("Normalised insert size density (y axis in log10 scale(bwa).pdf", onefile = TRUE)
fig2a
fig2alog
dev.off() #close the pdf

#geom_point(stat="identity", colour = "red") includes points in the plot