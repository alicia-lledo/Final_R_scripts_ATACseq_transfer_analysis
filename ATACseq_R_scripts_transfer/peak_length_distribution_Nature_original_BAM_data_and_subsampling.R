getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/Peak_calling")
getwd()



#To generate histograms of the Nature data peak calling using the original BAM file not fastq processed:

##Nature original BAM data no downsampling no --call-summits:

Nature.original.BAM.data.no.callsummits.table <- read.delim ("Nature_original_BAM_file_no_callsummits_peak_length.txt", header=FALSE)
head (Nature.original.BAM.data.no.callsummits.table)
dim (Nature.original.BAM.data.no.callsummits.table)
class (Nature.original.BAM.data.no.callsummits.table)
max(Nature.original.BAM.data.no.callsummits.table$V1) #report maximum value in the data frame column V1
min(Nature.original.BAM.data.no.callsummits.table$V1) #report minimum value in the data frame column V1
median (Nature.original.BAM.data.no.callsummits.table$V1)
x_axis_breaks_Nature_original.data <- seq(from = 200, to = 61000, by = 400) #numeric vector that starts at 200 and reaches 60400 by spacing the values every 400
#par(mfrow = c(1,2))
#par(mar = rep(2, 4))
par(lab = (c(21, 10, 3)))
hist (Nature.original.BAM.data.no.callsummits.table$V1, breaks=x_axis_breaks_Nature_original.data, col="blue", xlim=c(200,61000), ylim=c(0,27000), main="Histogram of Nature original BAM data peaks no --call-summits", xlab="length (bp)", ylab="Number of peaks")

#Remove values higher than 4035 and repeat the plot
#Use subset function
Nature_original_selected_values <- subset (Nature.original.BAM.data.no.callsummits.table, Nature.original.BAM.data.no.callsummits.table$V1 < 4035)
head (Nature_original_selected_values)
class (Nature_original_selected_values)
dim (Nature_original_selected_values)
length (which (Nature.original.BAM.data.no.callsummits.table$V1 > 4035))
par(lab = (c(21, 10, 3)))
hist (Nature_original_selected_values$V1, breaks=x_axis_breaks_mydata, col="blue", xlim=c(200,4200), ylim=c(0,25000), main="Histogram of Nature original BAM data peaks no --call-summits", xlab="length (bp)", ylab="Number of peaks")


##Nature original BAM data downsampling --call-summits:
Nature.original.BAM.data.subsampling.callsummits.table <- read.delim ("Nature_CD14_normalised_peak_length.txt", header=FALSE)
head (Nature.original.BAM.data.subsampling.callsummits.table)
dim (Nature.original.BAM.data.subsampling.callsummits.table)
class (Nature.original.BAM.data.subsampling.callsummits.table)
max(Nature.original.BAM.data.subsampling.callsummits.table$V1) #report maximum value in the data frame column V1
min(Nature.original.BAM.data.subsampling.callsummits.table$V1) #report minimum value in the data frame column V1
median (Nature.original.BAM.data.subsampling.callsummits.table$V1)
#x_axis_breaks_Nature_original.data <- seq(from = 200, to = 61000, by = 400) #numeric vector that starts at 200 and reaches 60400 by spacing the values every 400
#par(mfrow = c(1,2))
#par(mar = rep(2, 4))
par(lab = (c(21, 10, 3)))
hist (Nature.original.BAM.data.subsampling.callsummits.table$V1, breaks=200, col="blue", xlim=c(200,65000), ylim=c(0,18000), main="Histogram of Nature original BAM data downsampling peaks --call-summits", xlab="length (bp)", ylab="Number of peaks")


#Remove values higher than 4035 and repeat the plot
#Use subset function
Nature_original_subsampling_selected_values <- subset (Nature.original.BAM.data.subsampling.callsummits.table, Nature.original.BAM.data.subsampling.callsummits.table$V1 < 4035)
head (Nature_original_subsampling_selected_values)
class (Nature_original_subsampling_selected_values)
dim (Nature_original_subsampling_selected_values)
length (which (Nature_original_subsampling_selected_values$V1 > 4035))#from all the data only 249 peaks have a length greater than 4035
par(lab = (c(21, 10, 3)))
x_axis_breaks_mydata <- seq(from = 200, to = 4200, by = 200) #numeric vector that starts at 200 and reaches 4200 by spacing the values every 200
hist (Nature_original_subsampling_selected_values$V1, breaks=x_axis_breaks_mydata, col="blue", xlim=c(200,4200), ylim=c(0,18000), main="Histogram of Nature original BAM data peaks --call-summits", xlab="length (bp)", ylab="Number of peaks")

