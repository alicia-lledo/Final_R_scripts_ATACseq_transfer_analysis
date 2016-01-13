getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/Peak_calling")
getwd()


###Pulling the 2nd and 3rd columns with awk from the narrowPeak files:

##Nature data --call-summits
Nature.callsummits.table <- read.delim ("Nature_call_summits_length.txt", header=FALSE)
head (Nature.callsummits.table)
class (Nature.callsummits.table)
#Creates a numericvector to use as input for the hist break
x_axis_breaks_nature_right_scale <- seq(from = 200, to = 2300, by = 100) #numeric vector that starts at 200 and reaches 2300 by spacing the values every 100
x_axis_breaks_nature_changed_scale <- seq(from = 200, to = 4200, by = 100) #numeric vector that starts at 200 and reaches 2300 by spacing the values every 100
#par function includes many tasks and it is very useful for plot editing
par(mfrow = c(1,2))
#A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) 
par(mar = rep(4, 4)) #rep replicates the first element in this case 2 four times. It sets the margins of the plot so it fits.Otherwise I was getting an error
right.scale <-par(lab = (c(22, 10, 3))) #lab=c(x,y,n) The values of x and y give the (approximate) number of tickmarks on the x and y axes and len specifies the label length
histogram.right.scale.callsummits <-hist (Nature.callsummits.table$V1, breaks=x_axis_breaks_nature_right_scale, col="blue", xlim=c(200,2300), ylim=c(0,18000), main="Histogram of Nature data peak --call-summits", xlab="length (bp)", ylab="Number of peaks")
#To generate the same histogram but using the axis scale from my data histogram
changed.scale <-par(lab = (c(21, 10, 3))) #lab=c(x,y,n) The values of x and y give the (approximate) number of tickmarks on the x and y axes and len specifies the label length
histogram.changed.scale.callsummits <-hist (Nature.callsummits.table$V1, breaks=x_axis_breaks_nature_changed_scale, col="blue", xlim=c(200,4200), ylim=c(0,20000), main="Histogram of Nature data peak --call-summits changed scale", xlab="length (bp)", ylab="Number of peaks")
max(Nature.callsummits.table$V1) #report maximum value in the data frame column V1
min(Nature.callsummits.table$V1) #report minimum value in the data frame column V1
median (Nature.callsummits.table$V1)

##Nature data without --call-summits
Nature.no.callsummits.table <- read.delim ("Nature_no_call_summits_length.txt", header=FALSE)
head (Nature.no.callsummits.table)
dim (Nature.no.callsummits.table)
class (Nature.no.callsummits.table)
max(Nature.no.callsummits.table$V1) #report maximum value in the data frame column V1
min(Nature.no.callsummits.table$V1) #report minimum value in the data frame column V1
median (Nature.no.callsummits.table$V1)
par(mfrow = c(1,2))
par(mar = rep(4, 4))
right.scale 
histogram.right.scale.nocallsummits <-hist (Nature.no.callsummits.table$V1, breaks=x_axis_breaks_nature_right_scale, col="blue", xlim=c(200,2300), ylim=c(0,18000), main="Histogram of Nature data peak no --call-summits", xlab="length (bp)", ylab="Number of peaks")
#To generate the same histogram but using the axis scale from my data histogram
changed.scale 
histogram.changed.scale.nocallsummits <-hist (Nature.no.callsummits.table$V1, breaks=x_axis_breaks_nature_changed_scale, col="blue", xlim=c(200,4200), ylim=c(0,20000), main="Histogram of Nature data peak no --call-summits changed scale", xlab="length (bp)", ylab="Number of peaks")

#To create a plot with the Nature --call-summits and no --call-summits (in the right scale) side by side:
par(mfrow = c(1,2))
par(mar = rep(4, 4))
histogram.right.scale.callsummits <-hist (Nature.callsummits.table$V1, breaks=x_axis_breaks_nature_right_scale, col="blue", xlim=c(200,2300), ylim=c(0,18000), main="Histogram of Nature data peak --call-summits", xlab="length (bp)", ylab="Number of peaks")
histogram.right.scale.nocallsummits <-hist (Nature.no.callsummits.table$V1, breaks=x_axis_breaks_nature_right_scale, col="blue", xlim=c(200,2300), ylim=c(0,18000), main="Histogram of Nature data peak no --call-summits", xlab="length (bp)", ylab="Number of peaks")

##Mydata --call-summits
Mydata.callsummits.table <- read.delim ("Mydata_call_summits_length.txt", header=FALSE)
head (Mydata.callsummits.table)
class (Mydata.callsummits.table)
dim (Mydata.callsummits.table)
#Creates a numericvector to use as input for the hist break
x_axis_breaks_mydata <- seq(from = 200, to = 4200, by = 200) #numeric vector that starts at 200 and reaches 4200 by spacing the values every 200
par(mfrow = c(1,2))
par(mar = rep(4, 4))
par(lab = (c(21, 10, 3))) 
Mydata.histogram.callsummits <- hist (Mydata.callsummits.table$V1, breaks=x_axis_breaks_mydata, col="blue", xlim=c(200,4200), ylim=c(0,20000), main="Histogram of My data peak --call-summits", xlab="length (bp)", ylab="Number of peaks")
max(Mydata.callsummits.table$V1) #report maximum value in the data frame column V1
min(Mydata.callsummits.table$V1) #report minimum value in the data frame column V1
median (Mydata.callsummits.table$V1)

##Mydata without --call-summits
Mydata.no.callsummits.table <- read.delim ("Mydata_no_call_summits_length.txt", header=FALSE)
head (Mydata.no.callsummits.table)
class (Mydata.no.callsummits.table)
dim (Mydata.no.callsummits.table)
Mydata.histogram.no.callsummits <- hist (Mydata.callsummits.table$V1, breaks=x_axis_breaks_mydata, col="blue", xlim=c(200,4200), ylim=c(0,20000), main="Histogram of My data peak no --call-summits", xlab="length (bp)", ylab="Number of peaks")
max(Mydata.callsummits.table$V1) #report maximum value in the data frame column V1
min(Mydata.callsummits.table$V1) #report minimum value in the data frame column V1
median (Mydata.callsummits.table$V1)


#To generate the two histograms (side by side) of My data and Nature data in the same scale

#--call-summits:
par(mfrow = c(1,2))
par(mar = rep(4, 4)) #rep replicates the first element in this case 2 four times. It sets the margins of the plot so it fits.Otherwise I was getting an error
par(lab = (c(21, 10, 3)))
Mydata.histogram.callsummits <- hist (Mydata.callsummits.table$V1, breaks=x_axis_breaks_mydata, col="blue", xlim=c(200,4200), ylim=c(0,20000), main="Histogram of My data peak --call-summits", xlab="length (bp)", ylab="Number of peaks") 
histogram.changed.scale.callsummits <-hist (Nature.callsummits.table$V1, breaks=x_axis_breaks_nature_changed_scale, col="blue", xlim=c(200,4200), ylim=c(0,20000), main="Histogram of Nature data peak --call-summits changed scale", xlab="length (bp)", ylab="Number of peaks")

#no --call-summits:
#--call-summits:
par(mfrow = c(1,2))
par(mar = rep(4, 4)) #rep replicates the first element in this case 2 four times. It sets the margins of the plot so it fits.Otherwise I was getting an error
par(lab = (c(21, 10, 3)))
Mydata.histogram.no.callsummits <- hist (Mydata.callsummits.table$V1, breaks=x_axis_breaks_mydata, col="blue", xlim=c(200,4200), ylim=c(0,20000), main="Histogram of My data peak no --call-summits", xlab="length (bp)", ylab="Number of peaks")
histogram.changed.scale.nocallsummits <-hist (Nature.no.callsummits.table$V1, breaks=x_axis_breaks_nature_changed_scale, col="blue", xlim=c(200,4200), ylim=c(0,20000), main="Histogram of Nature data peak no --call-summits changed scale", xlab="length (bp)", ylab="Number of peaks")

