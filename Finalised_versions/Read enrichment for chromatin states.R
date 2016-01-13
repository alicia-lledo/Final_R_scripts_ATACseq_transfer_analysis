getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data analysis/QC_plots/pipeline")
getwd ()
library(ggplot2)
library(gplots)


#1.READ THE FILES
  #col1 is chr
  #col2 length of the fragment
  #col 3 annotation of the fragment using CD14 segments Rodmap file
read_fragments_annotated <- read.delim("./lesional_bowtie/psoriasis_bowtie_annotated_fragments.bed", header=FALSE)



#2. RENAME LABELS WITH CHROMATIN STATES
#Re-name the col 3 replacing numeric labels with states
  #For each row where col3 has a particular value replace it by the given name

read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 1 )]<- "TssA"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 2 )]<- "TssAFlnk"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 3 )]<- "TssAFlnkU"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 4 )]<- "TssFlnkD"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 5)]<- "Tx"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 6 )]<- "TxWk"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 7)]<- "EnhG1"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 8)]<- "EnhG2"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 9)]<- "EnhA1"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 10 )]<- "EnhA2"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 11)]<- "EnhWk"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 12 )]<- "ZNF/Rpts"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 13 )]<- "Het"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 14 )]<- "TssBiv"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 15 )]<- "EnhBiv"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 16 )]<- "ReprPC"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 17 )]<- "ReprPC (weak)"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 18 )]<- "Quies"




#3. DENSITY PLOT OF FRAGMENT SIZE FOR EACH CHROMATIN STATE

#OPTION A.To create a density plot I could build a frequency table using the function table()
  #table(read_fragments_annotated$V3) #to create a contingency table counting the number of events at each of the factors levels fed into the function

#OPTION B.To create a histogram/density plot of the states of chromatin at each fragment size
  ##ggplot geonmdensities is able to calculate frequencies from a table with single events, without need of generating first a frequency table.
ggplot(read_fragments_annotated, aes(x=V2)) + geom_density(aes(group=V3, col=V3)) + xlim(0, 800)+ 
  xlab("Fragment size (bp)") + ylab("Density")
      ###I need to modify the legen title




#4. TO FED THE DATA INTO A MATRIX

# create a matrix with a row for all possible fragment sizes, and a column for each chromatin category and filled with 0
my.matrix<-matrix (0, nrow=max(read_fragments_annotated$V2), ncol=18) #the number of sizes can be found by doing range (read_fragment_annotation$V2)
                  #data=0 to fill it with 0
                  #nrow and ncol to indicate the size of the matrix
range (read_fragments_annotated$V2) #it tells the numbers between which the fragment sizes range (19 and 3516)

#name the columns of the matrix
states <- c("TssA", "TssAFlnk","TssAFlnkU", "TssFlnkD", "Tx", "TxWk", "EnhG1", "EnhG2", 
          "EnhA1", "EnhA2", "EnhWk", "ZNF/Rpts", "Het", "TssBiv", "EnhBiv", "ReprPC", "ReprPC (weak)", "Quies")
colnames (my.matrix) <- states

#5.USE THE MATRIX TO BUILD A FREQUENCY TABLE TO LATER PLOT THE HEATMAP

#Matrix of frequencies of fragment size (row) for each chr state (columns)  and filled with frequencies (x,y and z variables)
  # Use the function table () to calculate frequencies


# For each category, count up how often each fragment length is seen by creating different table objects (which is like a vector witha factor that labelsthe contained values, in this case fragment sizes)
  #subset () Return subsets of vectors, matrices or data frames which meet conditions
  #$V2 indicates to use V2 as the factor level to report the frequency counts

table.TssA <- table(subset(read_fragments_annotated, V3 == "TssA")$V2) 
table.TssAFlnk <- table(subset(read_fragments_annotated, V3 == "TssAFlnk")$V2)
table.TssFlnkU <- table(subset(read_fragments_annotated, V3 == "TssAFlnkU")$V2)
table.TssFlnkD<- table(subset(read_fragments_annotated, V3 =="TssFlnkD")$V2)
table.Tx<- table(subset(read_fragments_annotated, V3 == "Tx")$V2)
table.TxWk<- table(subset(read_fragments_annotated, V3 =="TxWk")$V2)
table.EnhG1<- table(subset(read_fragments_annotated, V3 == "EnhG1")$V2)
table.EnhG2<- table(subset(read_fragments_annotated, V3 == "EnhG2")$V2)
table.EnhA1<- table(subset(read_fragments_annotated, V3 == "EnhA1")$V2)
table.EnhA2<- table(subset(read_fragments_annotated, V3 == "EnhA2")$V2)
table.EnhWk<- table(subset(read_fragments_annotated, V3 == "EnhWk")$V2)
table.ZNF_Rpts<- table(subset(read_fragments_annotated, V3 =="ZNF/Rpts")$V2)
table.Het <- table(subset(read_fragments_annotated, V3 == "Het")$V2)
table.TssBiv <- table(subset(read_fragments_annotated, V3 == "TssBiv")$V2)
table.EnhBiv<- table(subset(read_fragments_annotated, V3 =="EnhBiv" )$V2)
table.ReprPC<- table(subset(read_fragments_annotated, V3 == "ReprPC")$V2)
table.ReprPC_weak <- table(subset(read_fragments_annotated, V3 == "ReprPC (weak)")$V2)
table.Quies<- table(subset(read_fragments_annotated, V3 =="Quies")$V2)






# Insert the frequency counts into the matrix in the appropriate rows
  #names () extract the name of the object, in this case each of the levels from the table object to which the frequency count has been performed,which is the different fragment lengths
  #It looks at the row position of the matrix corresponsing to the the number stored in the level of the table (which is the fragment length)
  #For each chr state the count is filled in a different column
my.matrix[as.numeric((names(table.TssA))), 1] <- as.vector(table.TssA)                                                                  #fill in rows with the frquency counts from the tables
my.matrix[as.numeric((names(table.TssAFlnk))), 2] <- as.vector(table.TssAFlnk) 
my.matrix[as.numeric((names(table.TssFlnkU))), 3] <- as.vector(table.TssFlnkU)
my.matrix[as.numeric((names(table.TssFlnkD))), 4] <- as.vector(table.TssFlnkD)
my.matrix[as.numeric((names(table.Tx))), 5] <- as.vector(table.Tx)
my.matrix[as.numeric((names(table.TxWk))), 6] <- as.vector(table.TxWk)
my.matrix[as.numeric((names(table.EnhG1))), 7] <- as.vector(table.EnhG1)
my.matrix[as.numeric((names(table.EnhG2))), 8] <- as.vector(table.EnhG2)
my.matrix[as.numeric((names(table.EnhA1))), 9] <- as.vector(table.EnhA1)
my.matrix[as.numeric((names(table.EnhA2))), 10] <- as.vector(table.EnhA2)
my.matrix[as.numeric((names(table.EnhWk))), 11] <- as.vector(table.EnhWk)
my.matrix[as.numeric((names(table.ZNF_Rpts))), 12] <- as.vector(table.ZNF_Rpts)
my.matrix[as.numeric((names(table.Het))), 13] <- as.vector(table.Het)
my.matrix[as.numeric((names(table.TssBiv))), 14] <- as.vector(table.TssBiv)
my.matrix[as.numeric((names(table.EnhBiv))), 15] <- as.vector(table.EnhBiv)
my.matrix[as.numeric((names(table.ReprPC))), 16] <- as.vector(table.ReprPC)
my.matrix[as.numeric((names(table.ReprPC_weak))), 17] <- as.vector(table.ReprPC_weak)
my.matrix[as.numeric((names(table.Quies))), 18] <- as.vector(table.Quies)



#6.PLOT THE FREQUENCY DATA INTO A HEATMAP (PREIOR TO NORMALISATION)

#t() transposes col and rows. Now rows are chromatin states and col the fragment sizes
#colv and rowv  when NA avoids doing clustering of the rows and col

transposed.matrix<-t(my.matrix) #to have the chr states in the y axis

#To change colours
newcolours<-colorRampPalette(c("white","blue"))(256) #function colorRampPalette () pass it the name of colours following R conventions, and in brackets the number of tones in between
labels<- c(rep("",49),"50",rep("",49),"100", rep("",49),"150", rep("",49),"200",rep("",49),"250",
           rep("",49),"300",rep("",49),"350",rep("",49),"400", rep("",49),"450", rep("",49),"500", rep("",49),"550",
           rep("",49),"600", rep("",49),"650", rep("",49),"700", rep("",49),"750", rep("",49),"800", rep("",49),"850")


heatmap.2(transposed.matrix[, 1:850], Colv = FALSE, Rowv = FALSE, labCol = labels, col=newcolours, scale="none", trace="none")

#After normalisation I will have positive and negative values it will assign blue<0 white+0 and red>0
colorRampPalette(c("blue","white","red"))(256)




#NORMALISATION:

#1.Normalise each value of each category to the total number of reads of the category--> to make chr states comparable
  #e.g Count of each length annotates as TssA/ total number of reads across all lengths annotated as TssA

  #it calculates total num of fragments for each chr state (2 indicates calculating totals by col)
    #apply() iterates the same function over a data frame ; e.g here the sum()
totals <- apply(my.matrix, 2, sum)
class(totals)
str(totals)


  #Create a new matrix which will contain each of the counts for each fragment length and chr state divided by the total number of counts for that chr state
normalised.matrix.for.chr.states <- matrix (NA, ncol=18, nrow=max(read_fragments_annotated$V2))
normalised.matrix.for.fragment.length <-matrix (NA, ncol=18, nrow=max(read_fragments_annotated$V2))
  
#Perform the division by the totals for each of the values of the original matrix using a loop:
    #it uses the col variable to move across the different chr states columns and divide the original value by the total number of fragments in that chr state stored in the position col of the vector totals
for (col in 1:18){
  normalised.matrix.for.chr.states[,col] <-my.matrix[,col]/totals[col]
}

#2.Normalise to the total number of reads for a particular category and chr state--> to make different fragment length abundance comparable;to balance the different abundance
  #e.g divide each previously normalised count of a particular fragment length by the total number of reads of that length mapping to all the chr states
  
#Calculate total number of reads of a particular length in all the chr states using the function row.sum () which will add all the values for each row and store it in a vector which each sum value in the position of the fragment length

all.categories.sum <-rowSums (my.matrix[])
all.categories.sum[which(all.categories.sum==0)] <- NA



#Divide each value of the matrix by the total number of fragments of a length across all categories
  #use the following loop
for (row in 1:max(read_fragments_annotated$V2)){
  normalised.matrix.for.fragment.length[row,] <- normalised.matrix.for.chr.states[row,]/all.categories.sum[row]
}


#Divide by total number of annotated reads
total.annotated.reads <- nrow(read_fragments_annotated)
normalised.matrix.for.fragment.length[] <- normalised.matrix.for.fragment.length[]/total.annotated.reads



#log10 transform
log10.normalised.matrix<-log10(normalised.matrix.for.fragment.length)
log10.normalised.matrix[which(log10.normalised.matrix == "-Inf")] <- NA #to avoid the 0 to interfere and mess up the plotting
colnames (log10.normalised.matrix)<- states #to name the col of the final matrix and get the chr states

#plot the heatmap with these values (transpose the matrix to have chr states in Y and frag length in x):
newcolours<-colorRampPalette(c("blue", "white","red"))(100)
#breaks = seq(0,800,by=5) #to smooth the x axis
transposed.log10.normalised.matrix <- t(log10.normalised.matrix)
heatmap.2(transposed.log10.normalised.matrix[,50:550], labCol = labels, col=newcolours, scale = "row", 
          trace="none", Colv = FALSE, Rowv = FALSE, dendrogram = "none", na.rm=TRUE) #remove the rows where the NA are
  #heatmap2 () function scales colouring based on rows or colours assuming that tey follow a normal distribution with mean 0 and std 1







