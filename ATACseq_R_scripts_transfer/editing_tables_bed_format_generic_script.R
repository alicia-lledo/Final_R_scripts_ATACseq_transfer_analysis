getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/Peak_annotation")
getwd ()

#To create a txt file with the content replacing 4th column withan ID name

#Read table with table.delim which assumens "\t" separated file by default and indicate there is no header
table.transformation <- read.delim ("ENCODE_DUKE_DHS_Pk_bed_file.txt", header=FALSE)
#Print the first few rows and col from the table
head (table.transformation)
#Print the first few values  of all the raws for column 4
head (table.transformation [,4])
#Replace the values of all the rows of col 4 by an string ""
table.transformation [,4] <- "ENCODE_DUKE_CD14_DHS_Pk"
#Check the replacement has been done
head (table.transformation)
#To write the changes and create a file containing the modified BED table:
write.table(table.transformation, file = "/PhD_U/Year_1/ATAC seq/Data Analysis/Peak_annotation/GAT_ENCODE_DUKE_DHS_Pk_bed_file.txt", append = FALSE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "\t")
new.table <-read.table("GAT_ENCODE_DUKE_DHS_Pk_bed_file.txt", header=FALSE, sep="\t")
head (new.table)

