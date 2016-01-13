getwd ()
setwd ("U:/PhD_U/Year_1/ATAC seq/Data Analysis/Peak_calling")
getwd ()

#To create a txt file with the content of UCSC_Nature_CD14_alignment_without_duplicas_MAPQ_30_unique_additional_flags_all_no_mitochondria_blacklist_filtered_PC_peaks.narrowPeak.txt and the 4th column which contains the name of each peak changed to "Nature"

#Read table with table.delim which assumens "\t" separated file by default and indicate there is no header
Nature_data <- read.delim ("UCSC_Nature_CD14_alignment_without_duplicas_MAPQ_30_unique_additional_flags_all_no_mitochondria_blacklist_filtered_PC_peaks.narrowPeak.txt", header=FALSE)
#Print the first few rows and col from the table
head (Nature_data)
#Print the first few values  of all the raws for column 4
head (Nature_data [,4])
#Replace the values of all the rows of col 4 by the string "Nature"
Nature_data [,4] <- "Nature"
#Check the replacement has been done
head (Nature_data)
#To write the changes and create a file containing the modified BED table:
write.table(Nature_data, file = "/PhD_U/Year_1/ATAC seq/Data Analysis/Peak_calling/GAT_Nature_CD14_alignment_without_duplicas_MAPQ_30_unique_additional_flags_all_no_mitochondria_blacklist_filtered_PC_peaks.narrowPeak.txt", append = FALSE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "\t")
Nature.testing.track_file <-read.table("GAT_Nature_CD14_alignment_without_duplicas_MAPQ_30_unique_additional_flags_all_no_mitochondria_blacklist_filtered_PC_peaks.narrowPeak.txt", header=TRUE, sep="\t")
head (Nature.testing.track_file)


#To create a txt file with the content of UCSC_my_data_alignment_without_duplicas_MAPQ_30_unique_additional_flags_all_no_mitochondria_call_summits_blacklist_filtered_PC_peaks.narrowPeak.txt and the 4th column which contains the name of each peak changed to "Nature"

#Read table with table.delim which assumens "\t" separated file by default and indicate there is no header
My_data <- read.delim ("UCSC_my_data_alignment_without_duplicas_MAPQ_30_unique_additional_flags_all_no_mitochondria_call_summits_blacklist_filtered_PC_peaks.narrowPeak.txt", header=FALSE)
#Print the first few rows and col from the table
head (My_data)
#Print the first few values  of all the raws for column 4
head (My_data [,4])
#Replace the values of all the rows of col 4 by the string "Nature"
My_data [,4] <- "Mydata"
#Check the replacement has been done
head (My_data)
#To write the changes and create a file containing the modified BED table:
write.table(My_data, file = "GAT_my_data_alignment_without_duplicas_MAPQ_30_unique_additional_flags_all_no_mitochondria_call_summits_blacklist_filtered_PC_peaks.narrowPeak.txt", append = FALSE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "\t")
Mydata.testing.track_file <-read.table("GAT_my_data_alignment_without_duplicas_MAPQ_30_unique_additional_flags_all_no_mitochondria_call_summits_blacklist_filtered_PC_peaks.narrowPeak.txt", header=TRUE, sep="\t")
head (Mydata.testing.track_file)


#To create a txt file with the content of L5_TGGTCA_L005_001.trim.st.all.rmdup_Nature_CD14_normalised_PC_blacklist_filtered.narrowPeak.txt and the 4th column which contains the name of each peak changed to "Nature_normalised"

#Read table with table.delim which assumens "\t" separated file by default and indicate there is no header
Nature_normalised_data <- read.delim ("L5_TGGTCA_L005_001.trim.st.all.rmdup_Nature_CD14_normalised_PC_blacklist_filtered.narrowPeak.txt", header=FALSE)
#Print the first few rows and col from the table
head (Nature_normalised_data)
#Print the first few values  of all the raws for column 4
head (Nature_normalised_data [,4])
#Replace the values of all the rows of col 4 by the string "Nature_normalised"
Nature_normalised_data [,4] <- "Nature_normalised"
#Check the replacement has been done
head (Nature_normalised_data)
#To write the changes and create a file containing the modified BED table:
write.table(Nature_normalised_data, file = "/PhD_U/Year_1/ATAC seq/Data Analysis/Peak_calling/GAT_Nature_normalised_CD14_alignment_without_duplicas_MAPQ_30_unique_additional_flags_no_mitochondria_blacklist_filtered_PC_peaks.narrowPeak.txt", append = FALSE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "\t")
Nature.normalised.testing.track_file <-read.table("GAT_Nature_normalised_CD14_alignment_without_duplicas_MAPQ_30_unique_additional_flags_no_mitochondria_blacklist_filtered_PC_peaks.narrowPeak.txt", header=TRUE, sep="\t")
head (Nature.normalised.testing.track_file)
dim (Nature.normalised.testing.track_file)


#Similarly to before to create a txt file with the content of hg19.bed.chrom.sizes.UCSC.txt and an additional 4th col to make the required bed file to feed into GAT. 
#quote=FALSE avoids characters to be written in the output in quotation marks and col.name and row.name FALSE avoids inclusing the names of rows and col as one additional row and col in the output
chromosomes <-read.delim ("hg19.bed.chrom.sizes.UCSC.txt", header=FALSE)#To create an R object containing the BED format table
chromosomes [,4]<- "ws" #To change the value of the 4th column for all the rows 
head (chromosomes)#To see first few lines of the table
write.table(chromosomes, file = "/PhD_U/Year_1/workspace_table.txt", append = FALSE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "\t")
testing.workspace <-read.table("workspace_table.txt", header=TRUE, sep="\t")
head (testing.workspace)
