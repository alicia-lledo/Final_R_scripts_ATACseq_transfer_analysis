

#The script takes as input the file with all the fragments annotated for the 18 chromatin states and splits them in different tables


#1.READ THE FILES
#col1 is chr
#col2 length of the fragment
#col 3 annotation of the fragment using CD14 segments Rodmap file
#read_fragments_annotated <- read.delim("/well/jknight/ATACseq/ATACseq_001/pipeline_blood_bwa/blood_bwa_annotated_fragments.bed", header=FALSE)

read_fragments_annotated <- read.delim("/well/jknight/ATACseq/ATACseq_less_stringent_pipeline/new_pipeline/lesional_skin.atac_analysis/lesional_annotated_fragments.bed", header=FALSE)


#1. RENAME LABELS WITH CHROMATIN STATES
#Re-name the col 3 replacing numeric labels with states
#For each row where col3 has a particular value replace it by the given name

read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 1 )]<- "TssA"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 2 )]<- "TssFlnk"
read_fragments_annotated$V3[which(read_fragments_annotated$V3 == 3 )]<- "TssFlnkU"
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




#2.SPLIT THE ANNOTATION STATES IN DIFFERENT TABLES
TssA <- subset(read_fragments_annotated, V3 == "TssA") 
TssFlnk <- subset(read_fragments_annotated, V3 == "TssFlnk")
TssFlnkU <- subset(read_fragments_annotated, V3 == "TssFlnkU")
TssFlnkD<- subset(read_fragments_annotated, V3 =="TssFlnkD")
Tx<- subset(read_fragments_annotated, V3 == "Tx")
TxWk<- subset(read_fragments_annotated, V3 =="TxWk")
EnhG1<- subset(read_fragments_annotated, V3 == "EnhG1")
EnhG2<- subset(read_fragments_annotated, V3 == "EnhG2")
EnhA1<- subset(read_fragments_annotated, V3 == "EnhA1")
EnhA2<- subset(read_fragments_annotated, V3 == "EnhA2")
EnhWk<- subset(read_fragments_annotated, V3 == "EnhWk")
ZNF_Rpts<- subset(read_fragments_annotated, V3 =="ZNF/Rpts")
Het <- subset(read_fragments_annotated, V3 == "Het")
TssBiv <- subset(read_fragments_annotated, V3 == "TssBiv")
EnhBiv<- subset(read_fragments_annotated, V3 =="EnhBiv" )
ReprPC<- subset(read_fragments_annotated, V3 == "ReprPC")
ReprPC_weak <- subset(read_fragments_annotated, V3 == "ReprPC (weak)")
Quies<- subset(read_fragments_annotated, V3 =="Quies")

write.table(TssA, file="TssA.txt", sep="\t", row.names = FALSE)
write.table(TssFlnk, file="TssFlnk.txt", sep="\t", row.names = FALSE)
write.table(TssFlnkU, file="TssFlnkU.txt", sep="\t", row.names = FALSE)
write.table(TssFlnkD, file="TssFlnkD.txt", sep="\t", row.names = FALSE)
write.table(Tx, file="Tx.txt", sep="\t", row.names = FALSE)
write.table(TxWk, file="TxWk.txt", sep="\t", row.names = FALSE)
write.table(EnhG1, file="EnhG1.txt", sep="\t", row.names = FALSE)
write.table(EnhG2, file="EnhG2.txt", sep="\t", row.names = FALSE)
write.table(EnhA1, file="EnhA1.txt", sep="\t", row.names = FALSE)
write.table(EnhA2, file="EnhA2.txt", sep="\t", row.names = FALSE)
write.table(EnhWk, file="EnhWk.txt", sep="\t", row.names = FALSE)
write.table(ZNF_Rpts, file="ZNF_Rpts.txt", sep="\t", row.names = FALSE)
write.table(Het, file="Het.txt", sep="\t", row.names = FALSE)
write.table(TssBiv, file="TssBiv.txt", sep="\t", row.names = FALSE)
write.table(EnhBiv, file="EnhBiv.txt", sep="\t", row.names = FALSE)
write.table(ReprPC, file="ReprPC.txt", sep="\t", row.names = FALSE)
write.table(ReprPC_weak, file="ReprPC_weak.txt", sep="\t", row.names = FALSE)
write.table(Quies, file="Quies.txt", sep="\t", row.names = FALSE)


#3. PULL THE FRAGMENT SIZE COLUMN AND WRITE IT IN A DIFFERENT FILE

write.table(read_fragments_annotated$V2, file="../gw.fragments.txt", row.names = FALSE)
