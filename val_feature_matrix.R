library(dplyr)
library(reshape2)
library(tidyr)
library(stringr)

setwd("/Users/oliver/RStudio")
structure <- read.delim("val_MFE.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE) # minimum free energy from ViennaRNA
nuc <- read.delim("val_nuc.count.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE) # Tm file 
score <- read.delim("val_score.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE) # sgRNA cutting efficiency score file

# Subset score file to include only unique ID and cutting efficiency score
score.df <- score[,c(1:2)]
colnames(score.df) <- c("sgRNAID", "cut.score")

# Create data frames for structure, gc, and temperature values
structure.df <- data.frame(sgRNAID = structure[,1], structure = structure[,2])
gc.df <- data.frame(sgRNAID = nuc[,1], gc = nuc[,3])
temp.df <- data.frame(sgRNAID = nuc[,1], temp = nuc[,2])

# Merge data frames
structure.temp <- left_join(structure.df, temp.df, by=c("sgRNAID"))
structure.temp.gc <- left_join(structure.temp, gc.df, by=c("sgRNAID"))
score.structure.temp.gc <- left_join(score.df, structure.temp.gc, by=c("sgRNAID"))

# Rename columns for clarity
colnames(score.structure.temp.gc) <- c("sgRNAID", "cut.score", "structure", "temp", "gc")


######## Add one-hot features
onehot.ind1 <- read.delim("val_ind1.txt", header=T, sep=",")
onehot.ind2 <- read.delim("val_ind2.txt", header=T, sep=",")
onehot.dep1 <- read.delim("val_dep1.txt", header=T, sep=",")
onehot.dep2 <- read.delim("val_dep2.txt", header=T, sep=",")
colnames(onehot.dep1)[1] <- "sgRNAID"
colnames(onehot.dep2)[1] <- "sgRNAID"

onehot.ind <- left_join(onehot.ind1, onehot.ind2, by="sgRNAID")
onehot.dep <- full_join(onehot.dep1[,1:ncol(onehot.dep1)-1], onehot.dep2[,1:ncol(onehot.dep2)-1], by="sgRNAID")

onehot <- full_join(onehot.ind, onehot.dep, by="sgRNAID")

raw <- left_join(score.structure.temp.gc, onehot, by=c("sgRNAID"))

write.table(raw, "val_raw.matrix.txt", quote=F, row.names=F, sep="\t")




######## Calculate QCT features
# Read the val_data.csv file
val_data <- read.csv("val_data.csv", header = TRUE, stringsAsFactors = FALSE)

# Filter for sequences that are exactly 20 characters long
val_data <- val_data %>% filter(nchar(nucleotide.sequence) == 20)

# Transform the nucleotide sequence into individual positions (p1 to p20)
data_transformed <- val_data %>%
  mutate(nucleotide.sequence = str_trim(nucleotide.sequence)) %>%
  mutate(nucleotide.sequence = strsplit(nucleotide.sequence, "")) %>%
  unnest_wider(nucleotide.sequence, names_sep = "")
colnames(data_transformed)[2:21] <- paste0("p", 1:20)

# Save the transformed data to val_sgRNA.sequence.txt
write.table(data_transformed, "val_sgRNA.sequence.txt", quote = FALSE, row.names = FALSE, sep = " ")



# Monomer QCT
# Read sgRNA sequence data
seq <- read.delim("val_sgRNA.sequence.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

tensor <- read.delim("HL.Bond.Monomer.txt", header=T, sep="\t", stringsAsFactors = F)

tensor[, 1] <- make.unique(tensor[, 1])
tensor.features <- tensor[,1]
rownames(tensor) <- tensor[,1]
tensor.df <- tensor[,2:ncol(tensor)]
tensor.t <- as.data.frame(t(tensor.df))
tensor.t$base <- names(tensor[,2:ncol(tensor)])

rownames(seq) <- seq[,1]
seq.melt <- melt(seq, id="sgRNAID")
colnames(seq.melt) <- c("sgRNAID", "position", "base")

seq.tensor <- left_join(seq.melt, tensor.t, by="base")
seq.tensor.melt <- melt(seq.tensor, id=c("sgRNAID", "position", "base"))
monomer <- dcast(seq.tensor.melt, sgRNAID ~ position + variable, value.var="value")
monomer <- monomer %>% select_if(~!all(is.na(.)))



# Basepair QCT
# Read sgRNA sequence data
seq <- read.delim("val_sgRNA.sequence.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

# Function to replace nucleotides with base pairs
replace_with_base_pairs <- function(nucleotide) {
  base_pairs <- nucleotide
  base_pairs[nucleotide == "A"] <- "AT"
  base_pairs[nucleotide == "T"] <- "TA"
  base_pairs[nucleotide == "G"] <- "GC"
  base_pairs[nucleotide == "C"] <- "CG"
  return(base_pairs)
}

# Apply the function from p1 to p20
seq[paste0("p", 1:20)] <- lapply(seq[paste0("p", 1:20)], replace_with_base_pairs)

# Print the updated seq dataframe to verify
print(head(seq))

# Save the updated seq dataframe
write.table(seq, "sgRNA_with_basepairs.txt", quote = FALSE, row.names = FALSE, sep = "\t")

tensor <- read.delim("HL.Bond.Basepair.txt", header=T, sep="\t", stringsAsFactors = F)
tensor.features <- tensor[,1]
rownames(tensor) <- tensor[,1]
tensor.df <- tensor[,2:ncol(tensor)]
tensor.t <- as.data.frame(t(tensor.df))
tensor.t$base <- names(tensor[,2:ncol(tensor)])

rownames(seq) <- seq[,1]
seq.melt <- melt(seq, id="sgRNAID")
colnames(seq.melt) <- c("sgRNAID", "position", "base")

seq.tensor <- left_join(seq.melt, tensor.t, by="base")
seq.tensor.melt <- melt(seq.tensor, id=c("sgRNAID", "position", "base"))
basepair <- dcast(seq.tensor.melt, sgRNAID ~ position + variable, value.var="value")
basepair <- basepair %>% select_if(~!all(is.na(.)))


# Dimer QCT
# Read sgRNA sequence data
seq <- read.delim("val_sgRNA.sequence.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

tensor <- read.delim("HL.Bond.Dimer.txt", header=T, sep="\t", stringsAsFactors = F)
seq.dimer <- seq %>% unite("p1", p1:p2, remove=F, sep= "") %>% unite("p2", p2:p3, remove=F, sep= "") %>% unite("p3", p3:p4, remove=F, sep= "") %>% unite("p4", p4:p5, remove=F, sep= "") %>% unite("p5", p5:p6, remove=F, sep= "") %>% unite("p6", p6:p7, remove=F, sep= "") %>% unite("p7", p7:p8, remove=F, sep= "") %>% unite("p8", p8:p9, remove=F, sep= "") %>% unite("p9", p9:p10, remove=F, sep= "") %>% unite("p10", p10:p11, remove=F, sep= "") %>% unite("p11", p11:p12, remove=F, sep= "") %>% unite("p12", p12:p13, remove=F, sep= "") %>% unite("p13", p13:p14, remove=F, sep= "") %>% unite("p14", p14:p15, remove=F, sep= "") %>% unite("p15", p15:p16, remove=F, sep= "") %>% unite("p16", p16:p17, remove=F, sep= "") %>% unite("p17", p17:p18, remove=F, sep= "") %>% unite("p18", p18:p19, remove=F, sep= "") %>% unite("p19", p19:p20, remove=T, sep= "")

tensor.features <- tensor[,1]
rownames(tensor) <- tensor[,1]
tensor.df <- tensor[,2:ncol(tensor)]
tensor.t <- as.data.frame(t(tensor.df))
tensor.t$base <- names(tensor[,2:ncol(tensor)])

rownames(seq.dimer) <- seq.dimer[,1]
seq.df <- seq.dimer[,1:20]
seq.melt <- melt(seq.df, id="sgRNAID")
colnames(seq.melt) <- c("sgRNAID", "position", "base")

seq.tensor <- left_join(seq.melt, tensor.t, by="base")
seq.tensor.melt <- melt(seq.tensor, id=c("sgRNAID", "position", "base"))
dimer <- dcast(seq.tensor.melt, sgRNAID ~ position + variable, value.var="value")
dimer <- dimer %>% select_if(~!all(is.na(.)))


# Trimer QCT
# Read sgRNA sequence data
seq <- read.delim("val_sgRNA.sequence.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

tensor <- read.delim("HL.Bond.Trimer.txt", header=T, sep="\t", stringsAsFactors = F)
seq.trimer <- seq %>% unite("p1", p1:p3, remove=F, sep= "") %>% unite("p2", p2:p4, remove=F, sep= "") %>% unite("p3", p3:p5, remove=F, sep= "") %>% unite("p4", p4:p6, remove=F, sep= "") %>% unite("p5", p5:p7, remove=F, sep= "") %>% unite("p6", p6:p8, remove=F, sep= "") %>% unite("p7", p7:p9, remove=F, sep= "") %>% unite("p8", p8:p10, remove=F, sep= "") %>% unite("p9", p9:p11, remove=F, sep= "") %>% unite("p10", p10:p12, remove=F, sep= "") %>% unite("p11", p11:p13, remove=F, sep= "") %>% unite("p12", p12:p14, remove=F, sep= "") %>% unite("p13", p13:p15, remove=F, sep= "") %>% unite("p14", p14:p16, remove=F, sep= "") %>% unite("p15", p15:p17, remove=F, sep= "") %>% unite("p16", p16:p18, remove=F, sep= "") %>% unite("p17", p17:p19, remove=F, sep= "") %>% unite("p18", p18:p20, remove=F, sep= "")

tensor[, 1] <- make.unique(tensor[, 1])
tensor.features <- tensor[,1]
rownames(tensor) <- tensor[,1]
tensor.df <- tensor[,2:ncol(tensor)]
tensor.t <- as.data.frame(t(tensor.df))
tensor.t$base <- names(tensor[,2:ncol(tensor)])

rownames(seq.trimer) <- seq.trimer[,1]
seq.df <- seq.trimer[,1:19]
seq.melt <- melt(seq.df, id="sgRNAID")
colnames(seq.melt) <- c("sgRNAID", "position", "base")

seq.tensor <- left_join(seq.melt, tensor.t, by="base")
seq.tensor.melt <- melt(seq.tensor, id=c("sgRNAID", "position", "base"))
trimer <- dcast(seq.tensor.melt, sgRNAID ~ position + variable, value.var="value")
trimer <- trimer %>% select_if(~!all(is.na(.)))



# Tetramer QCT
# Read sgRNA sequence data
seq <- read.delim("val_sgRNA.sequence.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

tensor <- read.delim("HL.Bond.Tetramer.txt", header=T, sep="\t", stringsAsFactors = F)
seq.tetramer <- seq %>% unite("p1", p1:p4, remove=F, sep= "") %>% unite("p2", p2:p5, remove=F, sep= "") %>% unite("p3", p3:p6, remove=F, sep= "") %>% unite("p4", p4:p7, remove=F, sep= "") %>% unite("p5", p5:p8, remove=F, sep= "") %>% unite("p6", p6:p9, remove=F, sep= "") %>% unite("p7", p7:p10, remove=F, sep= "") %>% unite("p8", p8:p11, remove=F, sep= "") %>% unite("p9", p9:p12, remove=F, sep= "") %>% unite("p10", p10:p13, remove=F, sep= "") %>% unite("p11", p11:p14, remove=F, sep= "") %>% unite("p12", p12:p15, remove=F, sep= "") %>% unite("p13", p13:p16, remove=F, sep= "") %>% unite("p14", p14:p17, remove=F, sep= "") %>% unite("p15", p15:p18, remove=F, sep= "") %>% unite("p16", p16:p19, remove=F, sep= "") %>% unite("p17", p17:p20, remove=F, sep= "") 

tensor[, 1] <- make.unique(tensor[, 1])
tensor.features <- tensor[,1]
rownames(tensor) <- tensor[,1]
tensor.df <- tensor[,2:ncol(tensor)]
tensor.t <- as.data.frame(t(tensor.df))
tensor.t$base <- names(tensor[,2:ncol(tensor)])

rownames(seq.tetramer) <- seq.tetramer[,1]
seq.df <- seq.tetramer[,1:18]
seq.melt <- melt(seq.df, id="sgRNAID")
colnames(seq.melt) <- c("sgRNAID", "position", "base")

seq.tensor <- left_join(seq.melt, tensor.t, by="base")
seq.tensor.melt <- melt(seq.tensor, id=c("sgRNAID", "position", "base"))
tetramer <- dcast(seq.tensor.melt, sgRNAID ~ position + variable, value.var="value")
tetramer <- tetramer %>% select_if(~!all(is.na(.)))


# Combine all QCT kmers by sgRNAID
combined_qct <- full_join(monomer, basepair, by = "sgRNAID") %>%
  full_join(dimer, by = "sgRNAID") %>%
  full_join(trimer, by = "sgRNAID") %>%
  full_join(tetramer, by = "sgRNAID")

# Clean column names: remove or replace problematic characters
clean_colnames <- function(df) {
  colnames(df) <- iconv(colnames(df), to = "ASCII//TRANSLIT", sub = "")
  colnames(df) <- gsub("[^[:alnum:]_]", "_", colnames(df))
  colnames(df) <- make.names(colnames(df), unique = TRUE)
  return(df)
}

combined_qct <- clean_colnames(combined_qct)

# Save the combined dataframe as CSV
write.csv(combined_qct, "val_quantum.melt.csv", row.names = FALSE)
write.table(combined_qct, "val_quantum.melt.txt", quote = F, row.names = F, sep = " ")


### combine raw matrix and quantum matrix to generate final validation feature matrix
df <- read.delim("val_raw.matrix.txt", header=T, sep="\t", stringsAsFactors = F)

# quantum chemical tensors
tensor <- read.csv("val_quantum.melt.csv", header = TRUE, stringsAsFactors = FALSE)
tensor[is.na(tensor)] <- 0

raw.tensor <- inner_join(raw, tensor, by="sgRNAID")

write.table(raw.tensor, "val_final_quantum.txt", quote=F, row.names=F, sep="\t")
write.csv(raw.tensor, "val_final_quantum.csv", row.names = FALSE)






