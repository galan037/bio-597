### Install packages ONLY IF NOT ALREADY INSTALLED
# Install DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# Intall dplyr and RColorBrewer
install.packages("tidyverse")
install.packages("RColorBrewer")


### Load in packages at the start of every R session
# Load DESeq2 package
library(DESeq2)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)

### Set working directory!!!
setwd("~/Desktop/2024_596-15_VFDB_analyses/")


###############################################
# Combine Salmon VFDB Outputs into single table
###############################################

# Get a list of all .sf files in the directory
file_list <- list.files(pattern="\\.sf$")

# Use purrr::map() to read all files into a list of dataframes
data_list <- map(file_list, ~read.delim(.x, header = TRUE) %>% select(Name, TPM))

# Add names to data_list files
names(data_list) <- file_list

# Use purrr::reduce() to join dataframes by shared variable Name
joined_data <- reduce(data_list, full_join, by = "Name")

# Rename TPM columns to "TPM.N"
tpm_names <- str_extract(names(data_list), "^[^_]+_[^_]+")
names(joined_data)[grep("TPM", names(joined_data))] <- tpm_names

# Make column 1 row names and remove column 1
joined_data <- joined_data %>% remove_rownames %>% column_to_rownames(var="Name")

# Remove rows that have all zeros
joined_data_no0 <- joined_data[rowSums(joined_data[])>0,]

# separate into rhizosphere and soil
joined_data_no0R <- joined_data_no0[, grepl("_R$", names(joined_data_no0))]
joined_data_no0S <- joined_data_no0[, grepl("_S$", names(joined_data_no0))]

# Read in metadata
VFDB <- read.delim("VFDB_ref.fas", sep = "\t", row.names = 1, header = T,stringsAsFactors=FALSE)

# merge metadata and quant values
count_data_R <- merge(VFDB, joined_data_no0R, by = 'row.names', type = "inner", match = "all")
count_data_S <- merge(VFDB, joined_data_no0S, by = 'row.names', type = "inner", match = "all")

colnames(count_data_R)[1:2] <- c("Virulence_factor","Virulence_factor_description")
colnames(count_data_S)[1:2] <- c("Virulence_factor","Virulence_factor_description")

# write csv file with final output
write.table(count_data_R, "VFDB Quants 596-15 Rhizo Samples.tsv", row.names = FALSE, quote = F, sep="\t")
write.table(count_data_S, "VFDB Quants 596-15 Soil Samples.tsv", row.names = FALSE, quote = F, sep="\t")

###############################################
# DESeq2
###############################################
## Prepare abundance table

# Load count data matrix with VFDBs as rows and samples as columns
count_data_R <- read.delim("VFDB Quants 596-15 Rhizo Samples.tsv", sep = "\t", header = T,row.names = 1, stringsAsFactors=FALSE)
count_data_S <- read.delim("VFDB Quants 596-15 Soil Samples.tsv", sep = "\t", header = T,row.names = 1, stringsAsFactors=FALSE)

## Prepare abundance table
# remove remaining text column
count_mat_R = round(na.omit(count_data_R[,!(names(count_data_R) %in% "Virulence_factor_description")]))
count_mat_S = round(na.omit(count_data_S[,!(names(count_data_S) %in% "Virulence_factor_description")]))


######## THE FOLLOWING NEEDS TO BE CHANGED BASED ON SAMPLE NAMES AND ORDER OF COLUMNS!!!!!!!!
## Create a table with sample metadata
# Create a vector with the soil type names

# Function to remove the second character from each string in the list
remove_second_char <- function(x) {
  paste0(substr(x, 1, 1), substr(x, 3, nchar(x)))
}

remove_R <- function(string_list) {
  # Remove "_R" from each string in the list
  new_string_list <- lapply(string_list, function(x) sub("_R$", "", x))
  
  # Return the modified list
  return(new_string_list)
}

remove_S <- function(string_list) {
  # Remove "_R" from each string in the list
  new_string_list <- lapply(string_list, function(x) sub("_S$", "", x))
  
  # Return the modified list
  return(new_string_list)
}

remove_number <- function(string_list) {
  # Remove the number part from each string in the list
  new_string_list <- gsub("\\d", "", string_list)
  
  return(new_string_list)
}
# Create separate TPM folders for rhizo and soil
tpm_names_R <- tpm_names[grep("_R$", tpm_names)]
tpm_names_S <- tpm_names[grep("_S$", tpm_names)]


# Apply the function to each element of the list
soil_types_R <- unlist(lapply(tpm_names_R, remove_R))
soil_types_S <- unlist(lapply(tpm_names_S, remove_S))

soil_types_R <- unlist(lapply(soil_types_R, remove_number))
soil_types_S <- unlist(lapply(soil_types_S, remove_number))


# Create a data frame with the row names and the soil_type column
sample_metadata_R <- data.frame(Soil_type = soil_types_R, row.names = tpm_names_R)
sample_metadata_S <- data.frame(Soil_type = soil_types_S, row.names = tpm_names_S)

# # Change
# colnames(sample_metadata_R) <- c("Site")
# # Extract last letters from tpm_names
# soil_type <- as.factor(substr(tpm_names_R, nchar(tpm_names_R), nchar(tpm_names_R)))
# 
# # Add new column Soil_Type to the data frame
# sample_metadata_R$Soil_type <- soil_type


######################GOOD TO THIS POINT#####################

# Create DESeq2 dataset object
dds_R <- DESeqDataSetFromMatrix(countData = count_mat_R, colData = sample_metadata_R, design = ~ Soil_type)
dds_S <- DESeqDataSetFromMatrix(countData = count_mat_S, colData = sample_metadata_S, design = ~ Soil_type)

# Perform differential abundance analysis
dds_R <- DESeq(dds_R)
dds_S <- DESeq(dds_S)


# Get differential expression results for all virulence factors
# Can change this to compare 2 specific soil type results, ie commented line below
res_rhizo <- results(dds_R)
res_soil <- results(dds_S)

# res <- results(dds, contrast = c("Soil_type", "Chaparral", "Lawn"))***


# Filter for virulence factor genes only
VF_res_R <- res_rhizo[row.names(res_rhizo) %in% row.names(count_data_R), ]
VF_res_S <- res_soil[row.names(res_soil) %in% row.names(count_data_S), ]

# Filter for significant differentially abundant VF genes (using an adjusted p-value cutoff of 0.05)
sig_VF_rhizo <- VF_res_R[which(VF_res_R$padj < 0.01), ]
sig_VF_soil <- VF_res_S[which(VF_res_S$padj < 0.01), ]

# Sort by log2-fold change (from most positive to most negative)
sig_VF_sorted_rhizo <- sig_VF_rhizo[order(-sig_VF_rhizo$log2FoldChange), ]
sig_VF_sorted_soil <- sig_VF_soil[order(-sig_VF_soil$log2FoldChange), ]

# View top 10 differentially abundant VF genes
head(sig_VF_sorted_rhizo, 10)
head(sig_VF_sorted_soil, 10)

# convert results to data frame
de_genes_rhizo <- as.data.frame(sig_VF_sorted_rhizo)
de_genes_soil <- as.data.frame(sig_VF_sorted_soil)

write.table(de_genes_rhizo, "VFDB Differentially Expressed Genes and Stats 596-15 Rhizo Samples.tsv", row.names = TRUE, quote = F, sep="\t")
write.table(de_genes_soil, "VFDB Differentially Expressed Genes and Stats 596-15 Soil Samples.tsv", row.names = TRUE, quote = F, sep="\t")

### subset results 
de_row_names_R <- row.names(sig_VF_sorted_rhizo)
de_row_names_S <- row.names(sig_VF_sorted_soil)


# Subset count_data based on differentially expressed rows
count_data_de_R <- count_data_R[de_row_names_R, ]
count_data_de_S <- count_data_S[de_row_names_S, ]

write.table(count_data_de_R, "VFDB Differentially Expressed 596-15 Rhizo Samples.tsv", row.names = TRUE, quote = F, sep="\t")
write.table(count_data_de_S, "VFDB Differentially Expressed 596-15 Soil Samples.tsv", row.names = TRUE, quote = F, sep="\t")


#### Cluster data by rows

# remove remaining text column
count_data_demat_R = as.matrix(round(na.omit(count_data_de_R[,!(names(count_data_de_R) %in% "Virulence_factor_description")])))
count_data_demat_S = as.matrix(round(na.omit(count_data_de_S[,!(names(count_data_de_S) %in% "Virulence_factor_description")])))

# calculate distance between rows using Euclidean distance, then perform heirarchical clustering
hc <- hclust(dist(count_data_demat_R), method = "ward.D2")
hc <- hclust(dist(count_data_demat_S), method = "ward.D2")


# obtain vector of cluster assignments for each row
row_clusters <- cutree(hc, k = 5)


# order rows of matrix based on cluster assignments
count_data_ordered_R <- count_data_demat_R[order(row_clusters), ]
count_data_ordered_S <- count_data_demat_S[order(row_clusters), ]

write.table(count_data_ordered_R, "VFDB Differentially Expressed Clustered 596-15 Rhizo Samples.tsv", row.names = TRUE, quote = F, sep="\t")
write.table(count_data_ordered_S, "VFDB Differentially Expressed Clustered 596-15 Soil Samples.tsv", row.names = TRUE, quote = F, sep="\t")

# Check the structure of your data frames to identify the common column(s)
str(count_data_ordered_R)
str(VFDB)

# Convert row names to a column in both datasets
VFDB$common_column <- rownames(VFDB)

count_data_ordered_S <- data.frame(common_column = rownames(count_data_ordered_S), count_data_ordered_S, row.names = NULL)
count_data_ordered_R <- data.frame(common_column = rownames(count_data_ordered_R), count_data_ordered_R, row.names = NULL)

# Merge datasets based on the common column
merged_data_R <- merge(VFDB, count_data_ordered_R, by = "common_column", all = TRUE)
merged_data_S <- merge(VFDB, count_data_ordered_S, by = "common_column", all = TRUE)

merged_data_R <- na.omit(merged_data_R)
merged_data_S <- na.omit(merged_data_S)

file_path_R <- "Rhizosphere"
file_path_S <- "Soil"


# Export merged data as a TSV file
write.table(merged_data, file = file_path_R, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(merged_data_S, file = file_path_S, sep = "\t", quote = FALSE, row.names = FALSE)


Rhizo_VFDB <- read.delim("VFDB_analyses_full.txt", sep = "\t", header = T,row.names = 1, stringsAsFactors=FALSE)
colnames(Rhizo_VFDB)[19] <- c("Virulence_factor_description")

# Example data frame
Rhizo_VFDB_sorted <- data.frame(Lowland = c(1:215),
                 Riperian = c(1:215),
                 Sedge = c(1:215)),


# Sort by column A, then by column B, then by column C
Rhizo_VFDB_final <- Rhizo_VFDB[order(Rhizo_VFDB_sorted$Lowland, Rhizo_VFDB_sorted$Riperian, Rhizo_VFDB_sorted$Sedge), ]

write.table(Rhizo_VFDB_final, "Rhizo_VFDB_final_ordered.tsv", row.names = TRUE, quote = F, sep="\t")
### Visualize with ggplot2 heatmap

# Convert count_data_demat matrix to a data frame with row and column indices as separate columns
count_data_df_R <- as.data.frame.table(count_data_ordered_R)
colnames(count_data_df_R) <- c("row", "column", "value")

count_data_df_S <- as.data.frame.table(count_data_ordered_S)
colnames(count_data_df_S) <- c("row", "column", "value")

# Define the color scale for the heatmap
colors <- colorRampPalette(brewer.pal(8, "BuPu"))(25)

# Create the heatmap using ggplot2, scaled from 0:5000
pdf("DE VFDB Clustered Row Heatmap Rhizo.pdf",width = 8, height = 7)
ggplot(count_data_df_R, aes(x = column, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = colors, limits = c(0,2500)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 1)) +
  labs(x = "Samples", y = "Virulence Factors")
dev.off()

pdf("DE VFDB Clustered Row Heatmap Soil.pdf",width = 8, height = 7)
ggplot(count_data_df_S, aes(x = column, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = colors, limits = c(0,2500)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 1)) +
  labs(x = "Samples", y = "Virulence Factors")
dev.off()





