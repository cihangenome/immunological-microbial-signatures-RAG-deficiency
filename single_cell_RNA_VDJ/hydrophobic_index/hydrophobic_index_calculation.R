# This code was used to determine the hydrophobic index.
# Dr. Stephen Daley, Centre for Immunology and Infection Control at Queensland University of Technology (s5.daley@qut.edu.au)

# Install and/or load "tidyverse" and "stringr" packages.
if (!require("tidyverse"))
  install.packages("tidyverse")
library(tidyverse)

if (!require("stringr"))
  install.packages("stringr")
library(stringr)

# With the conserved TCR Vbeta germline-encoded cysteine in CDR3beta defined as position 1 of CDR3beta, 
# Stadinski et al NI 2016 (PMID: 27348411) studied effects on TCR self-reactivity 
# of amino acid doublets at positions 6 and 7 (P6-P7) of CDR3beta. 
# Figure 4i of Stadinski et al NI 2016 (PMID: 27348411) was used to generate a dataframe,
# which categorises each of the 400 possible doublets as: 
# "S" for "strong" or promoting TCR self-reactivity (red in Figure 4i)
# "N" for "neutral" (white in Figure 4i) or
# "W" for "weak" or reducing self-reactivity (blue in Figure 4i).

# Download the dataframe called "stadinski.txt" and save it in the folder corresponding to working_directory
# Read it in to R:
stad <- read.delim("stadinski.txt", header = T)

# Note: # As the 'NA' doublet might wrongly be designated as 'not available', 'NA' is represented as 'NX' in this dataframe.
# Make a list of the 175 amino acid doublets found to promote TCR self-reactivity.
stad_S <- filter(stad, reactivity == "S")

# Load scRNAseq data
# Set "toi" (tissue of interest)
toi <- "laminapropria_10X_scRNASeq_VDJ_metadata"

# Set working directory to the folder that contains data
#data_folder <- "" # **Instruction: Make "data_folder" the path to the folder containing the data.
# The preceding line is the path to the folder on the author's computer  (delete this line)
#setwd(data_folder)

data_folder <-"./"

# Load data
r <- read.delim(file = paste0(data_folder,toi,".txt"), header = T)

# Return to the original working directory
#setwd(working_directory)

# CDR3 aa sequences are in the column called "CTaa"
# Each row contains scRNAseq data from a single cell. 
# Exclude cells (rows) without any CDR3 sequence
r$nchar_aa <- nchar(r$CTaa)
r1 <- r %>% filter(nchar_aa > 0)
# In the column called "CTaa", CDR3A and CDR3B aa sequences are separated by "_".
# Extract CDR3B
r1$cdr3B <- substr(r1$CTaa, str_locate(r1$CTaa, "_") + 1, r1$nchar_aa)

# Get columns required for hydrophobic index calculations
dfs <- r1 %>% select(orig.ident,  cdr3B) %>%
  # and exclude rows that do not have a CDR3B sequence
  filter(!(cdr3B == "NA"))

# Change colname of cdr3B to match previously written code
colnames(dfs) <- c("orig.ident", "cdr3aa")
dfs$len <- as.integer(nchar(dfs$cdr3aa))

# Extract the P6-P7 doublet from each CDR3 sequence
dfs["p6"] <- substr(dfs$cdr3aa, 6, 6)
dfs["p7"] <- substr(dfs$cdr3aa, 7, 7)
dfs$p67 <- paste0(dfs$p6, dfs$p7)

# Change doublets wrongly designated as 'NA' for 'not available' to 'NX'
dfs$p67[dfs$p67 == "NA"] <- "NX"

# Identify doublets that correspond to the 175 doublets found to promote T cell self-reactivity
dfs$S_d <- dfs$p67 %in% stad_S$doub 

results <- dfs %>%
  group_by(orig.ident, S_d) %>%
  summarise(hyd_reads = n()) %>% 
  pivot_wider(names_from = S_d, values_from = hyd_reads, values_fill = 0) %>%
  rename(hyd_present = "TRUE",hyd_absent = "FALSE") %>%
  mutate(hyd_index = 100 * hyd_present / (hyd_present + hyd_absent)) 
