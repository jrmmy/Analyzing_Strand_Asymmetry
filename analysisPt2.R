## Analysis with different dataset
library(maftools)
library(dplyr)
library(tidyr)
library(TCGAbiolinks)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(tidyverse)


setwd("/Users/jeremydawe/Documents/School/GoutLab/HumanMutations")


query <- GDCquery(project = "TCGA-LIHC",
                   data.category = "Simple Nucleotide Variation",
                   data.type = "Masked Somatic Mutation",
                   workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

gr <- import("Homo_sapiens.GRCh38.87.chr.gff3.gz")

genes <- gr[mcols(gr)$type == "gene"]

GDCdownload(query, method = "api")

LIHCdataAll <- GDCprepare(query)
LIHCdata <- LIHCdataAll %>% 
  select(Chromosome, Start_Position, End_Position, Hugo_Symbol, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification)
LIHCdata <- LIHCdata %>% filter(Start_Position == End_Position,
                                nchar(Reference_Allele) == 1,
                                nchar(Tumor_Seq_Allele2) == 1)

LIHCgr <- GRanges(seqnames = LIHCdata$Chromosome,
                  ranges = IRanges(start = LIHCdata$Start_Position,
                                   end = LIHCdata$End_Position),
                  strand = "*")

mcols(LIHCgr) <- LIHCdata[, c("Chromosome", "Start_Position", "End_Position", "Hugo_Symbol", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification")]
colnames(mcols(LIHCgr))[colnames(mcols(LIHCgr)) == "Tumor_Seq_Allele2"] <- "varAllele"

# where varAllele == "-" there is a frameshift del, get rid of these
LIHCgr <- LIHCgr[mcols(LIHCgr)$varAllele != "-"]

# Add a mutation column to get frequency information before splitting
mcols(LIHCgr)$mut <- paste0(LIHCgr$Reference_Allele, ">", LIHCgr$varAllele)
wholeGenomeCounts <- table(mcols(LIHCgr)$mut)

# Now we need to find genes that are on the tx(+) and tx(-) strands 
# tx(+) = coding strand on the reference 

txPlusGenes <- genes[strand(genes) == "+"]
txMinusGenes <- genes[strand(genes) == "-"]

hitsPlus <- findOverlaps(LIHCgr, txPlusGenes)
hitsMinus <- findOverlaps(LIHCgr, txMinusGenes)

LIHCtxP <- LIHCgr[queryHits(hitsPlus)]
LIHCtxM <- LIHCgr[queryHits(hitsMinus)] # Sum is greater than length(LIHCgr) but this is probably due to duplicates from different datasets

txPcounts <- table(mcols(LIHCtxP)$mut)
txMcounts <- table(mcols(LIHCtxM)$mut)

txPdf <- as.data.frame(txPcounts)
txMdf <- as.data.frame(txMcounts)

# Create paired bar graph for both counts
pairs <- list(
  c("A>C","T>G"),
  c("A>G","T>C"),
  c("A>T","T>A"),
  c("C>A","G>T"),
  c("C>G","G>C"),
  c("C>T","G>A")
)

pair_colors <- c(
  "A>C" = "#08306b", "T>G" = "#9ecae1",
  "A>G" = "#08306b", "T>C" = "#9ecae1",
  "A>T" = "#08306b", "T>A" = "#9ecae1",
  "C>A" = "#08306b", "G>T" = "#9ecae1",
  "C>G" = "#08306b", "G>C" = "#9ecae1",
  "C>T" = "#08306b", "G>A" = "#9ecae1"
)

txPdf$Pair <- sapply(txPdf$Var1, function(m) {
  for (p in pairs) {
    if (m %in% p) {
      return(paste(p, collapse="/"))
    }
  }
  return(NA)  # just in case
})

txMdf$Pair <- sapply(txPdf$Var1, function(m) {
  for (p in pairs) {
    if (m %in% p) {
      return(paste(p, collapse="/"))
    }
  }
  return(NA)  # just in case
})

txPdf$Group <- "tx+"
txMdf$Group <- "tx-"

combined <- rbind(txPdf, txMdf)

ggplot(combined, aes(x = Pair, y = Freq, fill = Var1)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values = pair_colors) +
  facet_wrap(~ Group, ncol = 1) +   # stack vertically
  theme_minimal() +
  labs(y = "Number of Mutations",
       x = "Complementary Mutation Pair",
       title = "Mutation Spectra in tx+ vs tx- Genes")


###### 
# Check antisense transcription
mcols(LIHCgr)$inTxP <- countOverlaps(LIHCgr, txPlusGenes) > 0
mcols(LIHCgr)$inTxM <- countOverlaps(LIHCgr, txMinusGenes) > 0

txP_A <- LIHCgr[mcols(LIHCgr)$inTxP & mcols(LIHCgr)$inTxM]
txP_N <- LIHCgr[mcols(LIHCgr)$inTxP & !mcols(LIHCgr)$inTxM]
txM_A <- LIHCgr[mcols(LIHCgr)$inTxM & mcols(LIHCgr)$inTxP]
txM_N <- LIHCgr[mcols(LIHCgr)$inTxM & !mcols(LIHCgr)$inTxP]

txPAcnts <- (table(txP_A$mut)/length(txP_A))
txPNcnts <- (table(txP_N$mut)/length(txP_N))
txMAcnts <- (table(txM_A$mut)/length(txM_A))
txMNcnts <- (table(txM_N$mut)/length(txM_N))

df <- bind_rows(
  data.frame(mutation = names(txPAcnts), proportion = as.numeric(txPAcnts), group = "tx+ antisense"),
  data.frame(mutation = names(txPNcnts), proportion = as.numeric(txPNcnts), group = "tx+ non-antisense"),
  data.frame(mutation = names(txMAcnts), proportion = as.numeric(txMAcnts), group = "tx- antisense"),
  data.frame(mutation = names(txMNcnts), proportion = as.numeric(txMNcnts), group = "tx- non-antisense")
)

pair_map <- c(
  "A>T"="A>T / T>A", "T>A"="A>T / T>A",
  "A>C"="A>C / T>G", "T>G"="A>C / T>G",
  "A>G"="A>G / T>C", "T>C"="A>G / T>C",
  "C>A"="C>A / G>T", "G>T"="C>A / G>T",
  "C>T"="C>T / G>A", "G>A"="C>T / G>A",
  "C>G"="C>G / G>C", "G>C"="C>G / G>C"
)

df$Pair <- pair_map[df$mutation]

ggplot(df, aes(x = Pair, y = proportion, fill = mutation)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values = pair_colors) +
  facet_wrap(~ group, ncol = 1) +   # lowercase g here
  theme_minimal(base_size = 14) +
  labs(
    y = "Proportion of Mutations",
    x = "Complementary Mutation Pair",
    title = "Mutation Spectra in tx+ vs tx− Genes (Antisense vs Non-Antisense)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )


#######################################################################################################
######################## Combine with coverage information ############################################
#######################################################################################################
#** RUN ALL PREVIOUS CODE

ccs <- readRDS("cc.rds")

# expected = as expression increases, so does the asymmetry because there is more tcr
# Need to divide genes based on tx level in sense, then get the asymmetry for gene from analysis
ccr <- ccs[which(ccs$ratio_anti<10),]
ccr <- ccr %>% mutate(senseBin = ntile(sense, 4))
table(ccr$senseBin)

# fix seqlevels
seqlevels(genes) <- paste0("chr", seqlevels(genes))

# Map mutations to genes
hits <- findOverlaps(LIHCgr, genes)
LIHCgr$geneID <- NA
LIHCgr$geneID[queryHits(hits)] <- mcols(genes)$ID[subjectHits(hits)]
  
LIHCgr$geneID <- sub("^gene:", "", LIHCgr$geneID)

# Add bin information to LIHCgr 
LIHCgr$coverageBin <- ccr$senseBin[match(mcols(LIHCgr)$geneID, ccr$geneID)]
LIHCgr$antisenseRatio <- ccr$ratio[match(mcols(LIHCgr)$geneID, ccr$geneID)]

# Drop NAs in coverageBin 
LIHCgr <- LIHCgr[!is.na(mcols(LIHCgr)$coverageBin)] #*This will change the proportions in the graph

# ReRun antisense code with bin information:

# Make antisense ratio bins
LIHCgr$ratioBin <- 0
mcols(LIHCgr)$ratioBin[mcols(LIHCgr)$antisenseRatio <= 0.05] <- 1
mcols(LIHCgr)$ratioBin[mcols(LIHCgr)$antisenseRatio > 0.9] <- 2

# Choose Bins
wantedCoverageBin <-  4
wantedRatioBin <- 2

BINNEDgr <- LIHCgr[which(mcols(LIHCgr)$coverageBin == wantedCoverageBin ),]
BINNEDgr <- BINNEDgr[which(mcols(BINNEDgr)$ratioBin == wantedRatioBin),]

txP_A <- BINNEDgr[mcols(BINNEDgr)$inTxP & mcols(BINNEDgr)$inTxM]
txP_N <- BINNEDgr[mcols(BINNEDgr)$inTxP & !mcols(BINNEDgr)$inTxM]
txM_A <- BINNEDgr[mcols(BINNEDgr)$inTxM & mcols(BINNEDgr)$inTxP]
txM_N <- BINNEDgr[mcols(BINNEDgr)$inTxM & !mcols(BINNEDgr)$inTxP]

txPAcnts <- (table(txP_A$mut)/length(txP_A))
txPNcnts <- (table(txP_N$mut)/length(txP_N))
txMAcnts <- (table(txM_A$mut)/length(txM_A))
txMNcnts <- (table(txM_N$mut)/length(txM_N))

df <- bind_rows(
  data.frame(mutation = names(txPAcnts), proportion = as.numeric(txPAcnts), group = "tx+ antisense"),
  data.frame(mutation = names(txPNcnts), proportion = as.numeric(txPNcnts), group = "tx+ non-antisense"),
  data.frame(mutation = names(txMAcnts), proportion = as.numeric(txMAcnts), group = "tx- antisense"),
  data.frame(mutation = names(txMNcnts), proportion = as.numeric(txMNcnts), group = "tx- non-antisense")
)
df$Pair <- pair_map[df$mutation]

ggplot(df, aes(x = Pair, y = proportion, fill = mutation)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values = pair_colors) +
  facet_wrap(~ group, ncol = 1) +   # lowercase g here
  theme_minimal(base_size = 14) +
  labs(
    y = "Proportion of Mutations",
    x = "Complementary Mutation Pair",
    title = "Mutation Spectra in tx+ vs tx− Genes (Antisense vs Non-Antisense)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

######################## Check specificity with antisense regions ######################## 

# Get only antisense genes
hits <- findOverlaps(txPlusGenes, txMinusGenes, ignore.strand = TRUE)
antiPlusGenes <- txPlusGenes[unique(queryHits(hits))]
antiMinusGenes <- txMinusGenes[unique(subjectHits(hits))]

# Extract specific regions 
antisenseRegions <- pintersect(txPlusGenes[queryHits(hits)], txMinusGenes[subjectHits(hits)], ignore.strand = T) 

# Get GRanges of mutations within these regions
hits <- findOverlaps(LIHC, antisenseRegions)
LIHCregions <- LIHCgr[queryHits(hits)]

# Get Granges of mutations not in these regions, but in these genes
geneIDsInAntisenseRegions <- unique(antisenseRegions$gene_id)
genesInAntisenseRegions <- genes[genes$gene_id %in% geneIDsInAntisenseRegions]

hits <- findOverlaps(LIHCgr, genesInAntisenseRegions)
LIHCnonRegions <- LIHCgr[queryHits(hits)]

# Exclude genes in the overlapping regions
hits <- findOverlaps(LIHCnonRegions, LIHCregions)
LIHCnonRegions <- LIHCnonRegions[-queryHits(hits)]

# Graph barplots

# Function to compute mutation proportions
getMutationDF <- function(gr, groupName) {
  mut_counts <- table(gr$mut)
  df <- data.frame(
    mutation = names(mut_counts),
    proportion = as.numeric(mut_counts) / length(gr),
    group = groupName
  )
  df$Pair <- pair_map[df$mutation]
  return(df)
}

# Get dataframes
df_regions <- getMutationDF(LIHCregions, "Antisense Region")
df_nonRegions <- getMutationDF(LIHCnonRegions, "Non-Antisense Region")

# Combine
df <- bind_rows(df_regions, df_nonRegions)

# Plot
ggplot(df, aes(x = Pair, y = proportion, fill = pair)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = pair_colors) +
  facet_wrap(~ group, ncol = 1) +  # stack vertically
  theme_minimal(base_size = 14) +
  labs(
    y = "Proportion of Mutations",
    x = "Complementary Mutation Pair",
    title = "Mutation Spectra: Antisense Regions vs Non-Antisense Regions"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

#######################################################################################################
######################## Check lnc RNA w/ and w/out antisense #########################################
#######################################################################################################

lncGr <- import("gencode.v38.long_noncoding_RNAs.gff3.gz")

pc <- gr[which(gr$type=="gene" & gr$biotype=="protein_coding")]
seqlevels(pc) <- paste0("chr", seqlevels(pc))


fo <- findOverlaps(lncGr, pc, ignore.strand = T)

lncGr_isolated <- lncGr[-queryHits(fo)]
lncGr <- GenomicRanges::reduce(lncGr_isolated)
sum(width(lncGr))/1e9


# Label them with mutation information
hits <- findOverlaps(LIHCgr, lncGr, ignore.strand = T)
LIHClnc <- LIHCgr[queryHits(hits)]

# Get txPlus & txMinus for lnc
lncTxPlus <- lncGr[strand(lncGr) == "+"]
lncTxMinus <- lncGr[strand(lncGr) == "-"]

# Label the mutations GRanges
mcols(LIHClnc)$inTxP <- countOverlaps(LIHClnc, lncTxPlus) > 0
mcols(LIHClnc)$inTxM <- countOverlaps(LIHClnc, lncTxMinus) > 0

# Subset for the mutations on Plus/Minus strand
txP <- LIHClnc[mcols(LIHClnc)$inTxP]
txM <- LIHClnc[mcols(LIHClnc)$inTxM]

# Count each mutation
txPcnts <- (table(txP$mut)/length(txP))
txMcnts <- (table(txM$mut)/length(txM))

df <- bind_rows(
  data.frame(mutation = names(txPcnts), proportion = as.numeric(txPcnts), group = "tx+"),
  data.frame(mutation = names(txMcnts), proportion = as.numeric(txMcnts), group = "tx-")
)


pair_map <- c(
  "A>T"="A>T / T>A", "T>A"="A>T / T>A",
  "A>C"="A>C / T>G", "T>G"="A>C / T>G",
  "A>G"="A>G / T>C", "T>C"="A>G / T>C",
  "C>A"="C>A / G>T", "G>T"="C>A / G>T",
  "C>T"="C>T / G>A", "G>A"="C>T / G>A",
  "C>G"="C>G / G>C", "G>C"="C>G / G>C"
)

df$Pair <- pair_map[df$mutation]

ggplot(df, aes(x = Pair, y = proportion, fill = mutation)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values = pair_colors) +
  facet_wrap(~ group, ncol = 1) +   # lowercase g here
  theme_minimal(base_size = 14) +
  labs(
    y = "Proportion of Mutations",
    x = "Complementary Mutation Pair",
    title = "Mutation Spectra in tx+ vs tx− lncRNA"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
