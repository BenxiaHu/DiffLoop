rm(list=ls())
library(DESeq2)
library(tidyr)
library(writexl)
library(ggplot2)
library(optparse)

option_list <- list(
    make_option("--Pathid", type = "character"),
    make_option("--Rawcounts", type = "integer"),
    make_option("--Normalization", type = "integer"),
    make_option("--Treatid", type = "character"),
    make_option("--Ctrlid", type = "character"),
    make_option("--Treatnum", type = "integer"),
    make_option("--Ctrlnum", type = "integer"),
    make_option("--outfile", type = "character"),
    make_option("--log2FC", type = "integer"),
    make_option("--padj", type = "integer"),
)

opt <- parse_args(OptionParser(option_list = option_list))

countmatrix <- opt$Rawcounts
Treatid <- opt$Treatid
Ctrlid <- opt$Ctrlid

pathid <- opt$Pathid
Normalization <- opt$Normalization
output <- opt$outfile

log2FCvalue <- opt$log2FC
padjvalue <- opt$padj
output <- opt$outfile

Treatnum <- opt$Treatnum
Ctrlnum <- opt$Ctrlnum
setwd(pathid)
#Load gene count matrix and labels

countData <- as.matrix(read.table(countmatrix,header=T,sep="\t",row.names="loop_id"))

colData <- data.frame(condition = factor(c(rep(Treatid, Treatnum),rep(Ctrlid, Ctrlnum))),
                      row.names=c(paste0(Treatid, "_rep", seq_len(Treatnum)),paste0(Ctrlid, "_rep", seq_len(Ctrlnum))))

colData$condition <- relevel(colData$condition, ref = Ctrlid)

sample_ids <- rownames(colData)
normFactors <- read.table(Normalization,header=T,sep="\t",row.names='loop_id')
normFactors <- normFactors[rownames(countData), rownames(colData), drop = FALSE]

combined_scaling <- 1 / as.matrix(normFactors)

countData <- countData[, rownames(colData)]
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData, design = ~ condition)
normalizationFactors(dds) <- as.matrix(combined_scaling)
dds <- DESeq(dds)

res <- results(dds,contrast=c("condition",Treatid,Ctrlid))
resOrdered <- res[order(res$padj), ]
resOrdered <- na.omit(resOrdered)
df <- data.frame(x = rownames(resOrdered))

df2 <- df %>% separate(x, into = c("chrom1", "start1", "end1","chrom2", "start2", "end2"),sep = "[:|\\-]",convert = TRUE)

resOrdered <- data.frame(df2,resOrdered)
write.table(file=paste0(output,"_EPloop_DESeq2.txt"),resOrdered,sep="\t",quote=F,row.names=F)
write_xlsx(resOrdered, paste0(output,"_EPloop_DESeq2.xlsx"))

up <- nrow(resOrdered[resOrdered$padj <= padjvalue & resOrdered$log2FoldChange> log2FCvalue,])
down <- nrow(resOrdered[resOrdered$padj <= padjvalue & resOrdered$log2FoldChange< (-log2FCvalue),])

input <- na.omit(resOrdered)
input$colname <- c("gray")
input$pchname <- c(18)
#print(head(input))
input$colname[input$log2FoldChange > log2FCvalue & input$padj < padjvalue] <- "red"
input$pchname[input$log2FoldChange > log2FCvalue & input$padj < padjvalue] <- 18
input$colname[input$log2FoldChange < (-log2FCvalue) & input$padj < padjvalue] <- "blue"
input$pchname[input$log2FoldChange < (-log2FCvalue) & input$padj < padjvalue] <- 18

pdf(paste0(output,"_EPloop_DEseq2.pdf"))
p <- ggplot(input, aes(x=log2FoldChange, y=-log10(padj), color=colname)) +
  geom_point()+ ggtitle(paste0(output,"\nup=",up,";down=",down)) +
  scale_color_manual(labels = c("down", "non", "up"), values = c("blue","gray", "red"))+
  geom_hline(yintercept = -log10(padjvalue), color="black", linetype="dashed")+
  geom_vline(xintercept = log2FCvalue,color="#009999", linetype="dashed")+
  geom_vline(xintercept = -log2FCvalue, color="#009999", linetype="dashed")+
  theme_classic()
print(p)
dev.off()

