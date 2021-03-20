# ------------------------ DE_1.1 - analysis performed using R Version 4.0.3 ------------------------

# BC Differential Expression Analysis: second iteration (using Salmon) to quantify transcripts from
# 158 Breast Cancer cell lines, obtained through the NCBI SRA database. 
# (QC, trimming and alignment performed using STAR (first iteration) and Salmon at: http://usegalaxy.eu

# lncRNA_2.0 (Salmon, transcripts)
# DYT36T2, MCF75 missing salmon quant

x <- list.files(path = "t", pattern = "*.tabular"), full.names = TRUE)

# 0. Load required packages
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(RColorBrewer)
library(dplyr)
library(GO.db)
library(goseq)
library(biomaRt)

# 1. Import data: merge Salmon transcript quantification
DSal <- merge(`AU5656-S`[c("Name", "NumReads")], `BT20-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `BT201-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `BT202-S`[c("Name", "NumReads")], by = "Name")
DSal<- merge(DSal, `BT474-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `BT4741-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `BT483-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `BT549-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `BT5491-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `BT5492-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CAL120-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CAL1201-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CAL148-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CAL51-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CAL511-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CAL851-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CAL8511-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CAMA1-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CAMA11-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CL184A13-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CL184B52-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `CL600MPE-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `DU4475-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `EFM19-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `EFM192A-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `EVSAT-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HBL51-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC1143-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC11431-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC11432-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC11871-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC11872-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC1395-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC13951-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC1419-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC14191-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC1428-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC14281-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC1500-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC15001-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC15691-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC15692-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC1599-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC15991-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC1806-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC18061-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC1937-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC19371-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC19372-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC1954-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC19541-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC19542-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC202-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC2185-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC2218-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC2688-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC3153-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC38-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC381-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC382-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC518-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC70-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC701-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC702-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HCC712-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HDQP1-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HME1-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HS578T-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `HS578T1-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `JIMT1-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `KPL1-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `LY2-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MACLS2-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MB157-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MCF10A-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MCF10A1-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MCF10A2-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MCF10A3-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MCF10A4-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MCF12A-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MCF7-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MCF76-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB134-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB1341-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB134VI-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB157-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB1571-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB1572-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB175VII-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB2310-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB2311-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB2312-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB2315-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB2316-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB2317-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB2318-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB231d-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB231p-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB231t-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB330-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB361-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB3611-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB3612-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB415-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T1-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T10-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T11-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T12-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T13-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T14-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T2-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T3-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T4-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T5-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T6-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T7-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T8-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB436_T9-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB4361-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB453-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB4531-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB4532-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB468-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MDAMB4681-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MFM223-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `MX1-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `OCUBM-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PDX1258-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PDX12582-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PDX1328-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PDX13282-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PDXHCI002-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PDXHCI0021-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PXDTLB-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PXDTLB1-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PXDTLB2-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PXDTLB3-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PXDTLB4-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `PXDTLB5-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SKBR3-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SKBR310-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SKBR36-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SKBR37-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SKBR38-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SKBR39-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SUM102-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SUM1021-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SUM1315-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SUM149-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SUM1491-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `SUM159-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `T47D-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `T47D1-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `UACC8126-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `ZR7514-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `ZR75301-S`[c("Name", "NumReads")], by = "Name")
DSal <- merge(DSal, `ZR75302-S`[c("Name", "NumReads")], by = "Name")

colnames(DSal) <- c("ID", "AU5656","BT20", "BT201","BT202", "BT474", "BT474(1)", "BT483", "BT549", "BT549(1)", 
"BT549(2)", "CAL120", "CAL120(1)", "CAL148", "CAL51", "CAL51(1)", "CAL851", "CAL851(1)", "CAMA1", "CAMA1(1)", 
"CL184A1(3)", "CL184B5(2)", "CL600MPE", "DU4475", "EFM19", "EFM192A", "EVSAT", "HBL51", "HCC1143", "HCC1143(1)", 
"HCC1143(2)", "HCC1187(1)", "HCC1187(2)", "HCC1395", "HCC1395(1)", "HCC1419", "HCC1419(1)", "HCC1428", "HCC1428(1)", 
"HCC1500", "HCC1500(1)", "HCC1569(1)", "HCC1569(2)", "HCC1599", "HCC1599(1)", "HCC1806", "HCC1806(1)", "HCC1937", 
"HCC1937(1)", "HCC1937(2)", "HCC1954", "HCC1954(1)", "HCC1954(2)", "HCC202", "HCC2185", "HCC2218", "HCC2688", "HCC3153", 
"HCC38", "HCC38(1)", "HCC38(2)", "HCC518", "HCC70", "HCC70(1)", "HCC70(2)", "HCC712", "HDQP1", "HME1", "HS578T", "HS578T(1)", 
"JIMT1", "KPL1", "LY2", "MACLS2", "MB157", "MCF10A", "MCF10A1", "MCF10A2", "MCF10A3", "MCF10A4", "MCF12A", "MCF7", "MCF7(6)",
"MDAMB134", "MDAMB1341","MDAMB134VI", "MDAMB157", "MDAMB157(1)", "MDAMB157(2)", "MDAMB175VII", "MDAMB231", "MDAMB231(0)", 
"MDAMB231(1)", "MDAMB231(5)", "MDAMB231(6)", "MDAMB231(7)", "MDAMB231(8)", "MDAMB231d", "MDAMB231p", "MDAMB231t", 
"MDAMB330", "MDAMB361", "MDAMB361(1)", "MDAMB361(2)", "MDAMB415", "MDAMB436(1)", "MDAMB436_T1","MDAMB436_T10", 
"MDAMB436_T11", "MDAMB436_T12", "MDAMB436_T13", "MDAMB436_T14", "MDAMB436_T2", "MDAMB436_T3", "MDAMB436_T4", 
"MDAMB436_T5", "MDAMB436_T6", "MDAMB436_T7","MDAMB436_T8", "MDAMB436_T9", "MDAMB436(2)", "MDAMB453", 
"MDAMB453(1)", "MDAMB453(2)", "MDAMB468", "MDAMB468(1)", "MFM223", "MX1", "OCUBM", "PDX1258", "PDX1258(2)", 
"PDX1328", "PDX1328(2)", "PDXHCI002", "PDXHCI002(2)", "PXDTLB", "PXDTLB1", "PXDTLB2", "PXDTLB3", "PXDTLB4", 
"PXDTLB5", "SKBR3", "SKBR3(10)", "SKBR3(6)", "SKBR3(7)", "SKBR3(8)", "SKBR3(9)", "SUM102", "SUM102(1)", 
"SUM1315", "SUM149", "SUM149(1)", "SUM159", "T47D", "T47D(1)", "UACC8126", "ZR751(4)", "ZR7530(1)", "ZR75302")

# 1.1. Import metadata (to feed DESeqMatrix)
dsalmatrix <- read.csv("~/Desktop/Investigación/INCan/lncRNAs/dsalmatrix.csv")
rownames(dsalmatrix) <- dsalmatrix$Linea
dsalmatrix$Linea <- NULL

# 1.2. DESeq2 dataset generation (convert characters to factors; specify levels for desired experimental design)
dsalmatrix$Condicion = factor(dsalmatrix$Condicion)
dsalmatrix$Condicion <- relevel(dsalmatrix$Condicion, "NT")
dsalmatrix$Tratamiento = as.factor(dsalmatrix$Tratamiento)
dsalmatrix$Tratamiento <- relevel(dsalmatrix$Tratamiento, "N")
dsalmatrix$Lote = as.factor(dsalmatrix$Lote)
dsalmatrix$Subtipo = as.factor(dsalmatrix$Subtipo)
dsalmatrix$Subtipo <- relevel(dsalmatrix$Subtipo, "NM")
dsalmatrix$Replicas = as.factor(dsalmatrix$Replicas)
rownames(DSal) <- gsub("\\..*","", rownames(DSal))
rownames(DSal) <- DSal$ID
DSal$ID <- NULL

# 1.2.1. Check if colnames(raw_counts) match rownames(metadata)
all(colnames(DSal) == rownames(dsalmatrix))

# 1.3. DESeqDataSet creation (define experimental design)
DDSal <- DESeqDataSetFromMatrix(countData = round(DSal), colData = dsalmatrix, design = ~  Condicion + Subtipo)
ddsal <- DESeqDataSet(DDSal, ~ Condicion + Subtipo)
ddsal <- collapseReplicates(ddsal, ddsal$Replicas)

# 1.4. Run DESeq2 analysis (158 samples, ~30 min)
DDSal1 <- DESeq(ddsal)
rownames(DDSal1) <- gsub("\\..*","", rownames(DDSal1))

# 2.0. Display results - alter DE analysis by defining contrast
res_sal <- results(DDSal1, name = "Condicion_T_vs_NT")
# An alternative to name is contrast_bc = c("Condicion", "T", "NT")
contrast_bc = c("Condicion", "T", "NT")
# Where the second term is used as the control

# 2.1. LFC Shrinkage: removes noise from low-count genes with high variances relative to the mean
# (158 samples, ~30 min)
res_sal_shrink <- lfcShrink(DDSal1, contrast = contrast_bc, res = res_sal, type = "normal")
par(mfrow=c(1,2))

# 2.2. Compare dispersion and variance (shrinked vs normal)
plotMA(res_sal, ylim = c(-2,2))
plotMA(res_sal_shrink, ylim = c(-2,2))

# 2.3. Eliminate gene update version number to allow matching using the attribute "ensembl_transcript_id"
# alternatively, use getBM::attributes that includes gene update number (pending, need to add)
rownames(res_sal_shrink) <- gsub("\\..*","", rownames(res_sal_shrink))

# 3.0. Annotate transcripts; using biomaRt increases annotated transcripts significantly
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
listAttributes(mart)
trans.ids <- getBM(attributes = c("external_transcript_name","ensembl_transcript_id"), 
                   filters = "ensembl_transcript_id", values = rownames(res_sal_shrink), 
                   mart = ensembl)

# 3.1. Create new column for annotated transcripts
res_sal_shrink$transannot<- rownames(res_sal_shrink)
colnames(trans.ids)[2] <- "transannot"

# 3.2. Create df to merge annotation file with counts file
res_sal_shrink$transannot<- rownames(res_sal_shrink)
colnames(trans.ids)[2] <- "transannot"
res_sal_df <- as.data.frame(res_sal_shrink)
res_sal_df <- merge(res_sal_df, trans.ids, by = "transannot")

# 4.0. Subsetting DE genes by significance and L2FC
bc05 <- subset(res_sal_df, res_sal_df$padj < 0.05)
bc_cut <- subset(bc05, bc05$log2FoldChange > 1.5)
# 4.0.1. Order by highest L2FC
bc_l2fc <- head(bc_cut[order(-bc_cut$log2FoldChange),], 100)
bc_cut_mat <- as.matrix(bc_cut)
rownames(bc_cut_mat) <- bc_cut$external_transcript_name

# 4.1. Export .csv file
write.csv(bc_cut, "../Desktop/bc_cut.csv")

# DO NOT RUN: Gene Ontology (WIP)
gobc = bc_cut
rownames(gobc) = rownames(gobc[gobc$log2FoldChange!=0,])
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

# 5.0. vsd transformation for heatmap & PCA generation (158 samples, ~1 min)
vsd_sal <- vst(DDSal1, blind = FALSE)

# 5.0.1. Matrix generation (vsd_sal is still a DESeqObject)
mat_sal <- assay(vsd_sal)
mat_sal <- mat_sal - rowMeans(mat_sal)

# 5.0.2. Remove gene version - same as above
rownames(mat_sal) <- gsub("\\..*","", rownames(mat_sal))

# 5.1. Use ggplot2 for PCAplot - define groups and shapes
pcaData_sal <- plotPCA(vsd_sal, intgroup=c("Condicion", "Subtipo"), returnData=TRUE)
ggplot(pcaData_sal, aes(PC1, PC2, color=Subtipo, shape=Condicion)) +
  geom_point(size=3) +
  scale_shape_manual(values = c(20,10))

# 5.2. Filter transcripts with highest variance (in descending order)
# (add1). Filtering by (res_sal_shrink$log2FoldChange) -> gets the most DE genes, 
# i.e, genes with the biggest differences between specific, predefined groups
# (add2). Filtering by (rowVars(assay(vsd_sal)) -> gets highly variable genes, 
# i.e, genes most variable across all samples, regardless of which samples they are
topVarTranscripts <- head(order(-rowVars(assay(vsd_sal))),50)
mat_sal_var <- assay(vsd_sal)[topVarTranscripts, ]

# 5.2.1. Remove gene version - same as above (skip if not necessary)
rownames(mat_sal_var) <- gsub("\\..*","", rownames(mat_sal_var))

mat_sal_var <- mat_sal_var - rowMeans(mat_sal_var)
trans.mat <- as.matrix(trans.ids)
rownames(trans.mat) <- trans.ids$external_transcript_name
row.names(mat_sal) <- as.vector(trans.mat[,2])

# 6.0. Data visualization
# 6.0.1. Build ID'd matrix for easier heatmap annotation
df_sal <- as.data.frame(mat_sal)
df_sal$transannot <- rownames(df_sal)
trans.l2fc <- as.data.frame(bc_l2fc$transannot)
trans.l2fc$external_transcript_name <- bc_l2fc$external_transcript_name
colnames(trans.l2fc)[1] <- "transannot"
x <- merge(df_sal, trans.l2fc, by = "transannot")
rownames(x) <- x$external_transcript_name
x$transannot <- NULL
x$external_transcript_name <- NULL
x <- as.matrix(x)

# --- Create color palette
cols = colorRampPalette(c("green", "black","red"))(100)

# 6.2. Plot your heatmap
hm_sal1 <- pheatmap(x, annotation_col = dsalmatrix, color = cols, fontsize = 6, 
                    fontsize_row = 5, fontsize_col = 5, main = "DE transcripts: treated vs untreated")
