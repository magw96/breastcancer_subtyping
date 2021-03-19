# # lncRNA_1.0 (STAR, gene counts)

x <- list.files(path = "t", pattern = "*.tabular"), full.names = TRUE)

# 0. Load required packages
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(RColorBrewer)
library(dplyr)
library(GO.db)
library(goseq)

# 0. Merge gene expression analysis (.tabular files)

DS <- merge(`AU5656`, `BT20`, by = "V1")
DS <- merge(DS, `BT201`, by = "V1")
DS <- merge(DS, `BT202`, by = "V1")
DS <- merge(DS, `BT474`, by = "V1")
DS <- merge(DS, `BT4741`, by = "V1")
DS <- merge(DS, `BT483`, by = "V1")
DS <- merge(DS, `BT549`, by = "V1")
DS <- merge(DS, `BT5491`, by = "V1")
DS <- merge(DS, `BT5492`, by = "V1")
DS <- merge(DS, `CAL120`, by = "V1")
DS <- merge(DS, `CAL1201`, by = "V1")
DS <- merge(DS, `CAL148`, by = "V1")
DS <- merge(DS, `CAL51`, by = "V1")
DS <- merge(DS, `CAL511`, by = "V1")
DS <- merge(DS, `CAL851`, by = "V1")
DS <- merge(DS, `CAL8511`, by = "V1")
DS <- merge(DS, `CAMA1`, by = "V1")
DS <- merge(DS, `CAMA11`, by = "V1")
DS <- merge(DS, `CL184A13`, by = "V1")
DS <- merge(DS, `CL184B52`, by = "V1")
DS <- merge(DS, `CL600MPE`, by = "V1")
DS <- merge(DS, `DU4475`, by = "V1")
DS <- merge(DS, `DY36T2`, by = "V1")
DS <- merge(DS, `EFM19`, by = "V1")
DS <- merge(DS, `EFM192A`, by = "V1")
DS <- merge(DS, `EVSAT`, by = "V1")
DS <- merge(DS, `HBL51`, by = "V1")
DS <- merge(DS, `HCC1143`, by = "V1")
DS <- merge(DS, `HCC11431`, by = "V1")
DS <- merge(DS, `HCC11432`, by = "V1")
DS <- merge(DS, `HCC11871`, by = "V1")
DS <- merge(DS, `HCC11872`, by = "V1")
DS <- merge(DS, `HCC1395`, by = "V1")
DS <- merge(DS, `HCC13951`, by = "V1")
DS <- merge(DS, `HCC1419`, by = "V1")
DS <- merge(DS, `HCC14191`, by = "V1")
DS <- merge(DS, `HCC1428`, by = "V1")
DS <- merge(DS, `HCC14281`, by = "V1")
DS <- merge(DS, `HCC1500`, by = "V1")
DS <- merge(DS, `HCC15001`, by = "V1")
DS <- merge(DS, `HCC15691`, by = "V1")
DS <- merge(DS, `HCC15692`, by = "V1")
DS <- merge(DS, `HCC1599`, by = "V1")
DS <- merge(DS, `HCC15991`, by = "V1")
DS <- merge(DS, `HCC1806`, by = "V1")
DS <- merge(DS, `HCC18061`, by = "V1")
DS <- merge(DS, `HCC1937`, by = "V1")
DS <- merge(DS, `HCC19371`, by = "V1")
DS <- merge(DS, `HCC19372`, by = "V1")
DS <- merge(DS, `HCC1954`, by = "V1")
DS <- merge(DS, `HCC19541`, by = "V1")
DS <- merge(DS, `HCC19542`, by = "V1")
DS <- merge(DS, `HCC202`, by = "V1")
DS <- merge(DS, `HCC2185`, by = "V1")
DS <- merge(DS, `HCC2218`, by = "V1")
DS <- merge(DS, `HCC2688`, by = "V1")
DS <- merge(DS, `HCC3153`, by = "V1")
DS <- merge(DS, `HCC38`, by = "V1")
DS <- merge(DS, `HCC381`, by = "V1")
DS <- merge(DS, `HCC382`, by = "V1")
DS <- merge(DS, `HCC518`, by = "V1")
DS <- merge(DS, `HCC70`, by = "V1")
DS <- merge(DS, `HCC701`, by = "V1")
DS <- merge(DS, `HCC702`, by = "V1")
DS <- merge(DS, `HCC712`, by = "V1")
DS <- merge(DS, `HDQP1`, by = "V1")
DS <- merge(DS, `HME1`, by = "V1")
DS <- merge(DS, `HS578T`, by = "V1")
DS <- merge(DS, `HS578T1`, by = "V1")
DS <- merge(DS, `JIMT1`, by = "V1")
DS <- merge(DS, `KPL1`, by = "V1")
DS <- merge(DS, `LY2`, by = "V1")
DS <- merge(DS, `MACLS2`, by = "V1")
DS <- merge(DS, `MB157`, by = "V1")
DS <- merge(DS, `MCF10A`, by = "V1")
DS <- merge(DS, `MCF10A1`, by = "V1")
DS <- merge(DS, `MCF10A2`, by = "V1")
DS <- merge(DS, `MCF10A3`, by = "V1")
DS <- merge(DS, `MCF10A4`, by = "V1")
DS <- merge(DS, `MCF12A`, by = "V1")
DS <- merge(DS, `MCF7`, by = "V1")
DS <- merge(DS, `MCF75`, by = "V1")
DS <- merge(DS, `MCF76`, by = "V1")
DS <- merge(DS, `MDAMB134`, by = "V1")
DS <- merge(DS, `MDAMB1341`, by = "V1")
DS <- merge(DS, `MDAMB134VI`, by = "V1")
DS <- merge(DS, `MDAMB157`, by = "V1")
DS <- merge(DS, `MDAMB1571`, by = "V1")
DS <- merge(DS, `MDAMB1572`, by = "V1")
DS <- merge(DS, `MDAMB175VII`, by = "V1")
DS <- merge(DS, `MDAMB231`, by = "V1")
DS <- merge(DS, `MDAMB2310`, by = "V1")
DS <- merge(DS, `MDAMB2311`, by = "V1")
DS <- merge(DS, `MDAMB2315`, by = "V1")
DS <- merge(DS, `MDAMB2316`, by = "V1")
DS <- merge(DS, `MDAMB2317`, by = "V1")
DS <- merge(DS, `MDAMB2318`, by = "V1")
DS <- merge(DS, `MDAMB231d`, by = "V1")
DS <- merge(DS, `MDAMB231p`, by = "V1")
DS <- merge(DS, `MDAMB231t`, by = "V1")
DS <- merge(DS, `MDAMB330`, by = "V1")
DS <- merge(DS, `MDAMB361`, by = "V1")
DS <- merge(DS, `MDAMB3611`, by = "V1")
DS <- merge(DS, `MDAMB3612`, by = "V1")
DS <- merge(DS, `MDAMB415`, by = "V1")
DS <- merge(DS, `MDAMB436`, by = "V1")
DS <- merge(DS, `MDAMB436_T1`, by = "V1")
DS <- merge(DS, `MDAMB436_T10`, by = "V1")
DS <- merge(DS, `MDAMB436_T11`, by = "V1")
DS <- merge(DS, `MDAMB436_T12`, by = "V1")
DS <- merge(DS, `MDAMB436_T13`, by = "V1")
DS <- merge(DS, `MDAMB436_T14`, by = "V1")
DS <- merge(DS, `MDAMB436_T2`, by = "V1")
DS <- merge(DS, `MDAMB436_T3`, by = "V1")
DS <- merge(DS, `MDAMB436_T4`, by = "V1")
DS <- merge(DS, `MDAMB436_T5`, by = "V1")
DS <- merge(DS, `MDAMB436_T6`, by = "V1")
DS <- merge(DS, `MDAMB436_T7`, by = "V1")
DS <- merge(DS, `MDAMB436_T8`, by = "V1")
DS <- merge(DS, `MDAMB436_T9`, by = "V1")
DS <- merge(DS, `MDAMB4361`, by = "V1")
DS <- merge(DS, `MDAMB453`, by = "V1")
DS <- merge(DS, `MDAMB4531`, by = "V1")
DS <- merge(DS, `MDAMB4532`, by = "V1")
DS <- merge(DS, `MDAMB468`, by = "V1")
DS <- merge(DS, `MDAMB4681`, by = "V1")
DS <- merge(DS, `MFM223`, by = "V1")
DS <- merge(DS, `MX1`, by = "V1")
DS <- merge(DS, `OCUBM`, by = "V1")
DS <- merge(DS, `PDX1258`, by = "V1")
DS <- merge(DS, `PDX12582`, by = "V1")
DS <- merge(DS, `PDX1328`, by = "V1")
DS <- merge(DS, `PDX13282`, by = "V1")
DS <- merge(DS, `PDXHCI002`, by = "V1")
DS <- merge(DS, `PDXHCI0021`, by = "V1")
DS <- merge(DS, `PXDTLB`, by = "V1")
DS <- merge(DS, `PXDTLB1`, by = "V1")
DS <- merge(DS, `PXDTLB2`, by = "V1")
DS <- merge(DS, `PXDTLB3`, by = "V1")
DS <- merge(DS, `PXDTLB4`, by = "V1")
DS <- merge(DS, `PXDTLB5`, by = "V1")
DS <- merge(DS, `SKBR3`, by = "V1")
DS <- merge(DS, `SKBR310`, by = "V1")
DS <- merge(DS, `SKBR36`, by = "V1")
DS <- merge(DS, `SKBR37`, by = "V1")
DS <- merge(DS, `SKBR38`, by = "V1")
DS <- merge(DS, `SKBR39`, by = "V1")
DS <- merge(DS, `SUM102`, by = "V1")
DS <- merge(DS, `SUM1021`, by = "V1")
DS <- merge(DS, `SUM1315`, by = "V1")
DS <- merge(DS, `SUM149`, by = "V1")
DS <- merge(DS, `SUM1491`, by = "V1")
DS <- merge(DS, `SUM159`, by = "V1")
DS <- merge(DS, `T47D`, by = "V1")
DS <- merge(DS, `T47D1`, by = "V1")
DS <- merge(DS, `UACC8126`, by = "V1")
DS <- merge(DS, `ZR7514`, by = "V1")
DS <- merge(DS, `ZR75301`, by = "V1")
DS <- merge(DS, `ZR75302`, by = "V1")

# Select only columns of interest (eg, V1) total 159

DS <- DS [,c(1,2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95,98,101,104,107,110,113,116,119,122,125,128,131,134,137,140,143,146,149,152,155,158,161,164,167,170,173,176,179,182,185,188,191,194,197,200,203,206,209,212,215,218,221,224,227,230,233,236,239,242,245,248,251,254,257,260,263,266,269,272,275,278,281,284,287,290,293,296,299,302,305,308,311,314,317,320,323,326,329,332,335,338,341,344,347,350,353,356,359,362,365,368,371,374,377,380,383,386,389,392,395,398,401,404,407,410,413,416,419,422,425,428,431,434,437,440,443,446,449,452,455,458,461,464,467,470,473,476,479)]

# Rename columns

colnames(DS) <- c("ID", "AU5656","BT20", "BT201","BT202", "BT474", "BT474(1)", "BT483", "BT549", "BT549(1)", "BT549(2)", "CAL120", "CAL120(1)", "CAL148", "CAL51", "CAL51(1)", "CAL851", "CAL851(1)", "CAMA1", "CAMA1(1)", "CL184A1(3)", "CL184B5(2)", "CL600MPE", "DU4475", "DY36T2", "EFM19", "EFM192A", "EVSAT", "HBL51", "HCC1143", "HCC1143(1)", "HCC1143(2)", "HCC1187(1)", "HCC1187(2)", "HCC1395", "HCC1395(1)", "HCC1419", "HCC1419(1)", "HCC1428", "HCC1428(1)", "HCC1500", "HCC1500(1)", "HCC1569(1)", "HCC1569(2)", "HCC1599", "HCC1599(1)", "HCC1806", "HCC1806(1)", "HCC1937", "HCC1937(1)", "HCC1937(2)", "HCC1954", "HCC1954(1)", "HCC1954(2)", "HCC202", "HCC2185", "HCC2218", "HCC2688", "HCC3153", "HCC38", "HCC38(1)", "HCC38(2)", "HCC518", "HCC70", "HCC70(1)", "HCC70(2)", "HCC712", "HDQP1", "HME1", "HS578T", "HS578T(1)", "JIMT1", "KPL1", "LY2", "MACLS2", "MB157", "MCF10A", "MCF10A1", "MCF10A2", "MCF10A3", "MCF10A4", "MCF12A", "MCF7", "MCF7(5)", "MCF7(6)", "MDAMB134", "MDAMB1341","MDAMB134VI", "MDAMB157", "MDAMB157(1)", "MDAMB157(2)", "MDAMB175VII", "MDAMB231", "MDAMB231(0)", "MDAMB231(1)", "MDAMB231(5)", "MDAMB231(6)", "MDAMB231(7)", "MDAMB231(8)", "MDAMB231d", "MDAMB231p", "MDAMB231t", "MDAMB330", "MDAMB361", "MDAMB361(1)", "MDAMB361(2)", "MDAMB415", "MDAMB436(1)", "MDAMB436_T1","MDAMB436_T10", "MDAMB436_T11", "MDAMB436_T12", "MDAMB436_T13", "MDAMB436_T14", "MDAMB436_T2", "MDAMB436_T3", "MDAMB436_T4", "MDAMB436_T5", "MDAMB436_T6", "MDAMB436_T7","MDAMB436_T8", "MDAMB436_T9", "MDAMB436(2)", "MDAMB453", "MDAMB453(1)", "MDAMB453(2)", "MDAMB468", "MDAMB468(1)", "MFM223", "MX1", "OCUBM", "PDX1258", "PDX1258(2)", "PDX1328", "PDX1328(2)", "PDXHCI002", "PDXHCI002(2)", "PXDTLB", "PXDTLB1", "PXDTLB2", "PXDTLB3", "PXDTLB4", "PXDTLB5", "SKBR3", "SKBR3(10)", "SKBR3(6)", "SKBR3(7)", "SKBR3(8)", "SKBR3(9)", "SUM102", "SUM102(1)", "SUM1315", "SUM149", "SUM149(1)", "SUM159", "T47D", "T47D(1)", "UACC8126", "ZR751(4)", "ZR7530(1)", "ZR75302")

# Create matrix to feed DESeq2 (c(cell line, condition, gene expression))

dsmatrix <- read.csv("~/Desktop/dsmatrix.csv")
rownames(dsmatrix) <- dsmatrix$Linea
dsmatrix$Linea <- NULL

# DESeq2 dataset generation (may need to convert characters to factors)

dsmatrix$Condicion = factor(dsmatrix$Condicion)
dsmatrix$Tratamiento = as.factor(dsmatrix$Tratamiento)
dsmatrix$Lote = as.factor(dsmatrix$Lote)
dsmatrix$Subtipo = as.factor(dsmatrix$Subtipo)
dsmatrix$Subtipo <- relevel(dsmatrix$Subtipo, "NM")
dsmatrix$Replicas = as.factor(dsmatrix$Replicas)

rownames(DS) <- DS$ID
DS$ID <- NULL

DDS <- DESeqDataSetFromMatrix(countData = DS, colData = dsmatrix, design = ~  Condicion + Subtipo)
dds <- DESeqDataSet(DDS, ~ Condicion + Subtipo)
dds <- collapseReplicates(dds, dds$Replicas)

# DESeq2 analysis
DDS1 <- DESeq(dds)

# Display results
res <- results(DDS1, name = "Condicion_T_vs_NT")

res_shrink <- lfcShrink(DDS1, contrast = contrast_bc, res = res, type = "normal")
par(mfrow=c(1,2))
plotMA(res, ylim = c(-2,2))
plotMA(res_shrink, ylim = c(-2,2))

bcstar05 <- subset(res_shrink, res_shrink$padj < 0.05)
bcstar_cut <- subset(bcstar05, bcstar05$log2FoldChange > 1.5)

write.csv(bc_cut, "../Desktop/bc_cut.csv")

# Filter results with significant p-value
respv <- results(DDS1, alpha=0.01)

# Order results by p-value
resorden <- res2[order(res2$pvalue),]

# Compare cell types (or subtypes)
resultsNames(DDS1)
res.cell <- results(DDS1, contrast = c("Condicion", "T", "NT"))

# Plot results (basic)
plotMA(res, ylim=c(-5,5))

res1 <- subset(res, res$padj < 0.01)
res2 <- subset(res1, res1$log2FoldChange < -5 | res1$log2FoldChange > 5)
res3 <-subset(res_anno, res_anno$log2FoldChange < -5 | res_anno$log2FoldChange > 5)
res4 <- subset(res3, res3$padj < 0001)

# Similar to rlog transformation
vsd <- vst(DDS1, blind = FALSE)

# Determine differences between samples with pheatmap
sampleDists <- dist(t( assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Tratamiento, vsd$Subtipo,sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# And with PCA, in which the x-axis is the principal
# component that separates data the most

pcaData_sal <- plotPCA(vsd, intgroup=c("Condicion", "Subtipo"), returnData=TRUE)
ggplot(pcaData, aes(PC1, PC2, color=Subtipo, shape=Condicion)) +
       geom_point(size=3) +
       scale_shape_manual(values = c(20,10))
       xlab(paste0("PC1: ",percentVar[1],"% variance")) +
       ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
       coord_fixed()

mat <- assay(vsd)
mat <- mat - rowMeans(mat)

idx <- res2$external_gene_name
mat.hm <- mat[idx$idx]

topVarGenes <- head(order(-rowVars(assay(vsd))),35)
mat1 <- assay(vsd)[topVarGenes, ]
mat1 <- mat1 - rowMeans(mat1)
genmat <- as.matrix(geneannot)
rownames(genmat) <- geneannot$SYMBOL
row.names(mat1) <- as.vector(genmat[,2])

cols = colorRampPalette(c("green", "black","red"))(100)
hm <- pheatmap(mat1, annotation_col = dft, colors = cols, fontsize = 6, fontsize_row = 5, fontsize_col = 5, fontsize_number = 5, color = cols)

save.image("~/Desktop/lncRNAs/lncRNAs.RData")
