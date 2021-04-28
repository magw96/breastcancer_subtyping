# Breast Cancer Subtyping of RNA-Seq data using genefu (AIMS, SSP, PAM50)

library(genefu)
library(tidyverse)

# --- Working with salmon quant.sf files
# --- Import Salmon counts (using countsfromAbundance = LengthScaledTPM)
txi_f <- tximport(files, type="salmon", tx2gene=t2g[,c("ensembl_transcript_id", "ensembl_gene_id")], 
                  countsFromAbundance="lengthScaledTPM", ignoreTxVersion = T)

# --- Obtain Salmon abundance counts (length-scaled TPM)
# --- log10-scaled; + 1 to avoid problems with values equal to 0
bc_subtype <- log10(txi_f$abundance + 1)

# --- Assign column names to your samples
colnames(bc_subtype) <- c('AU565(6)', 'BT474', 'BT474(1)', 'BT483', 'BT549', 'CAL120', 'CAL851', 'CAMA1', 'CL184A1(3)', 'CL184B5(2)', 'EFM192A', 'EVSAT', 'HCC1143', 'HCC1395', 'HCC1419', 'HCC1428', 'HCC1500', 'HCC1569(1)', 'HCC1599', 'HCC1806', 'HCC1937', 'HCC1954', 'HCC202', 'HCC2185', 'HCC2218', 'HCC712', 'HME1', 'JIMT1', 'MCF10A', 'MCF10A1', 'MCF10A2', 'MCF10A3', 'MCF10A4', 'MCF12A', 'MCF7', 'MDAMB134', 'MDAMB175VII', 'MDAMB330', 'MDAMB361', 'MDAMB361(1)', 'MDAMB415', 'SKBR3', 'SKBR3(9)', 'SUM102', 'SUM1315', 'T47D', 'UACC8126', 'ZR751(4)', 'ZR7530(1)')

# --- Create the annotation file which will be used as input for genefu()
bc_annot <- as.data.frame(rownames(bc_subtype))
colnames(bc_annot) <- "ensembl_gene_id"

# --- Merge annotation file with gene_biotype (which contains c("ensembl_gene_id", "gene_biotype" & "entrez_gene_id")
# --- gene_biotype was previously obtained; annotated using biomaRt
bc_annot <- merge(bc_annot, gene_biotype, by = "ensembl_gene_id", no.dups = T)

# --- Remove duplicated EntrezGene IDs (as they can cause problems when used for identifying rows or columns)
bc_annot <- filter(bc_annot, !is.na(bc_annot$entrezgene_id))
bc_annot <- filter(bc_annot, !duplicated(bc_annot$entrezgene_id))

# --- Rename columns so they match the input genefu() expects; add "probe", which should be == to colnames of
# --- your data matrix
colnames(bc_annot) <- c("Gene.Symbol", "Gene.Biotype", "EntrezGene.ID")
bc_annot$probe <- bc_annot$Gene.Symbol

# --- Convert data frame to matrix (all values should be positive)
subtype_mat <- as.matrix(bc_subtype)

# --- Subset only genes that are contained in bc_annot (so nrow(data) == length(EntrezID))
subtype_mat <- bc_subtype[bc_annot$Gene.Symbol,]

# --- Use molecular.subtyping(); sbt.model can be changed according to your classification needs
# --- The data matrix must be transposed (genes [or probes] should be columns, samples should be rows)
# --- do.mapping = TRUE should be used when the genes (or probes) provided are not limited to those that 
# --- the function uses for its classification
pam50.subtypes <- molecular.subtyping(sbt.model = "pam50", data = t(subtype_mat), annot = bc_annot, do.mapping = T) 
ssp.2003.subtypes <- molecular.subtyping(sbt.model = "ssp2003", data = t(subtype_mat), annot = bc_annot, do.mapping = T)
ssp.subtypes <- molecular.subtyping(sbt.model = "ssp2006", data = t(subtype_mat), annot = bc_annot, do.mapping = T)
aims.subtypes <- molecular.subtyping(sbt.model = "AIMS", data = t(subtype_mat), annot = bc_annot, do.mapping = T)

# --- Visualize classifications

table(pam50.subtypes$subtype)
# Basal   Her2   LumB   LumA Normal 
# 17     13     13      2      4 

table(ssp.subtypes$subtype) # --- This subtyping most closely resembles cell line metadata
# Basal   Her2   LumB   LumA Normal 
# 17      9     11      9      3 

table(aims.subtypes$subtype)
# Basal  Her2  LumB 
# 22    21     6 

table(ssp.2003.subtypes$subtype)
# Basal   Her2   LumB   LumA Normal 
# 19      2     10     15      3 

table(ssp.subtypes$subtype, matrix$Subtipo)
#         NM HER2 LA LB TNBC
# Basal   5    3  0  0    9
# Her2    0    6  0  3    0
# LumB    0    1  8  2    0
# LumA    2    0  2  5    0
# Normal  2    0  0  0    1
