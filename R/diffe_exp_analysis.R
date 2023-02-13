# Load the RSE
load('data/rse_gene_SRP048604.RData')
library("edgeR")
library("ggplot2")
library("limma")
library("pheatmap")

# Expand the sra sample attributes
rse_gene_SRP048604 <- expand_sra_attributes(rse_gene_SRP048604)

# Look at sra sample attributes
colData(rse_gene_SRP048604)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP048604)))
]

#Generate a colData atribute with the interest atribute
rse_gene_SRP048604$culture <- factor(rse_gene_SRP048604$`sra_attribute.co-culture`)
rse_gene_SRP048604$culture <- relevel(rse_gene_SRP048604$culture, ref = 'OP9-GFP')
rse_gene_SRP048604$culture

# Make control quality assessment of samples
rse_gene_SRP048604$assigned_gene_prop <- rse_gene_SRP048604$recount_qc.gene_fc_count_all.assigned/rse_gene_SRP048604$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP048604$assigned_gene_prop)


# Visualize the quality of sample counts
# we can see that all samples are of high quality
hist(rse_gene_SRP048604$assigned_gene_prop, xlab = "Assigned gene proportion", main = 'Assigned gene proportion')

# Check the expression levels
gene_means <- rowMeans(assay(rse_gene_SRP048604, "counts"))
summary(gene_means)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0       0.0       1.5    2095.0     197.0 1355387.2

# Eliminate genes with low expression
rse_gene_SRP048604_unfiltered <- rse_gene_SRP048604
rse_gene_SRP048604 <- rse_gene_SRP048604[gene_means > 0.1, ]

# Check percentage of retained genes
round(nrow(rse_gene_SRP048604) / nrow(rse_gene_SRP048604_unfiltered) * 100, 2)
# 61.89

# Normalize data
dge <- DGEList(
  counts = assay(rse_gene_SRP048604, "counts"),
  genes = rowData(rse_gene_SRP048604)
)
dge <- calcNormFactors(dge)

# Visualize expression distribution in samples
ggplot(as.data.frame(colData(rse_gene_SRP048604)), aes(y = assigned_gene_prop, x = culture)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Culture grouo")
# Generate statistics model
mod <- model.matrix(~ culture + assigned_gene_prop,
                    data = colData(rse_gene_SRP048604)
)
colnames(mod)
mod
# Differential expression analysis
vGenes <- voom(dge, mod, plot = T)

eb_results <- eBayes(lmFit(vGenes))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP048604),
  sort.by = "p"
)
de_results
dim(de_results)
#39519    16

# Number of DEGs
table(de_results$P.Value < 0.05)
#FALSE  TRUE
#35529  3990

# Plots of DEGs
plotMA(eb_results, coef = 2)
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)

# Get the top 50 genes
exprs_heatmap <- vGenes$E[rank(de_results$logFC) <= 50, ]
exprs_heatmap


# Generate a heatmap of
df <- as.data.frame(colData(rse_gene_SRP048604)[, c("sra_attribute.co-culture")])
colData(rse_gene_SRP048604)[, "sra_attribute.co-culture"]

colnames(df) <- c("culture")
rownames(df) <- c('SRR1596221', 'SRR1596222','SRR1596223', 'SRR1596224')
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)

## Create plots MDS
plotMDS(vGenes$E, labels = df$culture, col = c("blue","red"))

