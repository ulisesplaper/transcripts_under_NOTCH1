# Load the RSE
load('data/rse_gene_SRP048604.RData')

# Expand the sra sample attributes
rse_gene_SRP048604 <- expand_sra_attributes(rse_gene_SRP048604)

# Look at sra sample attributes
colData(rse_gene_SRP048604)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP048604)))
]

#Generate a colData atribute with the interest atribute
rse_gene_SRP048604$culture <- factor(rse_gene_SRP048604$`sra_attribute.co-culture`)
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



