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
