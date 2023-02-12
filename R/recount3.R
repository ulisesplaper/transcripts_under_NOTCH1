## Download experiment data and generate
# RangeSummarizedObject

# Import packages
library("recount3")
library('iSEE')

## Current url
getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)

## Change URL
options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

## Check changes.
getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)

## Check human projects
human_projects <- available_projects(organism = 'human')

## Find project of interest
proj_info <- subset(
  human_projects,
  project == "SRP048604" & project_type == "data_sources"
)
## Create RangedSummarizedExperiment (RSE) object
## with project info
rse_gene_SRP048604 <- create_rse(proj_info)

# Explore RSE
rse_gene_SRP048604

#Scale coverage counts

assay(rse_gene_SRP048604, "counts") <- compute_read_counts(rse_gene_SRP048604)

# Save the scaled rse
save(rse_gene_SRP048604, file = 'data/rse_gene_SRP048604.RData')


# Explore the rse object
# Look at the column info
colData(rse_gene_SRP048604)
# Look at the row info
rowRanges(rse_gene_SRP048604)
# RSE sra sample attributes
rse_gene_SRP048604$sra.sample_attributes

#Visualize interactively the rse obj
iSEE::iSEE(rse_gene_SRP048604)

# RSE experiment title
rse_gene_SRP048604$sra.experiment_title
