# transcripts_under_NOTCH1

## Biological background

Genetic studies in T-cell acute lymphoblastic leukemia have uncovered a remarkable complexity of oncogenic and loss-of-function mutations. Amongst this plethora of genetic changes, NOTCH1 activating mutations stand out as the most frequently occurring genetic defect, identified in more than 50% of T-cell acute lymphoblastic leukemias, supporting an essential driver role for this gene in T-cell acute lymphoblastic leukemia oncogenesis. In this study, the transcriptional response of all protein coding genes and long non-coding RNAs upon pharmacological Notch inhibition in the human T-cell acute lymphoblastic leukemia cell line CUTLL1 was measured using RNA-sequencing

## Overall design:

CD34+ cells of 2 healthy donors are cultured on a OP9-GFP or OP9-DLL1 feeder layer. OP9-DLL1 activates NOTCH1.

![Experimental design](https://github.com/ulisesplaper/transcripts_under_NOTCH1/blob/master/data/experimental_design.png?raw=true)

## Script summary

### 1. recount3.R

Download RSE object, make exploratory analysis and save it in the data directory.

### 2. diffe_exp_analysis.R

Make differential expression analysis and generate plots.
