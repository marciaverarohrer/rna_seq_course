**PCA-Plot:** 
This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.

### **Axes meanings**

- **PC1 (x-axis):** The direction of greatest variance between samples
- **PC2 (y-axis):** The second-most variation direction
- The **percentages in parentheses**, e.g. _PC1: 37%_, mean:
 PC1 explains 37% of the total variance in gene expression.
### **How to interpret**

- Samples that cluster **closely together** → have **similar expression profiles**
- Samples **far apart** → large biological or technical differences
- Ideally, samples of the **same experimental group** should cluster together  
    (WT cases together, WT controls together, etc.)
If they don’t, it may indicate:
- Batch effects
- Sample mix-ups
- Strong outliers
- True biological heterogeneity

Q&A 11.12

condition one, condition 2
feature counts: how many reads do we have per gene?
gene 1 cond. 1 , 10 reads, cond. 2, 20 reads

normalization: how many counts per million of sequenced reads
CPM is counts per million
nr. of reads can be asked for in the sequencing facility

variance stabilizing t (vst) , we fit sample trying to stabilize what is already known, variability of data should be accounted for. 

importance of replication : how many samples do I have per condition?
biological replicates are needed to make more reliable conclusions.
more replicates -> more statistical power

counts is number of reads overlapping in a pathway

one circle is a pathway gene ontology with 5 genes, another circle DEG 200, overlap is number of 
more counts is more confortability that it is actually significant

gene or reads
