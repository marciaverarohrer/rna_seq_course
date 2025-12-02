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