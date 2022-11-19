Cancer Genomics
Download Bulk dataset from GEO dataset where I selected GSE130817 for leukemia cancer

Steps:
Normalize the dataset
Produce a Heatmap
Prepare a Differential Expression data from the main dataset
Plot the CircosPlot of those Upregulated Differntial Genes
Create SSGSEA function
Use the Logcpm data and plot the SSGSEA graph
Read the log2 normalized data into an object for further analysis
zscore the ssgsea output for comparative analysis
Heatmap plot for zscore

For Survival Analysis use another dataset and visualize the survival graphs

Single Cell Rna Seq:
Download the Singel Cell leukemia dataset from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212038	
Steps: Attaching the "dplyr," "Seurat," and "patchwork" packages

Load the dataset for prostate cancer


start the Seurat object with the unprocessed (non-normalized data).
using the Read10X() method, read a file
Plot QC metrics using a violin diagram.
FeatureScatter: FeatureScatter can be used for everything the object calculates, such as columns in object metadata, PC scores, etc. It is commonly used to illustrate feature-feature correlations.
Plotting labelled and unlabeled variations of variable features
graphing the nearest neighbour and compute SNN
Create clusters for the features
Plot the featureplots.

Interpretation:

I have created numerous popular plots for displaying single-cell data, including heatmaps and t-SNE plots (using dimensionality reduction technologies). These algorithms seek to understand the underlying data manifold in order to place
Using this I was able to identify the expressible genes.
