**Important note:** xCell was trained for use with bulk gene expression data. Due to many reasons we believe it is not a great solution for using it directly to infer cell types in scRNA-seq data. For such data we recommed using our method SingleR which was developed exactly for this purpose - https://github.com/dviraran/SingleR 

**User guide:** For more details on using xCell, please refer to Aran D. (2020) Cell-Type Enrichment Analysis of Bulk Transcriptomes Using xCell. In: Boegel S. (eds) Bioinformatics for Cancer Immunotherapy. Methods in Molecular Biology, vol 2120. https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_19 

# xCell - cell types enrichment analysis

xCell is a webtool that performs cell type enrichment analysis from gene expression data for 64 immune and stroma cell types. 

xCell is a gene signatures-based method learned from thousands of pure cell types from various sources. xCell applies a novel technique for reducing associations between closely related cell types. xCell signatures were validated using extensive in-silico simulations and also cytometry immunophenotyping, and were shown to outperform previous methods. xCell allows researchers to reliably portray the cellular heterogeneity landscape of tissue expression profiles. 

For more informations please refer to the xCell manuscript.

# Install

```R
devtools::install_github('dviraran/xCell')
```

Note: if the installation fails, try setting `options(unzip = "internal")` before calling the install function (thanks @hermidalc).

# Usage

```R
library(xCell)
exprMatrix = read.table(expr_file,header=TRUE,row.names=1, as.is=TRUE)
xCellAnalysis(exprMatrix)
```

This function performs all three steps in xCell, which can be performed seperately as well:

1. rawEnrichmentAnalysis
2. transformScores
3. spillOver

xCell loads the xCell.data object which is a list with the spill over and calibration parameters, the signatures and the list of genes it uses. However, one can use different signatures and different spill over functions to perform the analysis.

In addition xCell provides significance assessment of the null hypothesis that the cell type is not in the mixture:

1. xCellSignifcanceBetaDist - uses predefined beta distribution parameters from random mixtures generated using the reference data sets. Recommended. 
2. xCellSignifcanceRandomMatrix - uses random matrix and caculates beta distribution parameters. 

See documentation for more details on each function.

# Data input

The expression matrix should be a matrix with genes in rows and samples in columns. The rownames should be gene symbols. xCell uses the expression levels ranking and not the actual values, thus normalization does not have an effect, however normalizing to gene length is required. 

Importantly, xCell performs best with heterogenous dataset. Thus it is recommended to use all data combined in one run, and not break down to pieces (especially not cases and control in different runs).

# Notes for correct usage

1. xCell produces enrichment scores, not percentages. It is not a deconvolution method, but an enrichment method. That means that the main usage is for comparing across samples, not across cell types. xCell does an attempt to make the scores resemble percentages, but it is a hard problem, and is very platform and experiment specific. We have made some tests to compare the ability of xCell for cross-cell types analysis, and found that it generally performed better in that than other methods (on limited and comparable cell types), but this type of analysis should be performed carefully.

Regarding this issue, scaling the scores by samples is extremely dangerous and will inevitably will result in false interpretations.

2. For the linear transformation xCell uses a calibration parameter to make it resemble percentages. If your analysis produces high scores for a cell type that is clearly false, you may tweak the calibration parameters, which can be found in the xCell.data object (xCell.data$spill$fv$V3). In the paper, we used an automatic method to learn these parameters, but if you have prior knowledge about your mixture, this is very handy to get better results. 

3. xCell uses the variability among the samples for the linear transformation. xCell will only function with heterogenous mixtures. If there is no variability between the samples, xCell will not identify any signal. As noted above, it is highly recommended to use all data combined in one run. Failing to do so will again inevitably make xCell's results false. 

4. xCell is a method for detecting cell types enrichments from mixed samples, not to detect the cell of origin. xCell is still able to recognize the cell of origin many times, but it was not trained for this problem, and there is no expectation that it will perform better than other methods for that problem. In accordance, it is strongly not advised to be performed on single-cell data. 

# Vignettes

**11.6.2018 update:** Due to popular demand, we created a simple vignette to reproduce the analysis that was performed in the xCell paper of [correlating the xCell scores with the CyTOF immunoprofiling data available from ImmPort](http://comphealth.ucsf.edu/sample-apps/xCellView/xCell-Immport.html). The data files are available in the 'vignette' directory.

# Contributors

xCell is developed by the Butte lab. Please contact Dvir Aran: dvir.aran at technion.ac.il for any questions or suggestions.

To cite xCell: Aran, Hu and Butte, xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biology (2017) 18:220
