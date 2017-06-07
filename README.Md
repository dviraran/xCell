# xCell - cell types enrichment analysis

xCell is a webtool that performs cell type enrichment analysis from gene expression data for 64 immune and stroma cell types. 

xCell is a gene signatures-based method learned from thousands of pure cell types from various sources. xCell applies a novel technique for reducing associations between closley related cell types. xCell signatures were validated using extensive in-silico simulations and also cytometry immunophenotyping, and were shown to outperform previous methods. xCell allows researchers to reliably portray the cellular heterogeneity landscape of tissue expression profiles. 

For more informations please refer to the xCell manuscript.

# Usage

library(xCell)

exprMatrix = read.table(expr_file,header=TRUE,row.names=1, as.is=TRUE)

xCellAnalysis(exprMatrix)

# Data input

The expression matrix should be a matrix with genes in rows and samples in columns. The rownames should be gene symbols. xCell uses the expression levels ranking and not the actual values, thus normalization does not have an effect, however normalizing to gene length is required. 

Importantly, xCell performs best with heterogenous dataset. Thus it is recommended to use all data combined in one run, and not break down to pieces (especially not cases and control in different runs).

# Install

devtools::install_github('dviraran/xCell')

# Contributors

xCell is developed by the Butte lab. Please contact Dvir Aran: dvir.aran at ucsf edu for any questions or suggestions.

To cite xCell: Aran, Hu and Butte, xCell: Digitally portraying the tissue cellular heterogeneity landscape. bioRxiv, 2017
