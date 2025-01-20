 # Predictive Modeling and Network Analysis of Gene Expression in Breast Cancer

This repository contains the project for the Statistical Methods for High Dimensional Data course of the Master's Degree in Data Science at the University of Padua.

## Authors

- Eleonora Mesaglio
- Margherita Rigato
- Francesco Lollato


## Overview

The *dataset* was sourced from the [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70947) and contains:

- *289 samples* (143 cancerous and 146 normal tissues).
- *35,982 genes* with quantile-normalized expression values.
- The dataset is balanced with no missing values.

Our primary objective is to classify tissues as cancerous or normal based on gene expression data and analyze the network structure of gene relationships in both conditions.


## Methodology:

### Predictive Models

Several machine learning models were implemented to classify tissue samples:
1. *Lasso Logistic Regression*
2. *Support Vector Machines (SVM)* with Lasso regularization:
   - Squared Hinge Loss.
   - Hinge Loss.
3. *Elastic Net Logistic Regression*
4. *Random Forest*

### Network Analysis

We built gene interaction graphs to explore differences in correlation and dependency structures between normal and cancerous tissues:
1. *Correlation Networks*
2. *Graphical Lasso*

## Tools Used
We used R for all the analysis and Overleaf to write the report.

## How to Run the Code

1. Clone this repository:
   bash
   git clone https://github.com/your-repo-name.git
   
2. Install required R packages:
   R
   install.packages(c("glmnet", "gcdnet", "sparseSVM", "ranger", "igraph", "corrplot", "tidyverse", "broom", "Matrix", "randomcoloR"))
   
3. Download  the dataset from this link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70947.
4. Place the file (Breast_GSE70947.csv) in the working directory.
5. Open and run the R script codice.R in an R environment or RStudio.

## License
This project is licensed under the terms of the MIT license.
