# The ensembling method associates prepregnant BMI with metabolomic changes in a multi-ethnic cohort in Hawaii

## Description

This is the github repository for the project **The ensembling method associates prepregnant BMI with metabolomic changes in a multi-ethnic cohort in Hawaii** by Leyang Tao, Bowei Li, Yuheng Du, Shengyu Hung and Lana Garmire. This project provides a comprehensive pipeline for analyzing metabolomic profiles in relation to maternal obesity using continuous BMI values. It incorporates a stacking ensemble machine learning model that combines linear and nonlinear approaches to elucidate ethnicity-specific metabolic variations. It contains code and data for generating Figure 1,2 in the paper. 

## Getting Started

### Dependencies
* Linux Working Environment
* [R](https://www.R-project.org)
* [Anaconda3](https://www.anaconda.com/)
* [Jupyter](https://jupyter.org)
* Python libraries:
  * [Numpy](https://numpy.org/)
  * [pandas](https://pandas.pydata.org/docs/index.html)
  * [Scipy](https://scipy.org/)
  * [category_encoders](https://pypi.org/project/category-encoders/)
  * [Scikit-learn](http://scikit-learn.org/)
  * [matplotlib](https://matplotlib.org/)
  * [seanborn](https://seaborn.pydata.org/)
  * [shap](https://shap.readthedocs.io/en/latest/)
  * [xgboost](https://pypi.org/project/xgboost/)
* R libraries:
  * [bigsnpr](https://privefl.github.io/bigsnpr/)
  * [brglm2](https://ikosmidis.r-universe.dev/brglm2/)
  * [car](https://github.com/cran/car/)
  * [caret](https://topepo.github.io/caret/)
  * [corrplot](https://cran.r-project.org/package=corrplot)
  * [data.table](https://rdatatable.gitlab.io/data.table/)
  * [dummies](https://cran.r-project.org/package=dummies)
  * [effects](https://cran.r-project.org/package=effects)
  * [gdata](https://cran.r-project.org/package=gdata)
  * [dplyr](https://dplyr.tidyverse.org/)
  * [ggplot2](https://ggplot2.tidyverse.org/)
  * [ggrepel](https://cran.r-project.org/package=ggrepel)
  * [Hmisc](https://cran.r-project.org/package=Hmisc)
  * [impute](https://www.bioconductor.org/packages/release/bioc/html/impute.html)
  * [mice](https://cran.r-project.org/package=mice)
  * [limma](http://bioconductor.org/packages/release/bioc/html/limma.html)
  * [pathifier](https://cran.r-project.org/package=pathifier)
  * [pathview](https://www.bioconductor.org/packages/release/bioc/html/pathview.html)
  * [princurve](https://cran.r-project.org/package=princurve)
  * [reshape](https://cran.r-project.org/package=reshape)
  * [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/00Index.html)  <!-- Note: base package, link provided for reference -->
  * [stringr](https://stringr.tidyverse.org/)
  * [sva](https://bioconductor.org/packages/release/bioc/html/sva.html)
  * [tidyverse](https://www.tidyverse.org/)
  * [vsn](https://bioconductor.org/packages/release/bioc/html/vsn.html)

### Installing

Installing the R kernel on the jupyter
```R
install.packages('IRkernel')
IRkernel::installspec()  # to register the kernel in the current R installation
```
Use the [Bioconductor](https://www.bioconductor.org/install/) to install R packages.
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(*PackageName* = "*Version*")
```




Use the package manager [pip](https://pip.pypa.io/en/stable/) to install python packages.
```bash
pip install Numpy
```




## Repository directories & files

The directories are as follows:
+ [`Figures`](Figures) contains subfigures of Figure 1 to 6, excluding 4th, showing in the paper.
+ [`data`](data) empty folder to put data that required for preprocessing and further figure generation.
+ [`output`](output) Contains output from each script. 
+ [`lilikoipack`](lilikoipack) contains the R script of [lilikoipack](https://github.com/lanagarmire/lilikoi) for [`Step3_Pathway_Analysis.R`](Step3_Pathway_Analysis.R)

The other files are as follows.
+ [`Step1_Preprocessing.R`](Step1_Preprocessing.R) contains the preprocessing code for metabolites data and SOV figure generation. 
+ [`Step2_SHAP.ipynb`](Step2_SHAP.ipynb) contains the  the primary meta-learner and ensemble approach applied to the processed metabolites data, with results compared against a single linear regression model and interpreted using SHAP values.
+ [`Step3_Pathway_Analysis.R`](Step3_Pathway_Analysis.R) contains the pathway analysis with Kyoto Encyclopedia of Genes and Genomes (KEGG) database. 
+ [`Step4_Correlation.R`](Step4_Correlation.R) contains spearman correlation analysis to find out BMI associated metabolites.
+ [`Step5_Figure2BtoE.R`](Step5_Figure2BtoE.R) generates bar plots comparing metabolite levels across ethnicities, annotated with Wilcoxon p-values.


### Local execution
+ for `.ipynb` files: Using jupyter lab to execute the codes


## Data Availability
The data used in this study (supposedly located in the `data` folder) is available upon reasonable request. Metabolites data is publicly available in metabolomics workbench repository (https://www.metabolomicsworkbench.org/., study ID ST001114).


## Authors

### Contributors names and contact info

+ [@Lana Garmire](https://github.com/lanagarmire)
+ [@Leyang Tao](https://github.com/leyangt)
+ [@Bowei Li](https://github.com/GHBLA-NI)
+ [@Yuheng Du](https://github.com/yhdu36)
+ [@Shenyu Hung](https://github.com/Shengyu011314)

### Current Maintainer
* Bowei Li - https://github.com/GHBLA-NI

## License

This project is licensed under the `GNU General Public License v3.0` License - see the LICENSE.md file for details
