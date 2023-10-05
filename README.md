# Single-cell hematopoietic classification and feature engineering using single-cell DREMI scores 

In this study we used single-cell mass cytometry to characterise normal and leukemic hematopoiesis, using data from leukemia patients and  healthy donors. To annotate reliably single-cells, we profiled healthy donor peripheral blood and bone marrow, and we constructed a reference of hematopoietic cell types. Linear Discriminant Analysis was used to develop a single-cell hematopoietic classifier, and via Transfer Learning single-cells from patients with leukemia were mathced to their nearest healthy cell counterparts. Our analysis characterised the cellular composition of leukemia patients at time of diagnosis, while Machine Learning modeling with Feature Selection identified cell type-specific interactions between signalling rproteins that were associated with survival. 

We provide a bioinformatics pipeline that allows users to:
- implement a gating strategy based on density estimations
- develop a hematopoietic classifier using Machine Learning
- use the hematopoietic classifier to annotate cells from leukemia patients   
- engineer cell type-specific features using DREMI that allow for Machine Learning modelling

## Publication

Title: Identifying predictors of survival in patients with leukemia using single-cell mass cytometry and machine learning

Journal: The paper is submitted to Cell Reports Methods

Published: pre-print available at bioRxiv https://doi.org/10.1101/2022.08.13.503587

## Dependencies and System Requirements

Our bioinformatics pipeline is developed with R 4.1.3 (2022-03-10), running under: macOS Big Sur 11.6.3 

The analysis presented in our publication requires packages to be installed, and we refer to the following document about the R session info.

```
mySessionConfig.sh
```

Please refer to Execution_examples.md for more information. We provide step by step execution guidelines in the form of a tutorial.

More detailed vignette will appear soon. Please note that the examples presented here are minimal. 


#### Current Release

16-Mar-2022 : Beta version 1

## Contact

Comments and bug reports are welcome, please email: Dimitrios Kleftogiannis (dimitrios.kleftogiannis@uib.no)

We are also interested to know about how you have used our source code, including any improvements that you have implemented.
 
You are free to modify, extend or distribute our source code, as long as our copyright notice remains unchanged and included in its entirety. 

## License

This project is licensed under the MIT License.

Copyright 2022 Department of Informatics, University of Bergen (UiB) and the Centre of Cancer Biomarkers (CCBIO), Norway

You may only use the source code in this repository in compliance with the license provided in this repository. For more details, please refer to the file named "LICENSE.md".
