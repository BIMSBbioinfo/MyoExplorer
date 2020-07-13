# MyoExplorer
Companion code for Kim et al., Single-nucleus transcriptomics reveals functional compartmentalization in syncytial skeletal muscle cells
bioArxiv doi: https://doi.org/10.1101/2020.04.14.041665

This repository contains the companion code for recreating the single nucleus analysis figures from the Kim et al.

The markdown document contains all the code required to reproduce the main figures from the manuscript.

To run the markdown document you will need **R4.0**, and the following libraries:

    SingleCellExperiment_1.10.1
    ggplot2_3.3.2
    dplyr_0.8.5
    data.table_1.21.8
    ComplexHeatmap_2.4.2
    stringr_1.4.0
    tidyr_1.1.0
    magrittr_1.5


The MyoExplorer_Figures_Code.rmd markdown document contains the code to reproduce all computational figures.
It contains a section which installs the required libraries, and downloads the data into a prespecified folder.

Before running the code please specify the **params: input_file_location** parameter, in the header
of the .rmd document. Input files will be downloaded into the specified folder.



