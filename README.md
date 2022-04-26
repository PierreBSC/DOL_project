# Analysis of Disease associated Oligodendrocytes (DOLs) 

This repository contains the code and scripts used for the single-cell/spatial transcriptomic analysis of mouse oligodendrocytes in Alzheimers disease.
The data used in this project can be accessed using the following link :

- MARS-seq 2.0 data generated from whole WT and 5XAD mice brain : GSE....[]
- Chromium 10X data from the hippocampi of TauP301L, PS2APP/TauP301L, and PS2APP/TauP301L/TREM2KO mice (Lee et al. 2021) : [GSE153895](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153895).
- Smart-seq2 data from healthy and EAE mice (Falcão et al. 2018) : [GSE113973](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113973)
- Chromium 10X data from healthy and EAE mice at various stages of the disease (Wheeler et al. 2020) : [GSE129609](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129609).
- Chromium 10X data from Aging SVZ (Dulken et al. 2019) : [PRJNA450425](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA450425)
- Visium® 10X data from LPS injected mouse brains (Lee et al. 2021) : [GSE165098](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165098).
- Visium® 10X data from Non Dementia Control (NDC) and Alzheimer's Disease (AD) patients [(Chen et al. 2021)](https://www.biorxiv.org/content/10.1101/2021.07.07.451554v1.full): data are available upon request to the authors.

Additionnaly, two small files are needed for the analysis :
- DOL_genes.txt : a txt file containing the list of the DOL genes 
- tissue_positions_list.csv : a csv file with the correspondanc between Visium spot location and nucleotide barcode. It is needed for the analysis of the two Visium® datasets.

All those dataset were analysed using the [**pagoda**](https://github.com/kharchenkolab/pagoda2) pipeline. To install it please go on the corresponding [Github](https://github.com/kharchenkolab/pagoda2) page. Other packages need to be installed :
- [fifer](https://github.com/dustinfife/fifer)
- uwot, RColorBrewer, pheatmap and liger (CRAN packages)
- DESeq2, GEOquery and CountClust (Bioconductor packages)

This repository contains the following scripts :

- EAE_Falcao_analysis_script.R : script for the analysis of the Smart-Seq2 EAE dataset.
- EAE_Wheeler_analysis_script.R : script for the analysis of the 10X EAE dataset.
- Aging_SVZ_analysis_script.R : script for the analysis of the 10X aging SVZ analysis.
- 
