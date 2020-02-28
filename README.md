# Research_GRN
Gene regulatory networks for plant embryogenesis using WGCNA

This repository is based on my undergraduate thesis "Uncovering Gene Regulatory Networks Controlling Embryogenesis"

Temporal transcriptomic data were used to predict interactions among the genes involved in zygotic embryogenesis (ZE) and somatic embryogenesis (SE) of Arabidopsis using Weighted Gene Co-expression Network Analysis (WGCNA) package in R (Langfelder and Horvath, 2008).

The "data" folder contains the following datasets used for the study;
- Microarray expression data of SE [Becker et al. (2014)]
- RNA-Seq expression data of ZE [Gao et al. (2019)]
- Arabidopsis TFs [Plant Transcription Factor Database (http://planttfdb.cbi.pku.edu.cn/index.php)] 
- Embryo defective (EMB) genes [Meinke (2019)]
- SE marker genes [Magnani et al. (2017)]
- Genes related to hormones and enzymes expressed during embryogenesis [Xiang et al. (2011)]
- Arabidopsis gene descriptions [The Arabidopsis Information Resource (TAIR), 2019]

For each dataset of SE and ZE, the following codes are included;
- Preprocessing (Gene ID conversion, Gene filtering)
- Network Construction (Module detection, Moudle-trait relationship analysis, Validation)
- Analysis and Visualization (Functional analysis, ggplot2 visualizations)
- Module analysis (Module network construction)
