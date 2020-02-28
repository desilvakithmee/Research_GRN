# Research_GRN
Gene regulatory networks using WGCNA

This repository is based on my undergraduate thesis "Uncovering Gene Regulatory Networks Controlling Embryogenesis"

Temporal transcriptomic data were used to predict interactions among the genes involved in Arabidopsis ZE and SE using Weighted Gene Co-expression Network Analysis (WGCNA) package in R (Langfelder and Horvath, 2008).

The "data" folder contain the expression datasets used for the study and the following additional datasets;
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
