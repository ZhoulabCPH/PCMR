# PCMR
PCMR, A Comprehensive PreCancerous Molecular Repository and Online Analysis Platform (http://www.bio-data.cn/pcmr)
1. PCMR_analysis.R was used to calculate adjusted P-value between each paired groups. 
  - The P-value for comparing two groups in differential analysis was calculated using the Mann-Whitney U test and adjusted with the Bonferroni correction.
   - In the differential DNA methylation analysis, genes with an absolute methylation difference greater than 0.2 and an adjusted P-value less than 0.05 were considered differentially methylated.
   - For mRNA, microRNA, circRNA, and protein analysis, a fold change greater than 2 and an adjusted P-value below 0.05 were considered as differential.
2. This code will also compute the minimum, first quantile, median, mean, third quantile, and maximum for each group.
