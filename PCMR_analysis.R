# 加载必要的包
library(dplyr)
library(stringr)
library(rstatix)
library(tidyr)

omics = "microRNA"
omics_lower_case = tolower(omics)

gseID = "GSE20099"


#### 输入文件目录
dir_input = paste0("D:/pcmr/20240920/", omics,"/GSE10714/")

#### 输出文件目录
dir_output = paste0(dir_input,"/res/", gseID)
if(!file.exists(dir_output)){
  dir.create(dir_output, recursive = T)
}



#### 注释文件
## microRNA
mir_ref_bed = read.table( paste0("F:/project/database_pub/ref/","ref_miR.bed"),sep="\t",header = T)
colnames(mir_ref_bed) = c("chromosome","start","end","gene_symbol","strand","mature_id","mature_acc","precursor_id")

##mRNA
mir_ref_bed = read.table( paste0("F:/project/database_pub/ref/","ref_mRNA.bed"),sep="\t",header = T)
colnames(mir_ref_bed) = c("gene_symbol","chromosome","start","end","strand")


## GPL ID mapping
GPL_anno = read.table( paste0(dir_input,"GPL24158.txt"),sep="\t",header = T, row.names = 1 )


#### 临床信息表的表头：
#### organ_system	primary_site	cancer_type	disease_category	label	tissue	omics	omcis_sub	omcis_platform	platform	ncbi_gpl_id	ncbi_gse_id	is_raw_data	ncbi_gsm_id	geo_info	gender	age	cancer_stage	pmid
f_clinical = paste0(dir_input, "/clinical_info.txt")


#### 表达矩阵，行为RNA，列为样本
f_exp = paste0(dir_input, "/GSE211692_processed_data_matrix.txt")





########################################################################## 数据准备

ESCC_EV_clinical=read.table(f_clinical, header = T, #row.names = 2,
                            sep = "\t", check.names = F)
#rownames(ESCC_EV_clinical) = ESCC_EV_clinical$geo_info


#### 表达矩阵，行为RNA，列为样本
exp_total_sample <- read.table(f_exp, header = T, row.names = 1,
                               sep = "\t", check.names = F)


#### 提取有配对premalignant的样本，包括癌症，癌前，正常
dim(exp_total_sample)
list_all_genes = as.data.frame( t( exp_total_sample[, ESCC_EV_clinical$ncbi_gsm_id ] ) )
dim(list_all_genes)


#### 行名修改为GSM ID （上一步中进行了矩阵转置）
#rownames(list_all_genes) = ESCC_EV_clinical[rownames(list_all_genes), "ncbi_gsm_id"]
list_all_genes = round(list_all_genes, digits = 3)


#### problem1: 有的列包括多个miR id, 需要拆开，每一个miR赋值相同
dim(list_all_genes)
final_list =  as.data.frame( rep(as.data.frame(list_all_genes), str_count(colnames(list_all_genes), ",")+1))
colnames(final_list) = gsub(" ","", unlist(strsplit(as.character(colnames(list_all_genes)), ",")) )
rownames(final_list) = rownames(list_all_genes)
dim(final_list)


if(!is.numeric(final_list[1,1])){
  print("!!!!!!!!!!!!!!!!!!!!!!!!!!!! ERROR, data is not numeric type")
}



##########################################################################
#### 数据的标准化.
#### normalizeBetweenArrays输入的矩阵：行是基因，列是样本

library(limma)

final_list = t(final_list)
final_list <- normalizeBetweenArrays(final_list)
final_list = t(final_list)

###########################################################################




##################################################################################### 

final_list = final_list[ ,colSums(is.na(final_list))<1]
# final_list_stats = final_list[,1:5]
final_list_stats = as.data.frame( final_list )

max_frequencies <- numeric(ncol(final_list_stats))
# 遍历每一列，计算最大频数
for (i in seq_along(final_list_stats)) {
  freq <- table(as.character(final_list_stats[[i]]))
  max_frequencies[i] <- max(freq)
}


final_list_stats = final_list_stats[, max_frequencies/nrow(final_list_stats)<0.9]



###################### calculate p.val and p.adj
df_for_ggplot2 <- as.data.frame(as.table(as.matrix(final_list_stats)))
df_for_ggplot2$Var1 = as.character(df_for_ggplot2$Var1)
df_for_ggplot2$Var2 = as.character(df_for_ggplot2$Var2)

rownames(ESCC_EV_clinical) = ESCC_EV_clinical$ncbi_gsm_id
df_for_ggplot2$label = ESCC_EV_clinical[df_for_ggplot2$Var1, "label"]
df_for_ggplot2$cancer_type = ESCC_EV_clinical[df_for_ggplot2$Var1, "cancer_type"]

colnames(df_for_ggplot2) = c("gsm","GeneSymbol","exp","label","cancer_type")


#####=============================================================================================================
##### 如果一套数据集中有normal，disease1,disease2,disease3，将normal分别复制一份给对应的疾病
##### 临床数据表格中，normal样本对应的cancer_type 这一列为各种疾病名称用逗号分隔。如果无多种疾病，该数据集cancer_type这一列都是没有逗号的
all_cancer_type = unique(df_for_ggplot2$cancer_type)
sep_disease_name = ";"

if(sum(grepl(sep_disease_name, all_cancer_type ) )>0){
  all_cancer_type = all_cancer_type[!grepl(sep_disease_name, all_cancer_type)]
  
  normal = df_for_ggplot2[grepl(sep_disease_name, df_for_ggplot2$cancer_type), ]
  for(one_cancer_type in all_cancer_type  ){
    one_normal = normal
    one_normal$cancer_type = one_cancer_type
    df_for_ggplot2 = rbind(df_for_ggplot2, one_normal)
  }
  df_for_ggplot2 = df_for_ggplot2[!grepl(sep_disease_name, df_for_ggplot2$cancer_type), ]
}
#####=============================================================================================================


table(df_for_ggplot2$cancer_type, df_for_ggplot2$label)


df_p_val0 <- df_for_ggplot2 %>% group_by(cancer_type, GeneSymbol) %>%
  wilcox_test(exp  ~ label) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% as.data.frame



############################## 统计Minimum,First quartile,Median,Third quartile,Maximum
############################## 计算 p-val: HC vs pre, pre vs cancer, HC vs cancer
final_list_stats = as.data.frame(final_list)
final_list_stats = final_list_stats[, colSums(is.na(final_list_stats))<1]


rownames(ESCC_EV_clinical) = ESCC_EV_clinical$ncbi_gsm_id
final_list_stats = as.data.frame(final_list_stats)
final_list_stats$label = ESCC_EV_clinical[rownames(final_list_stats), "label"]
final_list_stats$cancer_type = ESCC_EV_clinical[rownames(final_list_stats), "cancer_type"]


# testData = final_list_stats[,c("hsa-miR-28-3p", "hsa-miR-27a-5p", "hsa-miR-518b", "cancer_type", "label")]
dat_boxplot = aggregate(final_list_stats[,1:(ncol(final_list_stats)-2)], by=list(cancer_type=final_list_stats$cancer_type, label=final_list_stats$label), summary)




