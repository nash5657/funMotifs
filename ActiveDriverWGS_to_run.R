library(ActiveDriverWGS)
library(tidyverse)
library(data.table)
#args
args = commandArgs(trailingOnly=TRUE)
#elements file
df_elem = args[1]
#mutations file
df_mut = args[2]
#min number of mutations in element
n_mut_ele = args[3]
#output file
out_file = args[4]
#output file  sig
out_file_sig = args[5]
#output file results
out_file_active = args[6]
#number of cores
n_cores = args[7]

#print(df_elem)
#print(n_mut_ele)
#elements
data_my <- read.csv(df_elem, sep = '\t', header = FALSE)
#head(data_my) %>%  print

#at least 2 mutations in an element
data_my_2mut <- data_my[which(data_my[,6] >= as.numeric(n_mut_ele)),] 

#elements to test
data_elem <- data_my_2mut[,c(1:3,15)] 
colnames(data_elem) <- c('chr',
                         'start',
                         'end',
                         'id')
data_elem$chr <- gsub('23','X',data_elem$chr)
data_elem$chr <- gsub('24','Y',data_elem$chr)  
data_elem$chr <- paste0('chr',data_elem$chr)
data_elem$start <- as.numeric(data_elem$start)
data_elem$end <- as.numeric(data_elem$end)
data_elem$id <- as.character(data_elem$id)


#mutations to test
data_mut <- fread(df_mut, header = FALSE)
mut <- data_mut
colnames(mut)[1:9] <- c('chr', 'pos1', 'pos2','ref', 'alt', 'cohort','var_type', 'sample_id','patient')
mut_in <- mut[,c('chr', 'pos1', 'pos2','ref', 'alt', 'patient')] %>%  unique
mut_in$chr <- gsub('23','X',mut_in$chr)
mut_in$chr <- gsub('24','Y',mut_in$chr)
mut_in$pos1 <- as.numeric(mut_in$pos1)
mut_in$pos2 <- as.numeric(mut_in$pos2)
if(dim(data_my_2mut)[1] >0){
  result = ActiveDriverWGS(mutations = mut_in,
                           elements = data_elem, mc.cores = as.numeric(n_cores))
  colnames(result)[1] = 'id'
  colnames(data_my)[15] = 'id'
  result_join <- left_join(data_my,result, by= 'id' ) %>%  dplyr::select(-c(element_muts_obs:site_enriched, fdr_site:has_site_mutations))
  no_pval_id<=which(is.na(result_join$pp_element))
  result_join$pp_element[no_pval_id]=1.0
  result_join$fdr_element[no_pval_id]=1.0
  result_join_sig=result_join[which(result_join$fdr_element<0.05),]
}else{
  result_join = NULL
  result_join_sig=NULL
  
} 

write.table(result, file = out_file_active, sep='\t', col.names = FALSE, row.names = FALSE, quote = FALSE, na = 'NA')

#colnames(result_join)[(length(result_join[1,])-1):length(result_join[1,])] <- c('ElementPval_ActiveDriver',   'ElementFDR_ActiveDriver')
write.table(result_join, file = out_file,sep='\t', col.names = FALSE, row.names = FALSE, quote = FALSE, na = 'NA')
write.table(result_join_sig, file = out_file_sig,sep='\t', col.names = FALSE, row.names = FALSE, quote = FALSE, na = 'NA')
