#Sleuth Analysis in R
library('biomaRt')
library('RCurl')
library('gtools')

setwd('/Users/HoltLab/Documents/R')
preparedata = function(base_dir, Ids, compare){
  library('biomaRt')
  library('RCurl')
  setwd('/Users/HoltLab/Documents/R')
  options(mc.cores = 6L) # use the 6 cores of the mac pro
  #prepare for loading
  base_dir =  base_dir  #'/Users/HoltLab/Documents/R/RNA_seq_ignacio_11.28'
  sample_id =dir(base_dir)
  
  kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, "kallisto"))
  kal_dirs
  
  sample_condition = t(data.frame('1WTDex','2SpidrinDex','3H2ADex','4deltaQDex', '5NullDex','6WTStarve','7SpidrinStarve','8H2AStarve','9deltaQStarve','9NullStarve','1WTDex','2SpidrinDex','3H2ADex','4deltaQDex','5NullDex','6WTStarve','7SpidrinStarve','8H2AStarve','9deltaQStarve','9NullStarve'))
    #(data.frame(c(rep('Dex',5), rep('Starve',5),rep('Dex',5), rep('Starve',5) ))) #for comparing WT DEX
  rownames(sample_condition) = 1:20
  
  print(Ids)
  s2c = data.frame(sample_id,sample_condition)
  colnames(s2c) = c("sample","condition")
  
  s2c <- dplyr::mutate(s2c, path = kal_dirs)
  growthcond = c(rep('Dex',5), rep('Starve',5),rep('Dex',5), rep('Starve',5) )
  s2c = cbind(s2c, growthcond)
  s2c_cut = s2c[grep(Ids, s2c[,2]),]
  
  s2c_cut$condition = factor(s2c_cut$condition)
  s2c_cut$sample = factor(s2c_cut$sample)
  
  
  # add annotations
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "scerevisiae_gene_ensembl",
                           host = 'ensembl.org')
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  
  #so <- sleuth_prep(s2c_cut, ~condition, target_mapping = t2g, aggregation_column = 'ens_gene')
  # to run it in gene based analysis mode
  
  #prep/fit
  so <- sleuth_prep(s2c_cut, ~ condition , target_mapping = t2g)
  so <- sleuth_fit(so)
  #return(so)
 # so
  
  #perform differential expression analysis
  so <- sleuth_wt(so, compare)
  return(so) # return sleuth object when not directly writing files
  
  ### write results 
  results_table <- sleuth_results(so, compare)
  results_ordered <- results_table[order(results_table$qval),]
  table(results_ordered$qval <= 0.001)

  # Writes tables for dif expressed genes
  write.table( subset(results_ordered, qval <= 0.001), file= paste('DE',Ids,'.qval_0.001.txt', sep = ''), sep="\t",row.names=F, quote=F)
  #write table with b-score less than -1
  write.table( subset(results_ordered, b <= -1), file=paste('DE',Ids,'b_-1.txt', sep = ''), sep="\t",row.names=F, quote=F)
  print('results written to working directory')
  #write table with b-score greater than 1
  write.table( subset(results_ordered, b >= 1), file=paste('DE',Ids,'b_1.txt', sep = ''), sep="\t",row.names=F, quote=F)
 
  ###### MA plotter
  ma = plot_ma(so, compare, test_type = "wt", which_model = "full",
          sig_level = 0.0001, point_alpha = 0.4, sig_color = "red",
          highlight = NULL, highlight_color = "green")
  pdf(paste('MA',Ids,'.pdf',sep = ""))
  print(ma)
  dev.off()
  ###### Volcano plotter
  volcano = plot_volcano(so, compare, test_type = "wt", which_model = "full",
                         sig_level = 0.0001, point_alpha = 0.4, sig_color = "red",
                         highlight = NULL) #not working right now
  
  pdf(paste('Volcano',Ids,'.pdf',sep = ""))
  print(volcano)
  dev.off()
 
}

so = preparedata('/Users/HoltLab/Documents/R/RNA_seq_ignacio_11.28','1WTDex|6WTStarve', 'condition6WTStarve')

#for running comparisions withing Dex and Starve groups
A = c('1WTDex','2SpidrinDex','3H2ADex','4deltaQDex', '5NullDex','6WTStarve','7SpidrinStarve','8H2AStarve','9deltaQStarve','9NullStarve','1WTDex','2SpidrinDex','3H2ADex','4deltaQDex','5NullDex','6WTStarve','7SpidrinStarve','8H2AStarve','9deltaQStarve','9NullStarve')
#least ratio test
a = combinations(n=5,r=2)
c = c(5,5)
D = rbind(c,c,c,c,c,c,c,c,c,c)
rownames(D) = 1:10
e = D + a

for(i in 1:length(e[,1])){
  B = A[e[i,]]
  maker = paste(B[1],'|',B[2], sep = '')
  preparedata('/Users/HoltLab/Documents/R/RNA_seq_ignacio_11.28',maker, paste('condition', B[2], sep = ''))
  }
 
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')


sleuth_live(so)
# I like butts
#----
##additional analysis

#Stor Plots as object
PLOT = plot_pca(so)

# for ven diagram stuff
ven = read.csv2('VennDiagram_up_in_starvation.txt', col.names = F, sep = '\t')
a = ven$FALSE.
for(i in 1:length(a)){
  b = as.character(a[i])
  c = unlist(strsplit(b,': '))
  d = c[2]
  e = data.frame(strsplit(d, ','))
  colnames(e) = 'genes'
  write.table(e, file = paste(c[1], '.txt', sep = ''), sep = '\n', row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}
