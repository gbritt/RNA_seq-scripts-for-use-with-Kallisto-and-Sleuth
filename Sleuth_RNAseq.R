#Sleuth Analysis in R
library('biomaRt')
library('RCurl')
library('gtools')

setwd('/Users/HoltLab/Documents/R/Experiments/RNA_seq_ADH2_Nacho_12.28.17/')
preparedata = function(base_dir, Ids, compare){
  library('biomaRt')
  library('RCurl')
  library('sleuth')
  setwd('/Users/HoltLab/Documents/R/Experiments/RNA_seq_ADH2_Nacho_12.28.17/')
  options(mc.cores = 6L) # use the 6 cores of the mac pro
  #prepare for loading
  base_dir =    '/Users/HoltLab/Documents/R/RNA_seq_ignacio_11.28' #base_dir
  #base_dir =    '/Users/HoltLab/Documents/R/RNA_Seq_with_tech_reps' #self explanatory
  #
  #sample_id =dir(base_dir)
  
  #kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, "kallisto"))
  #kal_dirs
  
  # sample_condition = t(data.frame('1WTDex','2SpidrinDex','3H2ADex','4deltaQDex', '5NullDex','6WTStarve','7SpidrinStarve','8H2AStarve','9deltaQStarve','9NullStarve','1WTDex','2SpidrinDex','3H2ADex','4deltaQDex','5NullDex','6WTStarve','7SpidrinStarve','8H2AStarve','9deltaQStarve','9NullStarve'))
  #   (data.frame(c(rep('Dex',5), rep('Starve',5),rep('Dex',5), rep('Starve',5) ))) #for comparing WT DEX
  # rownames(sample_condition) = 1:20
  # sample_conditonn = c('1WTDex','2SpidrinDex','3H2ADex','4deltaQDex', '5NullDex','6WTStarve','7SpidrinStarve','8H2AStarve','9deltaQStarve','9NullStarve','1WTDex','2SpidrinDex','3H2ADex','4deltaQDex','5NullDex','6WTStarve','7SpidrinStarve','8H2AStarve','9deltaQStarve','9NullStarve')
  # print(Ids)
  # s2c = data.frame(sample_id,sample_condition)
  # names4scatter = c('WTDex','SpidrinDex','H2ADex','deltaQDex', 'NullDex','WTStarve','SpidrinStarve','H2AStarve','deltaQStarve','NullStarve','WTDex_rep','SpidrinDex_rep','H2ADex_rep','deltaQDex_rep','NullDex_rep','WTStarve_rep','SpidrinStarve_rep','H2AStarve_rep','deltaQStarve_rep','NullStarve_rep')
  # names4scatter = c('WTDex','WTDex2','SpidrinDex','SpidrinDex2','H2ADex','H2ADex2','deltaQDex','deltaQDex2','NullDex','NullDex2','WTStarve','WTStarve2','SpidrinStarve','SpidrinStarve2','H2AStarve','H2AStarve2','deltaQStarve','deltaQStarve2','NullStarve','NullStarve2','WTDex_rep','WTDex_rep2','SpidrinDex_rep','SpidrinDex_rep2','H2ADex_rep','H2ADex_rep2','deltaQDex_rep','deltaQDex_rep2','NullDex_rep','NullDex_rep2','WTStarve_rep','WTStarve_rep2','SpidrinStarve_rep','SpidrinStarve_rep2','H2AStarve_rep','H2AStarve_rep2','deltaQStarve_rep','deltaQStarve_rep2','NullStarve_rep','NullStarve_rep2')
  # #if you have tech replicates
  # 
  #  s2c = data.frame(names4scatter,sample_condition)
  # 
  # 
  # colnames(s2c) = c("sample","condition")
  # 
  # s2c <- dplyr::mutate(s2c, path = kal_dirs)
  # growthcond = c(rep('Dex',5), rep('Starve',5),rep('Dex',5), rep('Starve',5) )
  # s2c = cbind(s2c, growthcond)
  # s2c_cut = s2c[grep(Ids, s2c[,2]),]
  # 
  # s2c_cut$condition = factor(s2c_cut$condition)
  # s2c_cut$sample = factor(s2c_cut$sample)
# for making metadata, will import instead.
  #metadata = read.table("ignacioMetaDatafull.txt", sep = "\t", header = TRUE)
  metadata = read.table("IgnaciometadataCombined.txt", sep = "\t", header = TRUE)
  metadata$path = as.character(metadata$path)
  metadata$Condition = as.factor(metadata$Condition)
  #factor(metadata$Strain)
  #relevel(metadata$sample, "")
  # add annotations
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "scerevisiae_gene_ensembl",
                           host = 'ensembl.org')
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name", "description"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  
  #so <- sleuth_prep(s2c_cut, ~condition, target_mapping = t2g, aggregation_column = 'ens_gene')
  # to run it in gene based analysis mode
  
  #set up for individual comparisons
  conA <- which(metadata$BiologicalReplicate == "WTDex")
  conB <- which(metadata$BiologicalReplicate == "WTStarve")
  
  conA <- which(metadata$BiologicalReplicate == "H2ADex")
  conB <- which(metadata$BiologicalReplicate == "H2AStarve")
  
  conA <- which(metadata$BiologicalReplicate == "deltaQDex")
  conB <- which(metadata$BiologicalReplicate == "deltaQStarve")
  
  conA = which(metadata$BiologicalReplicate == "SpidrinDex")
  conB = which(metadata$BiologicalReplicate == "SpidrinStarve")
  
  conA = which(metadata$BiologicalReplicate == "NullDex")
  conB = which(metadata$BiologicalReplicate == "NullStarve")
  
  metadata_subset = metadata[c(conA,conB),]
  metadata_subset = droplevels(metadata_subset)
  #prep/fit
  so <- sleuth_prep(metadata, target_mapping = t2g,
                    aggregation_column = 'ens_gene', 
                    extra_bootstrap_summary = TRUE)#Make a regression model for the factors in the metadata file
  so <- sleuth_fit(so, ~Condition,'full') # model gene expression as responding to this factor - (fits the full model)
  so = sleuth_fit(so, ~BiologicalReplicate + SequencingRun, "reduced")#create another model where expression is not dependent on any factor
  # Run a likelihood ratio test between the two models to see what transcripts appear to really be affected by these units
  so <- sleuth_lrt(so, 'reduced', 'full')
  results_table_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
  #return(so)
 # so
  
  #perform differential expression analysis
  #need to choose the comparison so I can list what they are with:
  so$fits[["full"]]$design_matrix
  
  #then choose to make a comparison
  so <- sleuth_wt(so, 'BiologicalReplicateNullStarve')
  results_table_wt <- sleuth_results(so, 'BiologicalReplicateNullStarve')
  
  
  #wald test doesnt give log fold change so here we can get that with
  #LOG(POWER(2.71828, Beta), 2)    #doesn't work?
  
  #--------------
  #  We want the Wald test results for those genes that also passed the more stringent (when there
  # 's only two factor levels) LRT.  Let's get the names of the transcripts that are in both lists 
  # (with multiple testing corrected p-value, a.k.a. q-value, less than 0.05).
  #  return(so) # return sleuth object when not directly writing files
  d5.lrt.sig_ids <- results_table_lrt$target_id[which(results_table_lrt$qval < 0.1)]
  d5.wt.sig_ids <- results_table_wt$target_id[which(results_table_wt$qval < 0.05)]
  shared_ids <- d5.wt.sig_ids[d5.wt.sig_ids %in% d5.lrt.sig_ids] 
  # now just grabbing the rows from the wald test with these transcripts
  shared_results <- results_table_wt[results_table_wt$target_id %in% shared_ids,]
  write.csv(shared_results, file="deltaQStarvevsdeltaQDexTranscripts.csv")
  #----------------
  # need to figure out how to test for differences between conditions - when do I correct for multiple testing and when do I not?
  
  ### write results 
  compare = ("BiologicalReplicateNullStarve")
  results_table <- sleuth_results(so, compare)
  
  results_ordered <- results_table[order(results_table$qval),]
  table(results_ordered$qval <= 0.01) #set to what you want, normally I use 0.001
  Ids = paste("NullStarve","NullDex", sep = "_versus_")
  # Writes tables for dif expressed genes
  write.table( subset(results_ordered, qval <= 0.01), file= paste('DE',Ids,'.qval_0.01.txt', sep = ''), sep="\t",row.names=F, quote=F)
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
  source("GB_Volcano_sleuth.R")
  volcano = GB_volcano(so, compare, test_type = "wt", which_model = "full",
                         sig_level = 0.0001, point_alpha = 0.4, sig_color = "red",
                         highlight = NULL) #not working right now
  volcano = volcano + scale_x_continuous(breaks=seq(-10,10), limits = c(-8,8)) + scale_y_continuous( limits = c(0,90)) 
  
  pdf(paste('Volcano',Ids,'.pdf',sep = ""))
  print(volcano)
  dev.off()
 
}

so = preparedata('/Users/HoltLab/Documents/R/RNA_seq_ignacio_11.28','1WTDex|6WTStarve', 'condition6WTStarve')
so = preparedata('/Users/HoltLab/Documents/R/RNA_Seq_with_tech_reps',A, 'condition6WTStarve')

#for running comparisions withing Dex and Starve groups
A = c('1WTDex','2SpidrinDex','3H2ADex','4deltaQDex', '5NullDex','6WTStarve','7SpidrinStarve','8H2AStarve','9deltaQStarve','9NullStarve','1WTDex','2SpidrinDex','3H2ADex','4deltaQDex','5NullDex','6WTStarve','7SpidrinStarve','8H2AStarve','9deltaQStarve','9NullStarve')
#least ratio test
a = combinations(n=5,r=2) # for selecting first half (1 - 5)
c = c(5,5)
D = rbind(c,c,c,c,c,c,c,c,c,c) # array of 5's to add to the 1-6 to select second half
rownames(D) = 1:10
e = D + a #added
#use a for first half and e for second half 

# for doing normal comparisons
N = matrix(c(1:5, 6:10),5) #dex vs starve


for(i in 1:length(e[,1])){
  w = e+10
  B = A[w[i,]] #use a for first half, e for second half crowsswise comparisons, and N for starve vs dex within conditions
  maker = paste(B[1],'|',B[2], sep = '')
  #print(maker) # to test
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
  #plot scatter plots
  
  
  
}

#least ratio test
a = combinations(n=10,r=2) # for selecting first half (1 - 5)
c = c(10,10)
D = t(rbind(rep(c,20))) # array of 5's to add to the 1-6 to select second half


for(i in 1:44){
  c = rbind(c,c(10,10))
}
rownames(c) = 1:45
twenty = c *2
thirty = c*3

tens = c + a #added
twenties = twenty + a
thirties = thirty + a
names4scatter = c('WTDex','WTDex2','SpidrinDex','SpidrinDex2','H2ADex','H2ADex2','deltaQDex','deltaQDex2','NullDex','NullDex2','WTStarve','WTStarve2','SpidrinStarve','SpidrinStarve2','H2AStarve','H2AStarve2','deltaQStarve','deltaQStarve2','NullStarve','NullStarve2','WTDex_rep','WTDex_rep2','SpidrinDex_rep','SpidrinDex_rep2','H2ADex_rep','H2ADex_rep2','deltaQDex_rep','deltaQDex_rep2','NullDex_rep','NullDex_rep2','WTStarve_rep','WTStarve_rep2','SpidrinStarve_rep','SpidrinStarve_rep2','H2AStarve_rep','H2AStarve_rep2','deltaQStarve_rep','deltaQStarve_rep2','NullStarve_rep','NullStarve_rep2')
# ^ with technical replicates
names4scatter = c('WTDex','SpidrinDex','H2ADex','deltaQDex', 'NullDex','WTStarve','SpidrinStarve','H2AStarve','deltaQStarve','NullStarve','WTDex_rep','SpidrinDex_rep','H2ADex_rep','deltaQDex_rep','NullDex_rep','WTStarve_rep','SpidrinStarve_rep','H2AStarve_rep','deltaQStarve_rep','NullStarve_rep')
for(i in 1:45){
  B = names4scatter[thirties[i,]] #use a for first half, e for second half crowsswise comparisons, and N for starve vs dex within conditions
  #scatt = plot_scatter(so,names4scatter[i],names4scatter[i+10])# for normal
  scatt = plot_scatter(so,B[1],B[2]) #for crosswise # need to add cartesian(ylim(y)) c
  pdf(paste('Scatter',paste(B[1],'_',B[2]),'.pdf',sep = ""))
  print(scatt)
  dev.off()
}

for(i in 1:length(e[,1])){
  B = A[e[i,]] #use a for first half, e for second half crowsswise comparisons, and N for starve vs dex within conditions
  maker = paste(B[1],'|',B[2], sep = '')
  #print(maker) # to test
  preparedata('/Users/HoltLab/Documents/R/RNA_seq_ignacio_11.28',maker, paste('condition', B[2], sep = ''))
}
#testing for rank of columns - dont want columns that are redundant
#The rank of a matrix columns is the number of columns that are independent of all the others. If the rank is smaller than the number of columns, then the LSE are not unique. In R. we can obtain the rank of matrix with the function qr, which we will describe in more detail in a following section.
Sex <- c(0,0,0,0,1,1,1,1)
A <-   c(1,1,0,0,0,0,0,0)
B <-   c(0,0,1,1,0,0,0,0)
C <-   c(0,0,0,0,1,1,0,0)
D <-   c(0,0,0,0,0,0,1,1)
X <- model.matrix(~Sex+A+B+C+D-1)
cat("ncol=",ncol(X),"rank=", qr(X)$rank,"\n")

