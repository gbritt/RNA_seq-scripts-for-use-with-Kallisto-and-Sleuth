#Sleuth Analysis in R
library('biomaRt')
library('RCurl')
preparedata = function(base_dir, Ids){
  

  base_dir = '/Volumes/HoltMacPro1Backup/Greg_Gem_Data' #base_dir  #
  sample_id =dir(file.path(base_dir,"strains"))

  kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "strains", id, "kallisto"))
  kal_dirs
  
  sample_condition = t(data.frame('WTDex','WTDex','SpidrinDex','SpidrinDex','H2ADex','H2ADex','deltaQDex','deltaQDex','NullDex', 'NullDex','WTStarve','WTStarve','SpidrinStarve','SpidrinStarve','H2AStarve','H2AStarve','deltaQStarve','deltaQStarve','NullStarve','NullStarve','WTDex','WTDex','SpidrinDex','SpidrinDex','H2ADex','H2ADex','deltaQDex','deltaQDex','NullDex','NullDex','WTStarve','WTStarve','SpidrinStarve','SpidrinStarve','H2AStarve','H2AStarve','deltaQStarve','deltaQStarve','NullStarve','NullStarve'))
                     rownames(sample_condition) = 1:40

                     
  s2c = data.frame(sample_id,sample_condition)
  colnames(s2c) = c("sample","condition")

  s2c <- dplyr::mutate(s2c, path = kal_dirs)
  s2c_cut = s2c[grep(Ids, s2c[,2]),]
  
  s2c_cut$condition = factor(s2c_cut$condition)
  s2c_cut$sample = factor(s2c_cut$sample)
  
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "scerevisiae_gene_ensembl",
                           host = 'ensembl.org')
  so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
  
  
  so <- sleuth_prep(s2c_cut, ~ condition)
  so <- sleuth_fit(so)
  #least ratio test
  so <- sleuth_fit(so, ~1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
 c\
}


  

sleuth_live(so)

----
#analysis via Dave Truong