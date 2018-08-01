#Ignacio RNAseq 1.2.18 - Heatmap

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

so <- sleuth_prep(metadata, target_mapping = t2g,
                  aggregation_column = 'ens_gene', 
                  extra_bootstrap_summary = TRUE)#Make a regression model for the factors in the metadata file
so <- sleuth_fit(so, ~Condition,'full') # model gene expression as responding to this factor - (fits the full model)
so = sleuth_fit(so, ~BiologicalReplicate + SequencingRun, "reduced")#create another model where expression is not dependent on any factor
# Run a likelihood ratio test between the two models to see what transcripts appear to really be affected by these units
so <- sleuth_lrt(so, 'reduced', 'full')

table = kallisto_table(so, use_filtered = FALSE, normalized = TRUE,
                       +                include_covariates = TRUE)

SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = FALSE)

top_genes = sleuth_significant[1]
heatmap_genes = top_genes[1:50,1]

so$est_counts_sf

plot_transcript_heatmap(so,heatmap_genes, units = 'tpm')

test_table = sleuth_to_matrix(so, 'obs_norm','est_counts')

plot_bootstrap(so, "YHR055C",
               units = "scaled_reads_per_base", color_by = "Condition")

tmp = installed.packages()

installedpackages = as.vector(tmp[is.na(tmp[,"Priority"]), 1])
save(installedpackages, file="~/Desktop/installed_packages.rda")

