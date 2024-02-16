library("TwoSampleMR")

##### read exposure file ################
exp0 <- read_exposure_data(
  filename = "CAMK1.txt",
  sep = "\t",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  #eaf_col = "Freq1",
  effect_allele_col = "a1",
  other_allele_col = "a2",
  pval_col = "p",
  min_pval = 0,
  id_col = "trait"
)

##### read outcome file #################

##### FinnGen file using R7 ###################
##### UKB file  download from UK biobank #############

PHD_out <- read_outcome_data(snps = exp0$SNP,
                                  filename = "***.txt",
                                  sep = "\t",
                                  snp_col = "rsids",
                                  beta_col = "beta",
                                  se_col = "sebeta",
                                  effect_allele_col = "alt",
                                  other_allele_col = "ref",
                                  eaf_col = "af_alt",
                                  pval_col = "pval",
                                  chr_col = "#chrom")

PHD_out$outcome <- "PHD"

#### perform MR analysis ################
harmonise <- harmonise_data(exposure_dat = exp0,outcome_dat = FINN_PHD_out)
clump <- clump_data(harmonise, clump_r2 = 0.01)
mr_results <- mr(clump,method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression","mr_wald_ratio"))
het <- mr_heterogeneity(clump)
ple <- mr_pleiotropy_test(clump)
