# Freebayes pipeline ####
# If you are a coder, this following probably looks like a mess. But it is functional and works
# GATK processed until snpFiltR in seperate R. file chp_3_analysis.GATK.R
# Load packages####

analysis_pkgs <- c("farff","tidyverse", "openxlsx", "furrr", "progressr", "here",
                   "tictoc", "scales", "glue",'plotly', 'readxl',# general packages
                   "vcfR","adegenet", "seqinr", "poppr", "mmod", 
                   "treemap", "ape", 'SNPfiltR', # pop gen packages
                   'geneHapR','pegas', 'dplyr', 'ggrepel')

#pak::pak(analysis_pkgs)

pacman::p_load(char = basename(analysis_pkgs), update = FALSE, install = FALSE)

# Load Ido's magic functions####
# These functions are for filtering samples from a vcf without erasing the format section
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", 
                      filename = "util.R")


# Read in VCF file for vcf comparisons ####
# No filtering
vcf1 <- read.vcfR(here("data/freebayes/A_rabiei_2024_Murdoch_WGRS_ArME14_v2.bwa2.fb.diploid.vcf.gz"))
# Filtered: no heterozygotes; quality > 30; bi-allelic polymorphic
vcf3 <- read.vcfR(here("data/freebayes/A_rabiei_2024_Murdoch_WGRS_ArME14_v2.bwa2.fb.diploid.Q30.noHet.poly.recode.vcf.gz"))
# Filtered: no heterozygotes; quality > 30; bi-allelic polymorphic; SNPs only
vcf5 <- read.vcfR(here("data/freebayes/A_rabiei_2024_Murdoch_WGRS_ArME14_v2.bwa2.fb.diploid.snps.Q30.noHet.poly.recode.vcf.gz"))

# Compare samples in separate vcfs
# Skip this step unless it is the first time processing data
# Compare vcf1 and vcf2, FB vs. GATK to determine where the extra FB sample is coming from
# 66 (FB) vs. 65 (GATK)
vcf1_samples <- colnames(vcf1@gt)[-1]
samples_names_tbl_1 <- tibble(sample=vcf1_samples)

vcf2_samples <- colnames(vcf2@gt)[-1]
samples_names_tbl_2 <- tibble(sample=vcf2_samples)

# Identify samples in vcf2 but not in vcf1
missing_in_vcf1 <- samples_names_tbl_2 %>%
  anti_join(samples_names_tbl_1, by = "sample")
# Identify samples in vcf1 but not in vcf2
missing_in_vcf2 <- samples_names_tbl_1 %>%
  anti_join(samples_names_tbl_2, by = "sample")
# Samples Ar23766 is missing
# It is likely that this genotype did not make it through GATK filtering
# Extract it from the FB VCFs and override them
# This can also be used to extract and override the VCFs once the extra data file from the replicates has been determined

# Filtering samples out of vcfs ####
# Apply to each vcf
# filtered_vcf <- read.vcfR('input/freebayes/A_rabiei_2024_Murdoch_WGRS_ArME14_v2.bwa2.fb.diploid (1).vcf.gz')
 

vcf_samples <- colnames(filtered_vcf@gt)[-1]
samples_to_exclude <- c("23755")
exclude_samples <- vcf_samples %in% samples_to_exclude
selected_columns <- c('FORMAT',vcf_samples[!exclude_samples]) 
filtered_vcf <- filtered_vcf[, selected_columns]

# filter vcf by sample function
filter_vcfr_samples <- function(vcf, include_samples, exclude_samples) {
  require(vcfR)
  rlang::check_exclusive(include_samples, exclude_samples)
  sample_names <- colnames(vcf@gt)[-1]
  if (missing(include_samples)) {
    include_samples <- sample_names[!sample_names %in% exclude_samples]
  }
  return(vcf[,c('FORMAT', include_samples)])
}

filtered_vcf <- filter_vcfr_samples(filtered_vcf, exclude_samples = exclude_samples)

# Write out the new VCF file
vcfR::write.vcf(vcf1, file = 'data/freebayes/vcf1_FB.vcf.gz' )
vcfR::write.vcf(vcf3, file = 'data/freebayes/vcf3_FB.vcf.gz')
vcfR::write.vcf(vcf5, file = 'data/freebayes/vcf5_FB.vcf.gz')

# Determine the error rate and heterogeneity of unfiltered vcf files ####
# This takes the VCFs directly out of the file, no-need to write it in
# Calculating the 
# filter vcfR samples by pattern
grep_vcfr_samples <- function(vcf, pattern, invert = FALSE) {
  require(vcfR)
  sample_names <- colnames(vcf@gt)[-1]
  return(vcf[,c('FORMAT', grep(pattern, x = sample_names, 
                               value = TRUE, invert = invert))])
}

#Approach: Create a function with the vcf a an input in addition to a singular sample name to define it, of which the objects entered must be standardized. Function can then be included to rotate through the entered [listed] vcfs to calculate the error rate (defines by samples in sample_names) and heterozygotes.as a proportion of each vcf 
sample_name <- "Ar23767"
sample_names <- c("Ar23767", "Ar23401", "Ar23653", "AR0242", "AR0052")

# Define the function
error_rate.prop <- function(gt,vcf, sample_name) {
  # vcf = read.vcfR("output/vcf6_filtered_inform.vcf")
  # gt = extract.gt(vcf)
  vcf_samples <- colnames(gt)[-1]
  selected_columns <- c('FORMAT', grep(sample_name, vcf_samples, value = TRUE))
  vcf_reps <- vcf[, selected_columns]
  
  gt_reps <- extract.gt(vcf_reps) %>% as_tibble() %>% rowwise() %>% 
    mutate(identical = c_across(everything()) %>% n_distinct(na.rm = TRUE) == 1) %>% 
    ungroup()
  error_rate <- table(gt_reps$identical) %>% prop.table()
  results <- tibble(error_rate = 1 - error_rate["TRUE"])
  
  return(results)
}
#This then enters in the vcf files and applies the above 
vcf_files <- list.files(path = 'data/freebayes/', pattern = 'FB.vcf.gz', full.names = TRUE) %>% 
  grep('indels',., invert = TRUE, value = TRUE)
#Designates the vcf files as a list
het_error_sum <- vector('list', length = length(vcf_files))
#Loops through the vcf list using the above defined function
for (f in vcf_files) {
  vcf <- read.vcfR(f)
  gt <- extract.gt(vcf)
  het_prop <- gt %>% is_het() %>% table() %>% prop.table()
  vcf_results <- sample_names %>% map(.f = ~error_rate.prop(gt, vcf, .x)) %>% 
    setNames(sample_names) %>% 
    list_rbind(., names_to = 'sample') %>% 
    mutate(het_prop = het_prop[2], file = f)
  het_error_sum[[basename(f)]] <- vcf_results
}

#To see the results show het_error_sum
#To get the error rate average calculate as below
#vcf1
mean(het_error_sum$vcf1_FB.vcf$error_rate)
#vcf3
mean(het_error_sum$vcf3_FB.vcf$error_rate)
#vcf5
mean(het_error_sum$vcf5_FB.vcf$error_rate)

# Calculating missingness for each vcf ####
# Function to calculate average missingness

calculate_average_missingness <- function(vcf) {
  gt_matrix <- extract.gt(vcf) %>%
    as.matrix() %>%
    {.[. == "."] <- NA; .}
  
  missing_per_genotype <- colSums(is.na(gt_matrix))
  total_genotypes <- nrow(gt_matrix)
  missing_rate_per_genotype <- missing_per_genotype / total_genotypes
  mean(missing_rate_per_genotype)
}
# List of VCF files
vcf_list <- list(vcf1, vcf3, vcf5)
# Initialize an empty list to store results
average_missingness_list <- list()
# Loop through each VCF in the list and calculate average missingness
for (vcf in vcf_list) {
  average_missingness <- calculate_average_missingness(vcf)
  average_missingness_list <- c(average_missingness_list, list(average_missingness))
}


# snpFiltR and vcfR post variant filtering
#Determining missingness for each vcf file####
#vcf1
dp <- extract.gt(vcf1, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

popmap <- data.frame(
  id = colnames(vcf1@gt)[2:length(colnames(vcf1@gt))],
  pop = substr(colnames(vcf1@gt)[2:length(colnames(vcf1@gt))], 3, 9)
)
#vcf3
dp <- extract.gt(vcf3, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

popmap3 <- data.frame(
  id = colnames(vcf3@gt)[2:length(colnames(vcf3@gt))],
  pop = substr(colnames(vcf3@gt)[2:length(colnames(vcf3@gt))], 3, 9)
)
#vcf5
dp <- extract.gt(vcf5, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

popmap5 <- data.frame(
  id = colnames(vcf5@gt)[2:length(colnames(vcf5@gt))],
  pop = substr(colnames(vcf5@gt)[2:length(colnames(vcf5@gt))], 3, 9)
)

#vcf1
#Visualize distributions
hard_filter(vcfR = vcf1)
vcf1 <-hard_filter(vcfR=vcf1, depth = 5, gq = 30)
#vcf1: 
#3.75% of genotypes fall below a read depth of 5 and were converted to NA
#2.22% of genotypes fall below a genotype quality of 30 and were converted to NA
#visualize and pick appropriate max depth cutoff
max_depth(vcf1)
vcf1 <- max_depth(vcf1, maxdepth = 120)
#1.37% of SNPs were above a mean depth of 120 and were removed from the vcf
#Trying to reduce the percentage of missing data to 0
missing_by_sample(vcfR = vcf1, popmap = popmap )
vcf1<-missing_by_sample(vcfR=vcf1, cutoff = .5)
#filter biallelic SNPs
vcf1 <- filter_biallelic(vcf1)
#remove invariant sites generated by dropping individuals. Work with one for the training data but then iwth two maybe for the 
vcf1 <- min_mac(vcf1, min.mac = 2)
#93.65% of SNPs fell below a minor allele count of 2 and were removed from the VCF
missing_by_snp(vcf1, cutoff = NULL)
vcf_sample_ids <- colnames(vcf1@gt)[2:length(colnames(vcf1@gt))]
filtered_popmap <- popmap[popmap$id %in% vcf_sample_ids, ]
miss.snp <-assess_missing_data_pca(vcfR=vcf1, popmap = filtered_popmap, thresholds = .99, clustering = FALSE)
filtered_vcf <- missing_by_snp(vcf1, cutoff = 0.99)
vcfR::write.vcf(filtered_vcf, "output/vcf1_filtered.vcf.gz")

#vcf3
hard_filter(vcfR = vcf3)
vcf3 <-hard_filter(vcfR=vcf3, depth = 5, gq = 30)
#vcf3: 
#7.59% of genotypes fall below a read depth of 5 and were converted to NA
#8.41% of genotypes fall below a genotype quality of 30 and were converted to NA
#visualize and pick appropriate max depth cutoff
max_depth(vcf3)
vcf3 <- max_depth(vcf3, maxdepth = 120)
#1.76% of SNPs were above a mean depth of 120 and were removed from the vcf
#Trying to reduce the percentage of missing data to 0
missing_by_sample(vcfR = vcf3, popmap = popmap3 )
vcf3<-missing_by_sample(vcfR=vcf3, cutoff = .5)
#5 samples from a cuttoff of 0.5
#filter biallelic SNPs
vcf3 <- filter_biallelic(vcf3)
#remove invariant sites generated by dropping individuals. Work with one for the training data but then iwth two maybe for the 
vcf3 <- min_mac(vcf3, min.mac = 2)
#93.65% of SNPs fell below a minor allele count of 2 and were removed from the VCF
missing_by_snp(vcf3, cutoff = NULL)
vcf_sample_ids <- colnames(vcf3@gt)[2:length(colnames(vcf3@gt))]
filtered_popmap3 <- popmap3[popmap3$id %in% vcf_sample_ids, ]
miss.snp <-assess_missing_data_pca(vcfR=vcf3, popmap = filtered_popmap3, thresholds = .99, clustering = FALSE)
filtered_vcf <- missing_by_snp(vcf3, cutoff = 0.99)
#12.58% of SNPs fell below a completeness cutoff of 0.95 and were removed from the VCF
vcfR::write.vcf(filtered_vcf, "output/vcf3_filtered.vcf.gz")

#vcf5
hard_filter(vcfR = vcf5)
vcf5 <-hard_filter(vcfR=vcf5, depth = 5, gq = 30)
#vcf5: 
#5.02% of genotypes fall below a read depth of 5 and were converted to NA
#7.53% of genotypes fall below a genotype quality of 30 and were converted to NA
#visualize and pick appropriate max depth cutoff
max_depth(vcf5)
vcf5 <- max_depth(vcf5, maxdepth = 120)
#2.47% of SNPs were above a mean depth of 120 and were removed from the vcf
#Trying to reduce the percentage of missing data to 0
miss <- missing_by_sample(vcfR = vcf5, popmap = popmap5 )
vcf5<-missing_by_sample(vcfR=vcf5, cutoff = .5)
miss$unfiltered.stats
write.xlsx(miss$unfiltered.stats, file = "output/vcf5_miss.stats.xlsx", rowNames = FALSE)
overall_mean_depth <-  mean(miss$unfiltered.stats$mean.depth, na.rm = TRUE)
#5 samples from a cuttoff of 0.5
#filter biallelic SNPs
vcf5 <- filter_biallelic(vcf5)
#remove invariant sites generated by dropping individuals. Work with one for the training data but then iwth two maybe for the 
vcf5 <- min_mac(vcf5, min.mac = 2)
#39.81% of SNPs fell below a minor allele count of 2 and were removed from the VCF
missing_by_snp(vcf5, cutoff = NULL)
vcf_sample_ids <- colnames(vcf5@gt)[2:length(colnames(vcf5@gt))]
filtered_popmap5 <- popmap5[popmap5$id %in% vcf_sample_ids, ]
miss.snp <-assess_missing_data_pca(vcfR=vcf5, popmap = filtered_popmap5, thresholds = .99, clustering = FALSE)
filtered_vcf <- missing_by_snp(vcf5, cutoff = 0.99)
#12.58% of SNPs fell below a completeness cutoff of 0.95 and were removed from the VCF
vcfR::write.vcf(filtered_vcf, "output/vcf5_filtered.vcf.gz")

#View new depth
vcf1_filtered <- read.vcfR("output/vcf1_filtered.vcf.gz")
dp <- extract.gt(vcf1_filtered, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

vcf3_filtered <- read.vcfR("output/vcf3_filtered.vcf.gz")
dp <- extract.gt(vcf1_filtered, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

vcf5_filtered <- read.vcfR("output/vcf5_filtered.vcf.gz")
dp <- extract.gt(vcf5_filtered, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)
extract.gt(vcf5_filtered, "DP")
#Output of dpeth stats
dp_matrix <- extract.gt(vcf5_filtered, element = "DP", as.numeric = TRUE)
average_depth_per_sample <- colMeans(dp_matrix, na.rm = TRUE) #%>% mean()


miss_df <- miss$unfiltered.stats %>%
  select(-samples, -samples.2, -AverageDepth) %>%  # Remove unwanted columns
  rename(samples = samples.1)
         


depth_df <- data.frame(samples = names(average_depth_per_sample),
                       AverageDepth = as.numeric(average_depth_per_sample))

# Merge with miss$unfiltered.stats by sample name
merged_df <- left_join(miss_df, depth_df, by = "samples")
                               
write.xlsx(merged_df, file = "output/vcf5_miss.stats.xlsx", rowNames = FALSE)
# Error rate after filtering#
# From this point, include the GATK samples
#This then enters in the vcf files and applies the above 
sample_names <- c("Ar23767", "Ar23401", "AR0242", "AR0052")
vcf_files <- list.files(path = 'output/', pattern = 'filtered_inform.vcf.gz$', full.names = TRUE) %>% 
  grep('indels',., invert = TRUE, value = TRUE)
#Designates the vcf files as a list
het_error_sum <- vector('list', length = length(vcf_files))
#Loops through the vcf list using the above defined function
for (f in vcf_files) {
  vcf <- read.vcfR(f)
  gt <- extract.gt(vcf)
  het_prop <- gt %>% is_het() %>% table() %>% prop.table()
  vcf_results <- sample_names %>% map(.f = ~error_rate.prop(gt, vcf, .x)) %>% setNames(sample_names) %>% list_rbind(., names_to = 'sample') %>% mutate(het_prop = het_prop[2], file = f)
  het_error_sum[[f]] <- vcf_results
}
#vcf1
mean(het_error_sum$`output/vcf1_filtered_inform`$error_rate)
#vcf2
mean(het_error_sum$`output/vcf2_filtered_inform`$error_rate)
#vcf3
mean(het_error_sum$`output/vcf3_filtered_inform`$error_rate)
#vcf4
mean(het_error_sum$`output/vcf4_filtered_inform`$error_rate)
#vcf5
mean(het_error_sum$`output/vcf5_filtered`$error_rate)
#vcf6
mean(het_error_sum$`output/vcf6_filtered`$error_rate)

# Determine informative loci####

informative_genind_vcf <- vcfR2genind(vcf5_filtered, return.alleles = TRUE)
ploidy(informative_genind_vcf) <- 1
gen.vcf.inform <- informloci(informative_genind_vcf, quiet = FALSE)
inform_loci <- locNames(gen.vcf.inform)
vcf5_filtered@fix 
keep_loci <- vcf5_filtered@fix %>% as_tibble() %>% mutate(lociname = paste(CHROM,POS,sep = "_")) %>% 
  pull(lociname) %in% inform_loci
filtered_vcf <- vcf5_filtered[keep_loci,]

vcfR::write.vcf(vcf5_filtered_inform, "output/vcf5_filtered_inform.vcf.gz")

# Poppr analysis####
#Converting vcfs to a genind onject for easier handling with poppr####

vcf1_filtered_inform <- read.vcfR("output/vcf1_filtered_inform.vcf.gz")
vcf3_filtered_inform <- read.vcfR("output/vcf3_filtered_inform.vcf.gz")
vcf5_filtered_inform <- read.vcfR("output/vcf5_filtered_inform.vcf.gz")

gen_obj1 <- vcfR2genind(vcf1_filtered_inform) 
gen_obj3 <- vcfR2genind(vcf3_filtered_inform)
gen_obj5 <- vcfR2genind(vcf5_filtered_inform)

ploidy(gen_obj1) <- 1
ploidy(gen_obj3) <- 1
ploidy(gen_obj5) <- 1
#Check that the genind objects have their correct ploidy designated Just check with one
tab(gen_obj1)[1:5,1:5]
#Convert genind object into a dataframe
X <- genind2df(gen_obj1)
head(X)
#Coverts is back to a genind object but ensures the ploidy is 1 and that it knows that only one character codes the allele so it does not assume
gen_obj1 <- df2genind(X, ploidy = 1, ncode=1)
#Retrieves the name of the individuals from the file to have a look at and make sure they are correct
indNames(gen_obj1)
#Prints genind object to the console
gen_obj

#Check vcf5 mean

#Converts the genind object to a genclone object. Genclone is used to handle clonal data
gen_obj1 <- poppr::as.genclone(gen_obj1)
gen_obj3 <- poppr::as.genclone(gen_obj3)
gen_obj5 <- poppr::as.genclone(gen_obj5)
#re-check that the names are correct
indNames(gen_obj)
#Provides the number of alleles per locus in the total dataset, so don't freak out if you have 2 for a haploid because if there are variants that is what will be present, especially if it has already been run through the informative loci
summary(gen_obj)
#View file to see what metadata to add to existing file
#names <- colnames(vcf1_filtered_inform@gt)[-1]

#Attach the metadata####
#do this one at a time with each of the objects 1,3,5
metadata <- read_excel("input/exp_pop_training_metadata.xlsx")
colnames(metadata)
filtered_metadata <- metadata[metadata$sample_ID %in% indNames(gen_obj1), ]
genclone_df <- data.frame(genotype = indNames(gen_obj3))

joined_data <- genclone_df %>%
  left_join(filtered_metadata, by = c("genotype" = "sample_ID"))

strata(gen_obj3) <- joined_data

#Check strata
gen_obj1
gen_obj3
gen_obj5

#View the separate strata but there should be no populations at this point
gen_obj

#Set the populations where necessary, I am going to go with pathogenicity for this
popNames(gen_obj)
setPop(gen_obj) <- ~genotype
setPop(gen_obj) <- ~range
setPop(gen_obj) <- ~row
setPop(gen_obj) <- ~rep
setPop(gen_obj) <- ~host
#Now when using popNames it will tell you path group 1-5 for example
gen_obj

#Get MLG for each VCF####
#Use genind for this
#Calculate mlg grid
#If replicates do not already match calculate by pairwise dist then bitwise dist
#Filter until technical replicates match

#vcf1
gen.ml_grid <- mlg.id(gen_obj1) 
gen.ml_grid 
gen <- poppr::as.genclone(gen_obj1)
#X <- genclone2genind(gen_obj1)
X <- genind2loci(gen)
#Pairwise genetic distances
pairwise_genetic_distance <- dist.gene(X, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
hist(pairwise_genetic_distance)
# Convert the matrix to a data frame
pairdist_df <- as.data.frame(as.matrix(pairwise_genetic_distance))
# Write the data frame to a CSV file
write.csv(pairdist_df, file = "output/vcf1_pairdist_matrix_bitwise.csv", row.names = TRUE)
#Bitwise distance 
#This is simply another method of comparing the differences but instead of using a model like pairwise distance it simply counts the 1 and 0 between samples and it commonly used for SNP data and is more efficienct
adist <- bitwise.dist(gen, mat=TRUE, euclidean = FALSE, percent = FALSE)
library(MASS)
# Convert the matrix to a data frame
adist_df <- as.data.frame(as.matrix(adist))
# Write the data frame to a CSV file
write.csv(adist_df, file = "output/vcf1_adist_matrix_bitwise.csv", row.names = TRUE)
#bitwise and pairwise give essentially the same results, just good to check both
#the greatest difference in snps between technical replicated was AR0052 with 7 SNPs the difference. Aim to have this sample in the same MLG
#To filter vcf1
mlg.filter(gen.2, threshold = 0.05, distance = "bitwise.dist", threads = 1L)

#Generates threshold statistics in which clusters are formed according to nei's
(thresholds <- mlg.filter(gen.2, distance = "bitwise.dist", stats = "THRESHOLDS",
                          threshold = 1))
(pcut <- cutoff_predictor(thresholds))
#Then filters according to the set threshold (note that idf you set it to <- .2 it will fitler to the 4 haplotypes)
mlg.filter(gen.2, distance = "bitwise.dist") <- 0.008 #pcut
mlg.table(gen.2)

gen.ml_grid <- mlg.id(gen.2) 
gen_obj1 <- gen.2

#vcf3
gen.ml_grid3 <- mlg.id(gen_obj3) 
gen.ml_grid3 
gen <- poppr::as.genclone(gen_obj3)
#X <- genclone2genind(gen_obj1)
X <- genind2loci(gen)
#Pairwise genetic distances
pairwise_genetic_distance <- dist.gene(X, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
hist(pairwise_genetic_distance)
# Convert the matrix to a data frame
pairdist_df <- as.data.frame(as.matrix(pairwise_genetic_distance))
# Write the data frame to a CSV file
write.csv(pairdist_df, file = "output/vcf3_pairdist_matrix_bitwise.csv", row.names = TRUE)
#There is no need to filter vcf3

#vcf5
gen.ml_grid5 <- mlg.id(gen_obj5) 
gen.ml_grid5 
gen <- poppr::as.genclone(gen_obj5)
#X <- genclone2genind(gen_obj1)
X <- genind2loci(gen)
#Pairwise genetic distances
pairwise_genetic_distance <- dist.gene(X, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
hist(pairwise_genetic_distance)
# Convert the matrix to a data frame
pairdist_df <- as.data.frame(as.matrix(pairwise_genetic_distance))
# Write the data frame to a CSV file
write.csv(pairdist_df, file = "output/vcf5_pairdist_matrix_bitwise.csv", row.names = TRUE)
#There is no need to filter vcf5

#PCA analysis 
#vcf1
gen1 <- genclone2genind(gen_obj1)
gen1<- gen_obj1
popNames(gen_obj1)
setPop(gen_obj1) <- ~genotype
gen_obj1

#vcf3
gen3 <- genclone2genind(gen_obj3)
gen3 <- gen_obj3
popNames(gen_obj3)
setPop(gen_obj3) <- ~genotype
gen_obj3

#vcf5
gen5 <- genclone2genind(gen_obj5)
gen5 <- gen_obj5
popNames(gen_obj5)
setPop(gen_obj5) <- ~genotype
gen_obj5

allelfreq <- scaleGen(gen, NA.method="mean", scale=TRUE)

pca <- dudi.pca(allelfreq, scale = FALSE, scannf = FALSE, nf = 6)

# Need the coordinated that each isolate exist on

pcadata <- pca$li %>% rownames_to_column('genotype') %>% 
  left_join(genind_df)
pcavar <- pca$li %>% summarise(across(starts_with('Axis'), var))

percent_var <- pcavar/sum(pcavar)

# Now we create the pca plot for principal components 1 and 2
library(scales)
library(paletteer)
library(viridis)
ggplot(pcadata, aes(x= Axis1, y= Axis2, colour = genotype)) + geom_point(size = 4) +  scale_colour_viridis_d(option = 'viridis') + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(colour = "Isolate", x= glue("C1: {percent(percent_var$Axis1)}"),
       y= glue("C2: {percent(percent_var$Axis2)}")) +
  theme_bw(14)
ggsave("output/PCA_isolates_TECAN_gendistdata.pdf", width = 7, height = 6)

# Plot points labeled instead
library(ggrepel)
ggplot(pcadata, aes(x = Axis1, y = Axis2, colour = genotype, label = genotype)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = Inf) +  # Set max.overlaps to Inf to allow all labels
  scale_colour_viridis_d(option = 'viridis') +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(colour = "Isolate", x = glue("C1: {percent(percent_var$Axis1)}"),
       y = glue("C2: {percent(percent_var$Axis2)}")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
ggsave("output/vcf1_PCA_exp_pop_trainin_infomloci_wlabs.pdf", width = 7, height = 6)


# For this data I actually have the 6, but it will be interesting to see the difference
# This will set and count the number of multilocus genotypes exist accross the population (from the strata) that you have set. In the first case this is genotypes ad should show the original six.

setPop(gen_obj) <- ~genotype
gen.ml_table <- mlg.table(gen_obj)
(gen_genotype_crosspop <- mlg.crosspop(gen_obj, df=TRUE))

# DAPC ####
# Compare each vcf DAPC, will better take into consideration the clonal nature of the data and define the level of clustering
# DAPC for MLGs
# Then DAPC for six original isolates, see if it matches the phylogeny and PCA

# vcf read in
vcf1_filtered_inform <- read.vcfR("output/vcf1_filtered_inform.vcf.gz")
vcf3_filtered_inform <- read.vcfR("output/vcf3_filtered_inform.vcf.gz")
vcf5_filtered_inform <- read.vcfR("output/vcf5_filtered_inform.vcf.gz")

# vcf5
gl <- vcfR2genlight(vcf5_filtered_inform)
ploidy(gl) <- 1
 
# Set strata
metadata <- read_excel("input/exp_pop_training_metadata.xlsx")
colnames(metadata)

filtered_metadata <- metadata[metadata$sample_ID %in% indNames(gl), ]
gl_df <- data.frame(genotype = indNames(gl))
joined_data <- gl_df %>%
  left_join(filtered_metadata, by = c("genotype" = "sample_ID"))

strata_data <- joined_data %>%
  dplyr::select(genotype, everything()) %>%
  column_to_rownames(var = "genotype")

strata(gl) <- strata_data

setPop(gl) <- ~host

#any(is.na(gl))

#DAPC exploration
# Define clusters by PCA a_score
set.seed(999)
grp <- find.clusters(tab(gl, NA.method = "mean"), max.n.clust = 20)
#7
#8
head(grp$Kstat, 20)
dapc_init <- dapc(tab(gl, NA.method = "mean"), grp$grp, n.pca = 10)
#5
opt <- optim.a.score(dapc_init)
#1
#not possible, need minimum 2
dapc_final <- dapc(tab(gl, NA.method = "mean"), grp$grp, n.pca = 3 , n.da = 2)
#Can also cluster with 11 and dapc_final <- dapc(tab(gl, NA.method = "mean"), grp$grp, n.pca = 10, n.da = 5), but it does not separate the controls well, clustering 0128 and 0020 together. 
# Scatterplot for DAPC for visualisation, not the final figure
scatter(dapc_final)
myCol<-c("darkblue","purple","darkgreen","orange","red","brown" ,"magenta" ,"gold")#,"green", "lightblue", "pink"#, "grey", "darkred", "black") 
scatter(dapc_final, legend = TRUE, cleg = 1, col = myCol, clabel = 0, cex = 2.5, scree.da = TRUE, scree.pca = TRUE, posi.pca = "topleft", posi.da = "bottomleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.01, inset.pca = 0.01, posi.leg = "bottomright")
dapc_final$var
xlim <- par("usr")[1:2]
mtext(paste0("Variance explained: ", round(dapc_final$var, 4)),
      side = 3, line = -2, at = xlim[2] - 0.3 * diff(xlim), cex = 0.9)

pca_eig <- dapc_final$eig # PCA eigenvalues
da_eig <- dapc_final$eig[1:dapc_final$n.da] # DA eigenvalues (subset of PCA eig)

# Plot PCA eigenvalues
plot(pca_eig, type = "b", pch = 19, col = "blue",
     xlab = "Principal Component", ylab = "Eigenvalue",
     main = "PCA Eigenvalues")

# Plot DA eigenvalues
plot(da_eig, type = "b", pch = 19, col = "red",
     xlab = "Discriminant Axis", ylab = "Eigenvalue",
     main = "DA Eigenvalues")

# Make a table of the isoalte groupings between groups
dapc_final
dapc_final$grp
# Assuming dapc$grp is a factor or named vector
grp_df <- data.frame(
  isolate = names(dapc_final$grp),
  group = as.character(dapc_final$grp),
  stringsAsFactors = FALSE
)
# Split isolates by group
grouped <- split(grp_df$isolate, grp_df$group)
# Find the maximum group size to align columns
max_len <- max(sapply(grouped, length))
# Pad each group with NA to make equal-length columns
grouped_padded <- lapply(grouped, function(x) {
  length(x) <- max_len
  return(x)
})
# Combine into a data frame
group_table <- as.data.frame(grouped_padded, stringsAsFactors = FALSE)
# View the result
View(group_table)

# Exporting to make prettier graphs in biorender
# Use the visulaisasion below instead, but this is also an option
# Need dapc_data for downstream things so do it anyway

# Extract coordinates and cluster assignments
coords <- dapc_final$ind.coord
# Match with group assignments from dapc_group
groups <- dapc_final$grp
# Combine into a data frame
dapc_data <- data.frame(ID = rownames(coords),
                        X = coords[,1],
                        Y = coords[,2],
                        Cluster = paste("Cluster", groups))

write.csv(dapc_data, "dapc_coordinates_by_group_8.csv", row.names = FALSE)
# Check for all clusters
table(dapc_data$Cluster)

# Assuming dapc_data is already created as:
# dapc_data <- data.frame(ID = rownames(coords),
#                         X = coords[,1],
#                         Y = coords[,2],
#                         Cluster = paste("Cluster", groups))

# Make publishible figure
# Define custom colors
myCol <- c("darkblue", "purple", "darkgreen", "orange", "red", "brown", "magenta", "gold")

inds <- colnames(vcf5_filtered_inform@gt)[2:length(colnames(vcf5_filtered_inform@gt))]

isolate_table <- tibble(
  id = inds,
  isolate = sub('^2Ar', 'AR', inds),
  control = grepl("^AR", sub('^2Ar', 'AR', inds)),
  replicate = grepl("2Ar", inds),
  comment = case_when(
    grepl("Ar0", id) ~ "Control",
    grepl("^AR", id) ~ "Control",
    TRUE ~ NA_character_
  )
)

# Define control isolates, excluding specific samples
exclude_samples <- c("AR0052", "AR0052_R", "AR0242", "AR0242_R", "Ar0023", "Ar0212", "Ar0020")

controls <- isolate_table %>%
  filter(control == TRUE, !isolate %in% exclude_samples) %>%
  pull(isolate) %>%
  unique()

# Prepare plot data
plot_data <- dapc_final$ind.coord %>%
  as.data.frame() %>%
  rownames_to_column("isolate") %>%
  mutate(is_control = isolate %in% controls)
# Identify control isolates from plot_data
# Rename isolate column to match dapc_data$ID
control_labels <- plot_data %>%
  filter(is_control) %>%
  select(ID = isolate)

# Filter dapc_data for control labels only
dapc_data_labeled <- dapc_data %>%
  filter(ID %in% controls)

# Plot
ggplot(dapc_data, aes(x = X, y = Y, color = Cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  geom_label_repel(
    data = filter(dapc_data, ID %in% controls),
    aes(label = ID),
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.40,
    point.padding = 0.40,
    segment.size = 0.4,
    show.legend = FALSE
  ) +
  scale_color_manual(values = myCol) +
  labs(x = "LD1", y = "LD2", color = "Sample Type") +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 13),
    plot.title = element_blank()
  )

#Captuing the variance explained in the dapc plot by the included DAPC PCs
dapc_final$var
summary.dapc(dapc_final)
pca_sd <- dapc_final$tab %>% as.data.frame() %>% summarise(across(everything(),sd))

pca_eig <- pca_sd^2
pca_var <- pca_eig/sum(pca_eig)
vardat <- pivot_longer(pca_var, everything(), names_to = 'pca', values_to = 'prop')
ggplot(vardat,aes(x=pca, y=prop)) + geom_col(width = 0.6) + scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = expansion(mult = c(0, .1))) + labs(y='Variance explained', x= '') + theme_bw(14) +theme(panel.grid.major.x = element_blank())

#The variance explained in find.clusteres in a bar plot
pca <- glPca(gl, nf = 100)  # nf = number of PCs to retain
eig <- pca$eig
pca_var <- eig / sum(eig)

vardat <- data.frame(
  pca = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  prop = pca_var[1:10]
)


vardat$fill <- ifelse(vardat$pca %in% c("PC1", "PC2", "PC3"), "gray40", "lightgray")

# Plot with custom fill
# Export and attach attached to the other final plot figure
ggplot(vardat, aes(x = pca, y = prop, fill = fill)) +
  geom_col(width = 0.3) +
  scale_fill_identity() + # Use the fill values as-is
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = expansion(mult = c(0, .1))) +
  labs(y = 'Variance Explained', x = '') +
  theme_bw(base_size = 14) +
  theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

# Calculating variance within clusters ####
# Step 1: Define your mapping from DAPC cluster to gl cluster
pop_map <- c("1" = "6", "2" = "4", "3" = "2", "4" = "1", "5" = "5", "6" = "7", "7" = "8", "8" = "3")

# Step 2: Apply the mapping to DAPC assignments
mapped_clusters <- pop_map[as.character(dapc_final$assign)]

# Step 3: Assign the mapped clusters back to the genlight object
pop(gl) <- mapped_clusters

# Step 4: Run your variance calculation function
variance_by_cluster <- cluster_variance(dapc_final, gl)

# Step 5: Convert to percentages
variance_percentages <- as.numeric(variance_by_cluster) / sum(as.numeric(variance_by_cluster)) * 100

# Step 6: Print results
print(variance_by_cluster)
print(variance_percentages)


# Distance tree calculation ####
# Export as a FASTA file for RAxML in geneious 
# Make sure 
genotype_table <- extract.gt(vcf5_filtered_inform, return.alleles = TRUE) %>% substr(., 1, 1) %>% as.data.frame() %>% t() %>% as.data.frame()

genotype_table[is.na(genotype_table)] <- "N"
# Create a FASTA-formatted string
genotype_table_string <- genotype_table %>% unite("fasta", everything(), sep = "") %>%
  rownames_to_column("isolate") %>% mutate(fasta = paste0(">", isolate, "\n", fasta))

# Write to FASTA file
writeLines(genotype_table_string$fasta, "output/training_data_8clust_2PC.fasta")

# Add the dapc clustering to the metadata ####
df_cluster <- data.frame(
  sample_ID = names(dapc_final$grp),
  cluster = as.vector(dapc_final$grp)
)
metadata_updated <- metadata %>%
  left_join(df_cluster, by = "sample_ID") 

write.xlsx(metadata_updated, file = "output/8_clust_training_metadata.xlsx", rowNames = TRUE)

strata_data <- strata_data %>%
  rownames_to_column(var = "sample_ID") %>%
  left_join(df_cluster, by = "sample_ID") %>%
  column_to_rownames(var = "sample_ID")
strata(gl) <- strata_data
gl

# Other DAPC validation calculations ####
# Pairwise of DAPC clusters based off genetic distance between the dapc clustering. 
# If possible, then compare clusters between loading plots, but there are not enough sample in this dataset#
# Must be done using a genind object. Re-import the data and line it up to 'updated_metadata"
library(hierfstat)
vcf_filtered <- read.vcfR("output/vcf5_filtered_inform.vcf.gz")

gl <- vcfR2genind(vcf_filtered)
metadata <- read_excel("output/8_clust_training_metadata.xlsx")
colnames(metadata)
filtered_metadata <- metadata[metadata$sample_ID %in% indNames(gl), ]
gl_df <- data.frame(genotype = indNames(gl))
joined_data <- gl_df %>%
  left_join(filtered_metadata, by = c("genotype" = "sample_ID"))

strata_data <- joined_data %>%
  dplyr::select(genotype, everything()) %>%
  column_to_rownames(var = "genotype")

strata(gl) <- strata_data

popNames(gl)
setPop(gl) <- ~cluster
gl
ploidy(gl) <- 1
#Get the eigenvalues or make clusters into factors
#eig.val <- dapc$eig
clusters <- dapc_final$assign
gl@pop <- as.factor(clusters)
hf_data <- genind2hierfstat(gl)

#Pairwise test using genetic distance
dist_matrix <- hierfstat::genet.dist(hf_data, diploid = FALSE, method = "Dch")
#Plot as a heatmap
dist_df <- as.data.frame(as.matrix(dist_matrix))
dist_df$Var1 <- rownames(dist_df)

dist_melt <- melt(dist_df, id.vars = "Var1")
pop_map <- c("1" = 6, "2" = 5, "3" = 4, "4" = 7, "5" = 3, "6" = 1, "7" = 2, "8" = 8)

# Apply mapping to Var1 and variable columns
dist_melt$Var1_mapped <- pop_map[as.character(dist_melt$Var1)]
dist_melt$variable_mapped <- pop_map[as.character(dist_melt$variable)]

ggplot(dist_melt, aes(x = Var1, y = variable, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 3)), size = 3) +
  scale_fill_viridis_c(name = "Genetic Distance (Dch)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Genetic Distance Heatmap", x = "Population", y = "Population")

library(reshape2)
#Calculate Weir & Cockerham Fst
clusters <- dapc_final$assign
gl@pop <- as.factor(clusters)
#Converting to a dataframe make it easier to make into the expected haploid format
gl@tab[1:5, 1:10]

df <- genind2df(gl, sep = "", usepop = TRUE)

df[-1] <- lapply(df[-1], function(x) {
  x[x == "00"] <- 1
  x[x == "11"] <- 2
  return(as.numeric(x))
})

hf_data <- df
colnames(hf_data)[1] <- "pop"

#Also assess how this works with the haploid dataset configured as a diploid set
hf_data <- genind2hierfstat(gl)

#Pairwise test using genetic distance Fst and bootsrapping
Fst_matrix <- pairwise.WCfst(hf_data, diploid=FALSE)

#Plot as a heatmap
fst_df <- as.data.frame(as.matrix(Fst_matrix))
fst_df$Var1 <- rownames(fst_df)
fst_melt <- melt(fst_df, id.vars = "Var1")
fst_melt$value[is.na(fst_melt$value)] <- 0

pop_map <- c("1" = 6, "2" = 5, "3" = 4, "4" = 7, "5" = 3, "6" = 1, "7" = 2, "8" = 8)

# Apply mapping to Var1 and variable columns
fst_melt$Var1_mapped <- pop_map[as.character(fst_melt$Var1)]
fst_melt$variable_mapped <- pop_map[as.character(fst_melt$variable)]


ggplot(fst_melt, aes(x = Var1, y = variable, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 3)), size = 3) +
  scale_fill_viridis_c(name = "Fst") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Genetic Distance Heatmap", x = "Population", y = "Population")

#Jaccards Index
library(vegan)

binary_matrix <- tab(gl, freq = FALSE, NA.method = "zero")

# Compute Jaccard distance (1 - similarity)
jaccard_dist <- vegdist(binary_matrix, method = "jaccard")
# Convert to similarity matrix if needed
jaccard_sim <- 1 - as.matrix(jaccard_dist)

# Step 1: Extract population assignments
pop_assignments <- pop(gl)
pop_assignments <- setNames(as.character(pop_assignments), indNames(gl))

# Step 2: Convert Jaccard similarity matrix to a data frame
jaccard_df <- as.data.frame(as.matrix(jaccard_sim))

# Step 3: Add population info
jaccard_df$Pop1 <- pop_assignments[rownames(jaccard_df)]
jaccard_df$Individual1 <- rownames(jaccard_df)

# Step 4: Pivot to long format
jaccard_long <- jaccard_df %>%
  pivot_longer(
    cols = -c(Individual1, Pop1),
    names_to = "Individual2",
    values_to = "JaccardSimilarity"
  )

# Step 5: Add Pop2 using named vector
jaccard_long <- jaccard_long %>%
  mutate(Pop2 = pop_assignments[Individual2]) %>%
  select(Individual1, Individual2, JaccardSimilarity, Pop1, Pop2)

# Step 6: Compute average similarity between populations
average_similarity <- jaccard_long %>%
  group_by(Pop1, Pop2) %>%
  summarise(AverageSimilarity = mean(JaccardSimilarity, na.rm = TRUE)) %>%
  ungroup()

# Step 7: Recode populations to match cluster labels
pop_map <- c("1" = 6, "2" = 5, "3" = 4, "4" = 7, "5" = 3, "6" = 1, "7" = 2, "8" = 8)
average_similarity_f <- average_similarity %>%
  mutate(
    Pop1 = recode(as.character(Pop1), !!!pop_map),
    Pop2 = recode(as.character(Pop2), !!!pop_map),
    Pop1 = as.factor(Pop1),
    Pop2 = as.factor(Pop2)
  )

# Step 8: Save to Excel
write.xlsx(average_similarity, file = "output/average_jaccard_similarity.xlsx", rowNames = TRUE)

# Step 9: Plot heatmap
ggplot(average_similarity_f, aes(x = Pop2, y = Pop1, fill = AverageSimilarity)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(AverageSimilarity, 2)), size = 3) +
  scale_fill_viridis(name = "Average Jaccard Similarity", option = "D") +
  labs(
    title = "Heatmap of Average Jaccard Similarity Between Populations",
    x = "Clusters",
    y = "Clusters"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


#MLG assessment
#First by assigning the MLGs according to the dapc cluster
gen <- read.vcfR('output/vcf5_filtered_inform.vcf.gz')
gen <- vcfR2genind(gen)
metadata <- read_excel("output/8_clust_training_metadata.xlsx")
colnames(metadata)
filtered_metadata <- metadata[metadata$sample_ID %in% indNames(gen), ]
gl_df <- data.frame(genotype = indNames(gen))
joined_data <- gl_df %>%
  left_join(filtered_metadata, by = c("genotype" = "sample_ID"))

strata_data <- joined_data %>%
  dplyr::select(genotype, everything()) %>%
  column_to_rownames(var = "genotype")

strata(gen) <- strata_data

popNames(gen)
setPop(gen) <- ~cluster
gen
# gen.ml_grid <- mlg.id(gen) 
# gen.ml_grid 
# gen@strata
#X <- genclone2genind(gen_obj1)
X <- genind2loci(gen)
#Pairwise genetic distances
pairwise_genetic_distance <- dist.gene(X, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
hist(pairwise_genetic_distance)
# Convert the matrix to a data frame
pairdist_df <- as.data.frame(as.matrix(pairwise_genetic_distance))
#check that the names match
# Step 2: Add population assignments
# Assuming 'gen' is your genind object
# Step 1: Extract population assignments and name them
pop_assignments <- pop(gen)

pop_assignments <- setNames(as.character(pop_assignments), indNames(gen))

str(pop_assignments)

# Step 2: Convert the pairwise distance matrix to a data frame
pairdist_df <- as.data.frame(as.matrix(pairwise_genetic_distance))

# Step 3: Add population info as row and column attributes
pairdist_df$Pop1 <- pop_assignments[rownames(pairdist_df)]
pairdist_df$Individual1 <- rownames(pairdist_df)

# Step 4: Pivot to long format
pairdist_long <- pairdist_df %>%
  pivot_longer(
    cols = -c(Individual1, Pop1),
    names_to = "Individual2",
    values_to = "GeneticDistance"
  )

# Step 5: Add Pop2 using named vector
pairdist_long <- pairdist_long %>%
  mutate(Pop2 = pop_assignments[Individual2])

pairdist_long <- pairdist_long %>%
  select(Individual1, Individual2, GeneticDistance, Pop1, Pop2)
#Fix the problematic NA values
# pairdist_long <- pairdist_long %>%
#   mutate(
#     Pop1 = as.numeric(replace_na(as.character(Pop1), "2")),
#     Pop2 = as.numeric(replace_na(as.character(Pop2), "2"))
#   )


average_distances <- pairdist_long %>%  group_by(Pop1, Pop2) %>%  summarise(AverageDistance = mean(GeneticDistance, na.rm = TRUE)) %>%  ungroup()
#Make the results into a matrix
distance_matrix <- average_distances %>% pivot_wider( names_from = Pop2, values_from = AverageDistance)
#Make it match what I already have in my otehr figures
pop_map <- c( "1" = 6, "2" = 5, "3" = 4, "4" = 7, "5" = 3, "6" = 1, "7" = 2, "8" = 8)

average_distances_f <- average_distances %>%  mutate( Pop1 = recode(as.character(Pop1), !!!pop_map), Pop2 = recode(as.character(Pop2), !!!pop_map)) 
write.xlsx(average_distances, file = "output/average_pairwise_dist.xlsx", rowNames = TRUE)

#Make into a factor for visulaisasion
average_distances_f <- average_distances_f %>% mutate(Pop1 = as.factor(Pop1), Pop2 = as.factor(Pop2))

#Plot

ggplot(average_distances_f, aes(x = Pop2, y = Pop1, fill = AverageDistance)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(AverageDistance, 1)), size = 3) +
  scale_fill_viridis(name = "Average Pairwise Distance", option = "D") +
  labs(
    title = "Heatmap of Average Genetic Distance Between Populations",
    x = "Clusters",
    y = "Clusters"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

library(viridis)

#Getting membership probabilities for the dapcs
dapc_final$posterior

# Define your custom colors
myCol <- c("darkblue", "purple", "darkgreen", "orange", "red", "brown", "magenta", "gold")

# Define your cluster label mapping
pop_map <- c("1" = 6, "2" = 5, "3" = 4, "4" = 7, "5" = 3, "6" = 1, "7" = 2, "8" = 8)

# Convert posterior matrix to data frame
posterior_df <- as.data.frame(dapc_final$posterior)
posterior_df$individual <- rownames(dapc_final$posterior)

# Convert to long format
posterior_long <- posterior_df %>%
  pivot_longer(cols = -individual, names_to = "cluster", values_to = "probability")

# Determine most likely cluster per individual
max_cluster <- posterior_df %>%
  select(-individual) %>%
  apply(1, function(x) names(x)[which.max(x)])

# Apply pop_map to remap cluster labels
remapped_cluster <- pop_map[max_cluster]

# Add remapped cluster info
posterior_long$remapped_cluster <- rep(remapped_cluster, each = 8)

# Order individuals by remapped cluster
posterior_long <- posterior_long %>%
  mutate(individual = factor(individual, levels = unique(individual[order(remapped_cluster)])))

# Create a named vector for legend labels based on pop_map
legend_labels <- paste("Cluster", pop_map)
names(legend_labels) <- names(pop_map)  # original cluster numbers

# Plot using ggplot2 with updated legend labels
ggplot(posterior_long, aes(x = individual, y = probability, fill = cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = myCol, labels = legend_labels) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Individuals", y = "Membership Probability", fill = "Cluster",
       title = "DAPC Posterior Membership Probabilities (Grouped by Remapped Cluster)")





