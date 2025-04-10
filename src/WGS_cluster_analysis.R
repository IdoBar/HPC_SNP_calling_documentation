devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")

# CRAN_packages <- c("tidyverse", "RColorBrewer", "ggrepel","GenABEL", "outliers", "SeqArray",
#                    "SNPRelate","adegenet", "SNPassoc", "vcfR", "glue", "paletteer")
CRAN_packages <- c("tidyverse","adegenet", "RColorBrewer",  "dendextend", "poppr","glue", #"ComplexHeatmap",
                   "pvclust", "colorspace", "gplots", "paletteer", "pheatmap", "vcfR", "here", "readxl")#, 
                   # "ggplotify", "cowplot", "pryr") # "GenABEL", "SNPRelate",
# pak::pak(CRAN_packages)
pacman::p_load(char=CRAN_packages, install = FALSE)


#### Read data files ####
# Read sample metadata and save to file
analysis_basename="A_rabiei_2023"
variant_method <- "Freebayes"
analysis_folder <- "../FB_all_22_05_2023"
analysis_outdir <- here(glue("output/{variant_method}"))


# read isolate sequencing info
sequencing_table <- readxl::read_excel(here("sample_info/A_rabiei_isolates_WGS_metadata.xlsx"), sheet = "merged_table") 

discard_samples <- sequencing_table %>% filter(grepl("Do not use|Use merged", Comments)) %>% 
  select(VCF_name) %>%  
  write_csv("sample_info/duplicated_samples.csv", col_names = FALSE) %>% unlist()

# isolate metadata 
samples_table <- readxl::read_excel(here("sample_info/Summary_tables_AR_2013-2022.xlsx"),
                                    "samples_db") #%>% filter(Isolate %in% sequencing_table$Isolate)

# read DArT-based clustering
dart_clusters <- readxl::read_excel(here("../../A_rabiei_DArT/output/A_rabiei_DArT_poppr.xlsx"), 
                                    "MLG_Cluster_table") %>%   select(-Pathotype)


vcf_file <- list.files(file.path(analysis_folder, "Filtered"), "A_rabiei_2018_2020.bwa.fb.stringent.Q20.GT95.noRep.poly.vcf", full.names = TRUE)

vcf <- read.vcfR(vcf_file)
vcf_samples <- colnames(extract.gt(vcf))
# ?vcfR::write.fasta
sum(!vcf_samples %in% discard_samples)
clean_vcf <- vcf#[,!vcf_samples %in% discard_samples]
sample_cols <- colnames(extract.gt(clean_vcf))
sample_names <- sub("_.+?$", "", colnames(extract.gt(clean_vcf)))

# phenotyping scores
raw_disease_scores <- read_excel("sample_info/raw_pheno_summary_scores.xlsx", 
                                 sheet = "isolate_genotype_trait_scores") %>% 
  inner_join(tibble(Isolate=sample_names), .)


# save phenotypes in wide format for GWAS
raw_disease_scores %>% group_by(Isolate, Genotype) %>% 
  summarise(mean_score=mean(mean_score)) %>% 
  pivot_wider(names_from = Genotype, values_from = mean_score) %>% 
  # left_join(tibble(Isolate=indNames(genind_obj)), .) %>% 
  filter(if_all(where(is.double), ~!is.na(.))) %>% 
  setNames(sub(" ", "_", names(.))) %>%
  write_xlsx(., "output/pheno_summary_scores_GWAS.xlsx", 
             sheet = "isolate_mean_scores_wide", overwritesheet=TRUE)

raw_disease_scores %>% group_by(Isolate, Genotype) %>% 
  filter(Trait=="Leaf Score") %>% 
  # summarise(mean_score=mean(mean_score)) %>% 
  pivot_wider(names_from = Genotype, values_from = mean_score) %>% 
  # left_join(tibble(Isolate=indNames(genind_obj)), .) %>% 
  filter(if_all(where(is.double), ~!is.na(.))) %>% 
  setNames(sub(" ", "_", names(.))) %>%
  write_xlsx(., "output/pheno_summary_scores_GWAS.xlsx", 
             sheet = "isolate_mean_leaf_scores_wide", overwritesheet=TRUE)

raw_disease_scores %>% group_by(Isolate, Genotype) %>% 
  filter(Trait=="Stem Score") %>% 
  # summarise(mean_score=mean(mean_score)) %>% 
  pivot_wider(names_from = Genotype, values_from = mean_score) %>% 
  # left_join(tibble(Isolate=indNames(genind_obj)), .) %>% 
  filter(if_all(where(is.double), ~!is.na(.))) %>% 
  setNames(sub(" ", "_", names(.))) %>%
  write_xlsx(., "output/pheno_summary_scores_GWAS.xlsx", 
             sheet = "isolate_mean_stem_scores_wide", overwritesheet=TRUE)


samples_table %>% 
  inner_join(tibble(Isolate=sample_names), .) %>% 
  select(Isolate, Patho_group=Path_rating) %>% 
  filter(!is.na(Patho_group)) %>%
  write_xlsx(., "output/pheno_summary_scores_GWAS.xlsx", 
             sheet = "isolate_Patho_group", overwritesheet=TRUE)

missing_pheno_samples <- tibble(sample_cols) %>% 
  filter(!grepl(paste(unique(raw_disease_scores$Isolate), 
                      collapse = "|"),sample_cols)) %>% 
  write_csv("sample_info/missing_pheno_samples.csv", col_names = FALSE)



#### adegenet ####
genind_obj <- vcfR2genind(clean_vcf) # clean_vcf
ploidy(genind_obj) <- 1
# summary(genind_obj)
indNames(genind_obj) <- sub("_.+?$", "", indNames(genind_obj))
# create a metadata file
samples_strata <- tibble(Isolate=indNames(genind_obj)) %>% 
  left_join(samples_table) %>% # slice(match(indNames(genind_obj), Isolate)) %>% 
  dplyr::select(Isolate, State, Year, Genotype, Path_rating, weighted_score) %>% 
  left_join(dart_clusters, by=c("Isolate"="Ind")) %>% rename(GBS_cluster = Cluster_fact) %>% 
  mutate(State=factor(State, levels = c("QLD", "NSW", "VIC", "SA", "WA")),
         GBS_cluster=factor(GBS_cluster), 
         MLG=factor(as.character(MLG)),
         Host=factor(Genotype), 
         # Year=factor(Year, levels = sort(unique(Year))),
         Patho.Group=factor(paste0("Group", as.character(Path_rating)), levels=paste0("Group", 0:5))) %>%
  mutate_if(is.factor, ~fct_na_value_to_level(., level = "Unknown"))#%>% 
  # mutate_all(~replace_na(..1, "Unknown"))


# save phenotypes in wide format for GWAS
# tibble(Isolate=indNames(genind_obj)) %>% 
#   left_join(samples_table) %>% # slice(match(indNames(genind_obj), Isolate)) %>% 
#   dplyr::select(Isolate, Path_rating, weighted_score) %>% 
#   # left_join(tibble(Isolate=indNames(genind_obj)), .) %>% 
#   filter(is.na(weighted_score))
#   filter(if_all(where(is.double), ~!is.na(.))) %>% 
#   setNames(sub(" ", "_", names(.))) %>%
#   write_xlsx(., "sample_info/raw_pheno_summary_scores.xlsx", 
#              sheet = "isolate_mean_scores_wide", overwritesheet=TRUE)


# Assign samples factors to Strata
genind_obj@strata <- samples_strata %>% as.data.frame() #%>% 
  # mutate(State=factor(State, levels = c("QLD", "NSW", "VIC", "SA", "WA")),
  #        Year=factor(Year, levels = sort(unique(Year))),
  #        Patho.Group=factor(Patho.Group, levels=paste0("Group", 0:5))) %>%
  # mutate_at(vars(Host, Haplotype, Patho.Group), ~fct_explicit_na(., na_level = "Unknown"))
                             # Pathogenicity=factor(Pathogenicity, levels = sort(unique(Pathogenicity))))
# set state as population factor
pop(genind_obj) <- genind_obj@strata$State
# pop(genind_obj) <- samples_table$State[match(indNames(genind_obj), samples_table$Isolate)]

str(strata(genind_obj))

# hist(distgenEUCL)

#### DAPC analysis ####
# adegenet::adegenetTutorial("dapc")
# https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-dapc.pdf
# https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
Arab_dapc <- dapc(genind_obj, var.contrib = TRUE, scale = FALSE, n.pca = 20, n.da = nPop(genind_obj) - 1)
scatter(Arab_dapc, cell = 0, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)
# cross-validating the number of PCs to retain
myInset <- function(dapc_obj, text_x=5, text_y=90, text_cex=0.85){
  # dapc_obj=Year_dapc$DAPC
  temp <- dapc_obj$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc_obj$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated variance (%)",
       cex=1, pch=20, type="h", lwd=2)
  text(text_x, text_y,  sprintf("PCs=%s", dapc_obj$n.pca),
       cex=text_cex, pos=4)
}

# now repeat 100 times
# plot_dapc <- function(gen_obj, strata_var, plot_filename, plot_widt=8, plot_heig=7,
#                       pal_pack="rcartocolor", pal="Bold",
#                       pca_range=6:18, reps=500, shapes=15:19, label_groups=FALSE, leg_pos="topright",
#                       pca_pos="bottomleft", ncores=parallel::detectCores()-1){
#   # calculate optimal number of PCs to retain
#   Arabx <- xvalDapc(tab(gen_obj, NA.method = "mean"), strata(gen_obj)[[strata_var]],
#                         n.pca = pca_range, n.rep = reps,
#                         parallel="snow", ncpus=ncores)
# 
#   # save plot
#   pdf(file = plot_filename, width = plot_widt, height=plot_heig)
#   scatter(Arabx$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
#           col = paletteer_d(!!pal_pack, !!pal), scree.da=FALSE, 
#           clabel = label_groups, posi.leg = leg_pos, scree.pca = FALSE, 
#           cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)
#   add.scatter(myInset(Arabx$DAPC), posi=pca_pos,
#               inset=c(-0.03,-0.01), ratio=.25,
#               bg=transp("white"))
#   dev.off()
#   return(invisible(Arabx))
# }

# set constants
pca_range=6:18
reps=1000
shapes=c(15:18, 9,25)
ncores=parallel::detectCores()-1
# cross-validating the number of PCs to retain (rep=1000)
# set consistent colours

favourite_pals <- setNames(c("awtools", "rcartocolor", "RColorBrewer", "ggsci"), 
                           c("mpalette", "Bold", "Set1", "category10_d3"))  
my_colours = list(
  State = levels(samples_strata$State) %>% setNames(c(as.character(paletteer_d("RColorBrewer::Set1", length(.)-1)),"grey1") , .),
  Year = sort(unique(samples_strata$Year)) %>% setNames(as.character(paletteer_d("ggsci::default_uchicago", length(.))), .),
  Patho.Group = levels(samples_strata$Patho.Group) %>% 
    setNames(c(as.character(paletteer_d("RColorBrewer::RdYlGn", direction = -1))[round(seq(4,11, length.out = length(.)-1) ,0)], "grey1"), .),
  #Host = levels(samples_strata$Host) %>% setNames(as.character(paletteer_d("RColorBrewer::Dark2", length(.))), .),
  GBS_cluster = levels(samples_strata$GBS_cluster) %>% 
    setNames(c(as.character(paletteer_d("rcartocolor::Bold", n=length(.)-1)), "grey1"), .)
)

# Cluster by year
strata_var <- "Year"
cols <- my_colours[[strata_var]]
# pal <- "Set1"
Year_dapc <- xvalDapc(tab(genind_obj, NA.method = "mean"), strata(genind_obj)[[strata_var]],
                  n.pca = pca_range, n.rep = reps,
                  parallel="snow", ncpus=ncores)

# save plot
# pdf(file = filedate(sprintf("%s_DAPC_%s_analysis", variant_method,  strata_var), 
#                       ".pdf", outdir = glue("{analysis_outdir}/plots")), width = 7, height=6)
scatter(Year_dapc$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
        col = my_colours[[strata_var]], scree.da=FALSE, # paletteer_d(!!favourite_pals[[pal]], !!pal)
        clabel = FALSE, posi.leg = "topright", scree.pca = FALSE, 
        cleg = 0.85, xax = 1, yax = 2, inset.solid = 1)
add.scatter(myInset(Year_dapc$DAPC), posi="bottomleft",
            inset=c(0.01,-0.01), ratio=.2,
            bg=transp("white"))
# dev.off()

# find unique alleles contributing to the variance of the year
contrib <- loadingplot(Year_dapc$DAPC$var.contr, axis=2,
                       thres=.005, lab.jitter=10)
contrib$var.names # marker names above the threshold

# cluster by state
strata_var <- "State"
# pal <- "Bold"
state_dapc <- xvalDapc(tab(genind_obj, NA.method = "mean"), strata(genind_obj)[[strata_var]],
                          n.pca = pca_range, n.rep = reps,
                          parallel="snow", ncpus=ncores)

# save plot
pdf(file = filedate(sprintf("%s_DAPC_%s_analysis", variant_method,  strata_var), 
                    ".pdf", outdir = glue("{analysis_outdir}/plots")), width = 7, height=6)
scatter(state_dapc$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
        col = my_colours[[strata_var]], scree.da=FALSE, 
        clabel = FALSE, posi.leg = "topright", scree.pca = FALSE, 
        cleg = 0.85, xax = 1, yax = 2, inset.solid = 1)
add.scatter(myInset(state_dapc$DAPC), posi="bottomleft",
            inset=c(0.01,-0.01), ratio=.2,
            bg=transp("white"))

dev.off()
  
  
# cowplot::plot_grid(p1.year,p1.state ,labels = c('A', 'B' ), label_size = 12)
# cluster by Host
strata_var <- "Host"
# pal <- "category10_d3"
host_dapc <- xvalDapc(tab(genind_obj, NA.method = "mean"), strata(genind_obj)[[strata_var]],
                      n.pca = pca_range, n.rep = reps,
                      parallel="snow", ncpus=ncores)

# save plot
pdf(file = filedate(sprintf("%s_DAPC_%s_analysis", variant_method,  strata_var), 
                    ".pdf", outdir = glue("{analysis_outdir}/plots")), width = 7, height=6)
scatter(host_dapc$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
        col = my_colours[[strata_var]], scree.da=FALSE, 
        clabel = FALSE, posi.leg = "topright", scree.pca = FALSE, 
        cleg = 0.85, xax = 1, yax = 2, inset.solid = 1)
add.scatter(myInset(host_dapc$DAPC), posi="bottomright",
            inset=c(0.01,-0.01), ratio=.15,
            bg=transp("white"))
dev.off()

# plot_grid(p1.base,p2.host ,labels = c('A', 'B' ), label_size = 12)  
  
# cluster by Pathogenicity
strata_var <- "Patho.Group"
# pal <- "mpalette"
patho_dapc <- xvalDapc(tab(genind_obj, NA.method = "mean"), strata(genind_obj)[[strata_var]],
                       n.pca = pca_range, n.rep = reps,
                       parallel="snow", ncpus=ncores)

# save plot
pdf(file = filedate(sprintf("%s_DAPC_%s_analysis", variant_method,  strata_var), 
                    ".pdf", outdir = glue("{analysis_outdir}/plots")), width = 8, height=7)
scatter(patho_dapc$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
        col = my_colours[[strata_var]], scree.da=FALSE, solid = 1,
        clabel = FALSE, posi.leg = "topright", scree.pca = FALSE, 
        cleg = 0.85, xax = 1, yax = 2, inset.solid = 1)
add.scatter(myInset(patho_dapc$DAPC), posi="bottomleft",
            inset=c(0.01,-0.03), ratio=.15,
            bg=transp("white"))
dev.off()
# DAPC by GBS clusters
strata_var <- "GBS_cluster"
# pal <- "Bold"
cluster_dapc <- xvalDapc(tab(genind_obj, NA.method = "mean"), genind_obj@strata[[strata_var]],
                      n.pca = pca_range, n.rep = reps,
                      parallel="snow", ncpus=ncores)

# save plot
pdf(file = filedate(sprintf("%s_DAPC_%s_analysis", variant_method,  strata_var), 
                    ".pdf", outdir = glue("{analysis_outdir}/plots")), width = 7, height=6)
scatter(cluster_dapc$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
        col = my_colours[[strata_var]], scree.da=FALSE, 
        clabel = FALSE, posi.leg = "bottomright", scree.pca = FALSE, 
        cleg = 0.85, xax = 1, yax = 2, inset.solid = 1)
add.scatter(myInset(cluster_dapc$DAPC), posi="bottomleft",
            inset=c(0.01,-0.01), ratio=.25,
            bg=transp("white"))
dev.off()

save.image(filedate(sprintf("%s_DAPC_analysis", variant_method), 
                    ".RData", outdir = analysis_outdir))

# Check if distance is different by cluster
# calculate AMOVA (with clone correction)
agc <- as.genclone(genind_obj)
# agc@strata <- strata_factors %>% filter(id %in% indNames(agc)) %>% mutate_if(is.factor, ~fct_drop(.))
amova.result <- poppr.amova(agc, ~GBS_cluster, clonecorrect = TRUE)
amova.cc.test <- randtest(amova.result)
plot(amova.cc.test)
amova.cc.test$pvalue

#### Manual Clustering ####
# Heatmap ####
# correct clone
Arab_cc <- clonecorrect(genind_obj, strata = NA)#, keep = 1:2)
Arab_gl <- dartR::gi2gl(Arab_cc)

strata(Arab_gl) <- genind_obj@strata
# indNames(Arab_gl)
ploidy(Arab_gl) <- 1

# export distance for GWAS
# calculate distance
bitwise.euc <- function(x) bitwise.dist(x, euclidean = TRUE, missing_match = FALSE)
gwas_gl_mat <- Arab_gl %>% dartR::gl.keep.ind(ind.list = unique(raw_disease_scores$Isolate), 
                                          recalc = TRUE, mono.rm = TRUE, verbose = 3) %>% 
  bitwise.dist(Arab_gl, euclidean = TRUE, missing_match = FALSE,  mat = TRUE)
write.table(gwas_gl_mat, "output/gwas_samples_dist_mat.txt", sep = "\t")  
Arab_cc.mat <- 

# pop(gen_clone)
# calculate distance
bitwise.euc <- function(x) bitwise.dist(x, euclidean = TRUE, missing_match = FALSE)
Arab_cc.mat <- bitwise.dist(Arab_gl, euclidean = TRUE, missing_match = FALSE,  mat = TRUE)
Arab_cc.dist <- bitwise.dist(Arab_gl, euclidean = TRUE, missing_match = FALSE, mat = FALSE)
# test <- Arab_gl %>% aboot(cutoff = 50, quiet = TRUE, sample = 10, distance = bitwise.euc)
pop_dist <- Arab_gl   %>%
  # genind2genpop(pop = ~pop) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = bitwise.euc)
# pop_tree_data <- as_tibble(tree) %>% left_join(strata_factors, by = c("label"="id")) %>% 
#   as.treedata()
pop_tree_data <- as_tibble(pop_tree) %>% 
  left_join(strata_factors %>% group_by(pop) %>% 
              slice(1), by = c("label"="pop")) %>% 
  # mutate(id=label, label=pop) %>% 
  as.treedata()


heatmap_metadata <- Arab_gl@strata %>% 
  mutate(Year=factor(Year, levels = sort(unique(Year)))) %>% 
  mutate_if(is.character, ~factor(replace_na(.x, "NA"))) %>% 
column_to_rownames("Isolate")
# data.frame(, row.names = sample_cols) %>% 
   
# row.names(heatmap_metadata) <- sample_cols
# specify colours
# ann_colors = list(
#   State = setNames(paletteer_d(!!favourite_pals[3]::!!names(favourite_pals[3]), n=4), 
#                    levels(heatmap_metadata$State)),
#   Host = setNames(paletteer_d(!!favourite_pals[4]::!!names(favourite_pals[4]), n=4), 
#                   levels(heatmap_metadata$Host)),
#   Year = 
# )

# specify colours
# my_colours = list(
#   State = levels(heatmap_metadata$State) %>% 
#     setNames(as.character(paletteer_d("RColorBrewer::Set1", length(.))), .),
#   Year = levels(heatmap_metadata$Year) %>% 
#     setNames(as.character(paletteer_d("ggsci::default_uchicago", length(.))), .),
#   Pathogenicity = levels(heatmap_metadata$Pathogenicity) %>% 
#     setNames(as.character(paletteer_d("RColorBrewer::RdYlGn", length(.), direction = -1)), .),
#   Host = levels(heatmap_metadata$Host) %>% 
#     setNames(as.character(paletteer_d(glue::glue("{favourite_pals[4]}::{names(favourite_pals[4])}"), 
#                                            n=length(.))), .), 
#   GBS_cluster = levels(heatmap_metadata$GBS_cluster) %>% 
#   setNames(c(as.character(paletteer_d(glue::glue("{favourite_pals[2]}::{names(favourite_pals[2])}"),
#                                     n=length(.)-1)), "gray80"), .)
# )
cluster_num <- 6

Heatmap(Arab_cc.mat, col = paletteer_c("viridis::viridis", 50, direction = -1), name = "Genetic distance", 
            )
dist_heatmap <- pheatmap(Arab_cc.mat, show_rownames = FALSE, show_colnames = FALSE, 
                         color = paletteer_c("viridis::viridis", 50, direction = -1), 
                         annotation_colors = my_colours,
                         # clustering_distance_rows = dist_matrix, clustering_distance_cols = dist_matrix,
                         annotation_row = heatmap_metadata[c("State")], #"Year", "Host"
                         annotation_col = heatmap_metadata[c("Patho.Group")], # , "GBS_cluster"
                         cutree_rows = cluster_num, cutree_cols = cluster_num , 
             filename = filedate(glue("A_rabiei_WGS_heatmap_{cluster_num}clusters"),ext="pdf", outdir = glue("{analysis_outdir}/plots/")),
                         width = 8, height = 6.5)


dist_heatmap <- pheatmap(Arab_cc.mat, annotation_col = heatmap_metadata[c("Host", "State")], cutree_rows = 2,
         cutree_cols = 2, annotation_colors = ann_colors,
         annotation_row = heatmap_metadata[c("Patho.Group", "Year")])
plot(hclust(Arab_cc.dist))
# distgenEUCL <- dist(genind_obj, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)




# clustering with pvclust
ibs.pv <- pvclust(Arab_cc.mat, nboot=500, parallel=TRUE, )
plot(ibs.pv)
ibs.clust <- pvpick(ibs.pv, alpha=0.99)
# Create a heatmap using the same data
Colv  <- as.dendrogram(Arab_cc.dist) %>%
  branches_attr_by_labels(ibs_meta$Sample, TF_values = c(2, Inf), attr = c("lwd")) %>%
  ladderize




X <- tab(Arab_cc, freq = TRUE, NA.method = "mean")
pca1 <- ade4::dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 20)
# define pathotypes
path_levels <- c("Low","Medium","Moderate","High","Very High", "Extreme") %>% factor(., levels = .)
pca_data <- pca1$li %>% rownames_to_column("Isolate")  %>%
  inner_join(samples_table) %>% # sequenced_isolates
  mutate(seq_fontface=if_else(Sequenced=="Illumina", "plain", "bold"),
         Pathotype=factor(Pathotype, levels=path_levels)) %>% #,
  # Coverage=round(Coverage, 2)) %>%
  left_join(sequencing_table %>% filter(Sequencing_Centre=="Macrogen") %>% 
              dplyr::select(Isolate, Sequencing_Centre)) %>% 
  mutate(Sequencing_Centre=if_else(is.na(Sequencing_Centre), "AG", Sequencing_Centre)) %>% 
  arrange(Year)





#### SNPRelate ####
# set thresholds
call_rate = 1
LD_rate = 0.5
# convert, save in "tmp.gds" with the default lzma compression algorithm
gds_file <- sprintf("output/data_files/%s_core.gds", gsub("[\\./]", "",analysis_folder))
seqVCF2GDS(vcf_file, gds_file)
# open a GDS file
genofile <- SeqArray::seqOpen(gds_file)
# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=LD_rate, autosome.only = FALSE, 
                          missing.rate = (1-call_rate))

#### IBS ####
ibs <- snpgdsIBS(genofile, num.thread=2, autosome.only = FALSE, missing.rate = (1-call_rate))
# snprelate_ibs <- snpgdsIBS(genofile, autosome.only = FALSE, maf = maf, missing.rate = (1-call_rate), sample.id = idnames(data1), snp.id = snpnames(data1))
pv_sig=0.95
pv_boot=1000



# ibs_mat <- snprelate_ibs$ibs
ibs_title <- "A. rabiei samples IBS-based clustering"
dendro_file <- filedate(sprintf("%s_pvclust%s_dendrogram",  variant_method,  pv_sig ), 
                        ".pdf", outdir = glue::glue("{analysis_outdir}/plots"))
mat_file <- filedate(sprintf("%s_pvclust%s_heatmap", variant_method,  pv_sig ), 
                     ".pdf", outdir = glue::glue("{analysis_outdir}/plots"))

ids <- read.gdsn(index.gdsn(genofile, "sample.id"))

data1_ibs_table <- process_ibs(ibs_mat = ibs$ibs,ids = idnames(data1), sig_alpha = pv_sig, bootNum = pv_boot, plotTitle = ibs_title, outDendroFile = dendro_file, rotateDendro = TRUE, outMatFile = mat_file)

# Process IBS matrix (make dendrogram and similarity matrix)
process_ibs <- function(ibs_mat, ids, bootNum=5000, sig_alpha=0.95, outDendroFile=FALSE, rotateDendro = TRUE, outMatFile=FALSE, clust_pal = "Paired", plotTitle, dendWidth=15, dendHeight=8, matWidth=12, matHeight=12, dendDimUnit="in",matDimUnit="in", plotRes=100){
  dimnames(ibs_mat) <- list(ids, ids)
  # clustering with pvclust
  ibs.pv <- pvclust(ibs_mat, nboot=bootNum, parallel=TRUE)
  plot(ibs.pv)
  ibs.clust <- pvpick(ibs.pv, alpha=sig_alpha)
  
  names(ibs.clust$clusters) <- paste0("Cluster", 1:length(ibs.clust$clusters))
  # Choose a colour palette
  pal <- brewer.pal(length(ibs.clust$clusters), clust_pal)
  # Transform the list to a dataframe
  ibs_meta <- bind_rows(lapply(names(ibs.clust$clusters),
                               function(l) data.frame(Cluster=l, Sample = ibs.clust$clusters[[l]])))
  if (length(ids)>nrow(ibs_meta)){
    non_clustered <- data.frame(Cluster = "Cluster0", Sample = ids[!ids %in% ibs_meta$Sample])
  } else non_clustered <- NULL
  # Add the rest of the non-clustered samples (and assign them as Cluster0), add colour to each cluster
  ibs_table <- ibs_meta %>%
    rbind(., non_clustered) %>%
    mutate(Cluster_int=as.numeric(sub("Cluster", "", Cluster))) %>%
    mutate(Cluster_col=ifelse(Cluster_int==0, "#000000", pal[Cluster_int])) %>% # Change to #000000 for black
    .[match(ibs.pv$hclust$labels[ibs.pv$hclust$order], .$Sample),]
  
  # Create a plot and save to file
  if (!outDendroFile==FALSE) {
    dendFmt <- sub("^.+\\.(\\w+{3})$", "\\1", outDendroFile)
    if (dendDimUnit=="in" &&  tolower(dendFmt)=="png") {
      dendWidth <- dendWidth*plotRes
      dendHeight <- dendHeight*plotRes
      dendDimUnit <- "px"
    }
    dendCmd <- sprintf("%s('%s', width = %d, height = %d, units = '%s')", 
                       tolower(dendFmt), outDendroFile, dendWidth,dendHeight, tolower(dendDimUnit) )
    eval(parse(text = dendCmd))
  } 
  hcd <- as.dendrogram(ibs.pv) %>%
    #pvclust_show_signif(ibs.pv, show_type = "lwd", signif_value = c(2, 1),alpha=0.25) %>%
    set("leaves_pch", ifelse(ibs_table$Cluster_int>0,19,18)) %>%  # node point type
    set("leaves_cex", 1) %>%  # node point size
    #set("leaves_cex", ifelse(ibs_table$Cluster_int>0,1,0.75)) %>%  # node point size
    set("leaves_col", ibs_table$Cluster_col) %>% #node point color
    branches_attr_by_labels(ibs_meta$Sample, TF_values = c(2, Inf), attr = c("lwd")) %>% # change branch width
    # rect.dendrogram(k=12, cluster = ibs_table$Cluster_int, border = 8, lty = 5, lwd = 1.5,
    #                 lower_rect = 0) %>%  # add rectangles around clusters
    plot(main=plotTitle)
  # pvrect(ibs.pv, alpha=0.95, lwd=1.5)
  if (!outDendroFile==FALSE) {
    dev.off()
    if (rotateDendro) imageFlip(pngName = outDendroFile, plotFormat = dendFmt)
  } 
  # Create a heatmap using the same data
  Colv  <- as.dendrogram(ibs.pv) %>%
    branches_attr_by_labels(ibs_meta$Sample, TF_values = c(2, Inf), attr = c("lwd")) %>%
    ladderize
  # creates my own color palette from red to green
  my_palette <- colorRampPalette(c("green", "black", "red"))(n = 299)
  heatFmt <- sub("^.+\\.(\\w+{3})$", "\\1", outMatFile)
  
  if (!outMatFile==FALSE) {
    if (matDimUnit=="in" &&  tolower(heatFmt)=="png") {
      matWidth <- matWidth*plotRes
      matHeight <- matHeight*plotRes
      matDimUnit <- "px"
    }
    heatCmd <- sprintf("%s('%s', width = %d, height = %d, units = '%s')", 
                       tolower(heatFmt), outMatFile, matWidth,matHeight, tolower(matDimUnit) )
    eval(parse(text = heatCmd))
  } 
  gplots::heatmap.2(ibs_mat, dendrogram = "col", Rowv=Colv, Colv=Colv,trace='none',
                    ColSideColors=sub("#000000", "ivory2", 
                                      ibs_table$Cluster_col[match(idnames(data_clean), ibs_table$Sample)]),  # pal[gr.row]
                    cexRow=1, sepcolor="black", sepwidth=c(0.01,0.01), keysize = 0.8,
                    labCol = FALSE, density.info="none", col = my_palette, scale = "row")
  if (!outMatFile==FALSE)  dev.off()
  return(ibs_table)
}


snpgdsClose(genofile)
