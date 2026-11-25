# ==============================================================================
# 1. SIMULATE FAMILY PEDIGREE AND GENOTYPES (simfam)
# ==============================================================================
library(simfam)
library(kinship2)
library(popkin)
library(gaston) # Required for BED matrix conversion and saving

set.seed(2026)

# data dimensions
n <- 100 # number of individuals per generation
G <- 3  # number of generations

# draw the random pedigree with those properties
data <- sim_pedigree(n, G)

# save these objects separately
            
fam <- data$fam
ids <- data$ids
kinship_local_G <- data$kinship_local

# Convert the `fam` table into a `pedigree` object required by the `kinship2` package.
obj <- pedigree(
  id = fam$id,
  dadid = fam$pat,
  momid = fam$mat,
  sex = fam$sex,
  affected = fam$pheno
)

# plot the pedigree
plot(obj, cex = 0.7)

m <- 10000
# draw uniform ancestral allele frequencies
p_anc <- runif(m)

# draw independent (unstructured) genotypes for founders
X_1 <- matrix(
  rbinom(n * m, 2, p_anc),
  nrow = m,
  ncol = n
)

# `simfam` requires names for the founders on `X_1`, to validate identities and align
colnames(X_1) <- ids[[1]]
kinship_est_1 <- popkin(X_1)

# true/expected kinship
kinship_1 <- diag(n) / 2
colnames(kinship_1) <- colnames(X_1)
rownames(kinship_1) <- colnames(X_1)

# Propagate genotypes down the pedigree
X <- geno_fam(X_1, fam)

kinship_local = kinship_fam(kinship_1, fam)
kinship_local_est = popkin(X)

svd_total = svd(kinship_local_est)
plot(svd_total$d)



plot_popkin(
  list(kinship_local, kinship_local_est),
  titles = c("Expected", "Estimated"),
  names = TRUE,
  names_cex = 0.8,
  leg_width = 0.2,
  mar = c(3,2)
)

# ==============================================================================
# 2. CONVERT AND SAVE TO PLINK BED FORMAT (gaston)
# ==============================================================================
message("Converting simulated data to BED format...")

# simfam outputs X as (loci x individuals). gaston expects (individuals x loci).
X_t <- t(X)

# Create a gaston bed.matrix object
bm <- as.bed.matrix(X_t)

# Update the pedigree information (.fam data)
# PLINK requires missing parents to be "0" instead of NA
bm@ped$famid  <- fam$fam
bm@ped$id     <- fam$id
bm@ped$father <- ifelse(is.na(fam$pat), "0", fam$pat)
bm@ped$mother <- ifelse(is.na(fam$mat), "0", fam$mat)
bm@ped$sex    <- fam$sex
bm@ped$pheno  <- fam$pheno

# Update the SNP information (.bim data)
# Since these are simulated, we create dummy SNP IDs, positions, and alleles
bm@snps$chr  <- 1
bm@snps$id   <- paste0("SNP", 1:m)
bm@snps$dist <- 0
bm@snps$pos  <- 1:m
bm@snps$A1   <- "A"
bm@snps$A2   <- "G"

# Save to disk
# I am naming it "yri_pruned_subset" so it seamlessly feeds into your next script
output_dir <- "family_data_synthetic"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
file_prefix <- file.path(output_dir, "family_data_syn")
write.bed.matrix(bm, file_prefix)

message(paste("Successfully saved to", file_prefix, "(.bed, .bim, .fam)"))

# ==============================================================================
# 3. DOWNSTREAM ANALYSIS (Your provided snippet)
# ==============================================================================
# source("inference_short.R")
# source("helper.R")

X_chr2 <- geno_fam(X_1, fam)

kinship_local = kinship_fam(kinship_1, fam)
kinship_local_est = popkin(X)

plot_popkin(
  list(kinship_local),
  titles = c("Kinship"),
  names = TRUE,
  names_cex = 0.8,
  leg_width = 0.2,
  mar = c(1,1)
)

# ==============================================================================
# 2. CONVERT AND SAVE TO PLINK BED FORMAT (gaston)
# ==============================================================================
message("Converting simulated data to BED format...")

# simfam outputs X as (loci x individuals). gaston expects (individuals x loci).
X_t <- t(X_chr2)

# Create a gaston bed.matrix object
bm <- as.bed.matrix(X_t)

# Update the pedigree information (.fam data)
# PLINK requires missing parents to be "0" instead of NA
bm@ped$famid  <- fam$fam
bm@ped$id     <- fam$id
bm@ped$father <- ifelse(is.na(fam$pat), "0", fam$pat)
bm@ped$mother <- ifelse(is.na(fam$mat), "0", fam$mat)
bm@ped$sex    <- fam$sex
bm@ped$pheno  <- fam$pheno

# Update the SNP information (.bim data)
# Since these are simulated, we create dummy SNP IDs, positions, and alleles
bm@snps$chr  <- 2
bm@snps$id   <- paste0("SNP", 1:m)
bm@snps$dist <- 0
bm@snps$pos  <- 1:m
bm@snps$A1   <- "A"
bm@snps$A2   <- "G"

# Save to disk
# I am naming it "yri_pruned_subset" so it seamlessly feeds into your next script
file_prefix <- file.path(output_dir, "family_data_syn_chr2")
write.bed.matrix(bm, file_prefix)

message(paste("Successfully saved to", file_prefix, "(.bed, .bim, .fam)"))



message("Loading genotype data...")
file_prefix <- file.path(output_dir, "family_data_syn")
x <- read.bed.matrix(file_prefix)
x <- select.snps(x, maf > 0.05)
