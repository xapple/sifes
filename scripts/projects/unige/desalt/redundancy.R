#!/usr/bin/env R

# Lib
library(vegan)
library(phyloseq)
library(metagenomeSeq)

# Expl
explanatory = read.table('/Users/sinclair/Desktop/explanatory.tsv', header=TRUE, sep='\t', row.names='X')

explanatory$custom_grouping  <- as.factor(explanatory$custom_grouping)
explanatory$custom_attribute <- as.factor(explanatory$custom_attribute)
explanatory$replicate_id     <- as.factor(explanatory$replicate_id)

explanatory_numeric = explanatory[,c("depth", "grain_size", "salinity", "temperature")]
explanatory_struct = explanatory[,c("custom_grouping", "custom_attribute", "replicate_id")]

# Resp and normalization
response      = read.table('/Users/sinclair/Desktop/response.tsv',    header=TRUE, sep='\t', row.names='X')
response_norm = response / rowSums(response)

response_dist_moris = vegdist(response,      index="morisita")
response_dist_bray  = vegdist(response_norm, index="bray")

# Test 1
answer_rda = rda(response_norm, explanatory, index='morisita')

plot(answer, scaling=1)
plot(answer, scaling=2)

# Test 2
bioe_rda = bioenv(response_dist_bray, explanatory_numeric, index='morisita')










### normalization with metagenomeSeq package
# building a phyloseq object

# molecular data
otu_to_norm  <- otu
comp_to_norm <- comp
samples_names <- comp[,s]

comp_ <- comp_to_norm[,c("Locality", "Station", "Grab")]
comp_$Station <- as.factor(comp_$Station)
comp_$Grab <- as.factor(comp_$Grab)
rownames(comp_) <- samples_names
rownames(otu_to_norm) <- samples_names

obj <- phyloseq(otu_table(otu_to_norm, taxa_are_rows = F), sample_data(comp_))
obj_m <- phyloseq_to_metagenomeSeq(obj)
p = cumNormStatFast(obj_m)
dat_mol_norm = cumNorm(obj_m, p = p)
otu_NORM <- t(MRcounts(dat_mol_norm, norm = TRUE, log = TRUE))