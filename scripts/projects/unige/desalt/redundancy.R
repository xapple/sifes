#!/usr/bin/env R

# Lib
library(vegan)

# Explanatory in two parts #
explanatory = read.table('/Users/sinclair/Desktop/explanatory.tsv', header=TRUE, sep='\t', row.names='X')

explanatory$custom_grouping  = as.factor(explanatory$custom_grouping)
explanatory$custom_attribute = as.factor(explanatory$custom_attribute)
explanatory$replicate_id     = as.factor(explanatory$replicate_id)

explanatory_numeric = explanatory[,c("depth", "grain_size", "salinity", "temperature")]
explanatory_struct  = explanatory[,c("custom_grouping", "custom_attribute", "replicate_id")]

# Response and normalization #
response      = read.table('/Users/sinclair/Desktop/response.tsv',    header=TRUE, sep='\t', row.names='X')
response_norm = response / rowSums(response) # Possilbit?? de prendre le log(1+x), Anscombe

# Distance matrices
#response_dist_moris = vegdist(response,      index="morisita")
#response_dist_bray  = vegdist(response_norm, index="bray")

# Test 1: Redundancy analysis #
answer_rda = rda(response_norm, explanatory_numeric)
answer_rda

plot(answer_rda, scaling=1)
plot(answer_rda, scaling=2)

# Test 2: Bioenv procedure #
answer_be = bioenv(response_norm, explanatory_numeric)
answer_be

# Test 3: ANOVA #
answer_anova = adonis(response_norm ~ explanatory$custom_grouping * explanatory$salinity)
#answer_anova = adonis(response_norm ~ explanatory$custom_grouping * explanatory$temperature)