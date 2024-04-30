################################################################################
# Packages
################################################################################

library(stringr)
library(parallel)

################################################################################
# Functions
################################################################################


one.sample.ttest = function(m,s,mu,n){
  tstat = (m-mu)/s
  2*pt(abs(tstat), n-1, lower.tail=FALSE)
}


################################################################################
# Load data
################################################################################

wt = read.delim("000-data/TableS1_wt_partners_barcodes.txt")
wt = wt[wt$aa.sequence != "" & !(wt$Name %in% c("CEBPG","CREB1")),]


load("005-Mochi/input_Mochi.Rdata")

binding_scores = read.delim("005-Mochi_output/mochi_model_pymochi_modabs/task_1/predictions/predicted_phenotypes_all.txt")
binding_scores$id = all_variants$id[match((binding_scores$aa_seq), all_variants$aa_seq)]

dG = read.delim("005-Mochi_output/mochi_model_pymochi_modabs/task_1/weights/weights_Binding.txt",header=T)

mult = read.delim("005-Mochi_output/mochi_model_pymochi_modabs/task_1/weights/weights_BindingMod_bZ.txt",header=T)


################################################################################
# Merge datasets
################################################################################

# keep only expected variants
binding_scores = binding_scores[substr(binding_scores$id,2,2) %in% c("x",letters[1:7]),]

# convert Mochi's position to bZ nomenclature
pos = dG$Pos
pos[pos < 28 | pos > 59] = NA
hept = floor((pos - 28)/ 7) + 1
hp = letters[(pos - 28) %% 7 + 1]
mut = str_sub(dG$id,-1,-1)
dG$id = paste0(hept,hp,mut)
dG$id[is.na(pos)] = NA
dG$id[1] = "0xx"

dG_mut = dG[!is.na(dG$id) & dG$id != "0xx",]
dG_bZ = tail(dG,51)
dG_bZ$bZ = wt$Name[wt$Name != "FOS"] # Fos doesn't have a ddG since it is the reference. It's ddG can be considered being 0


binding_scores$ddG_mut = dG_mut$mean[match(binding_scores$id, dG_mut$id)]
binding_scores$ddG_mut[binding_scores$id == "0xx"] = 0

binding_scores$ddG_bZ = dG_bZ$mean[match(binding_scores$bZ, dG_bZ$bZ)]
binding_scores$ddG_bZ[binding_scores$bZ == "FOS"] = 0

binding_scores$dG_mult = mult$mod_term_mean[match(binding_scores$bZ, mult$bZ)]

# we don't want to consider low count variants for further statistical analysis
binding_scores$filt = binding_scores$count_e1_s0 <= 10 | binding_scores$count_e2_s0 <= 10 | binding_scores$count_e3_s0 <= 10 | binding_scores$count_e4_s0 <= 10 | binding_scores$count_e5_s0 <= 10 | binding_scores$count_e6_s0 <= 10 | binding_scores$count_e1_s1 == 0 | binding_scores$count_e2_s1 == 0 | binding_scores$count_e3_s1 == 0 | binding_scores$count_e4_s1 == 0 | binding_scores$count_e5_s1 == 0 | binding_scores$count_e6_s1 == 0 


################################################################################
# Calculate relative, global and specific effects
################################################################################

# relative binding scores
wt_fitness = binding_scores$fitness[binding_scores$id == "0xx"]
names(wt_fitness) = binding_scores$bZ[binding_scores$id == "0xx"]

wt_fitness_sigma = binding_scores$sigma[binding_scores$id == "0xx"]
names(wt_fitness_sigma) = binding_scores$bZ[binding_scores$id == "0xx"]

binding_scores$wt_fitness = wt_fitness[binding_scores$bZ]
binding_scores$wt_fitness_sigma = wt_fitness_sigma[binding_scores$bZ]

binding_scores$rel = binding_scores$fitness - binding_scores$wt_fitness

binding_scores$rel_pval = one.sample.ttest(binding_scores$rel, sqrt(binding_scores$sigma^2 + binding_scores$wt_fitness_sigma^2), 0, 6)
binding_scores$rel_pval[binding_scores$id == "0xx" | binding_scores$filt] = NA
binding_scores$rel_fdr = p.adjust(binding_scores$rel_pval, method="fdr")


# global effects
binding_scores$glob = binding_scores$mean - binding_scores$wt_fitness
binding_scores$glob_pval = one.sample.ttest(binding_scores$glob, sqrt(binding_scores$std^2 + binding_scores$wt_fitness_sigma^2), 0, 6)
# here we don't want to consider ATF4 because of the fitting issue
binding_scores$glob_pval[binding_scores$id == "0xx" | binding_scores$bZ == "ATF4" | binding_scores$filt] = NA
binding_scores$glob_fdr = p.adjust(binding_scores$glob_pval, method="fdr")


# specific effects
binding_scores$spe = binding_scores$fitness - binding_scores$mean
binding_scores$spe_pval = one.sample.ttest(binding_scores$spe, sqrt(binding_scores$std^2 + binding_scores$sigma^2), 0, 6)
# here we don't want to consider ATF4 because of the fitting issue
binding_scores$spe_pval[binding_scores$id == "0xx" | binding_scores$bZ == "ATF4" | binding_scores$filt] = NA
binding_scores$spe_fdr = p.adjust(binding_scores$spe_pval, method="fdr")

save(binding_scores,file="000-data/binding_scores.Rdata")
write.table(binding_scores, file="000-data/TableS3_full_dataset.txt", row.names=F, sep="\t")


