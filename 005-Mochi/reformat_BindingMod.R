#!/usr/bin/env Rscript

###########################
### R packages
###########################

library(data.table)

###########################
### COMMAND-LINE OPTIONS
###########################

###########################
### Functions
###########################

###########################
### Globals
###########################

modweights_file <- "mochi_model_pymochi_modabs/task_1/weights/weights_BindingMod.txt"
fitness_file <- "GD_JUNX.RData"

###########################
### Main
###########################

#Multiplicative term weights
mw_dt <- fread(modweights_file)[order(id!='WT',Pos)]
for(i in 1:10){
  mw_dt[, paste0('mod_term_', i) := abs(.SD[[1]]+mw_dt[id=='WT',.SD[[1]],,.SDcols = paste0("fold_", i)]),,.SDcols = paste0('fold_', i)]
  mw_dt[id_ref=='WT', paste0('mod_term_', i) := .SD[[1]]/2,,.SDcols = paste0('mod_term_', i)]
}
mw_dt[, mod_term_mean := rowMeans(.SD),,.SDcols = paste0('mod_term_', 1:10)]
mw_dt[, mod_term_sd := apply(.SD, 1, sd),,.SDcols = paste0('mod_term_', 1:10)]

#bZ identity
load(fitness_file)
bz_dt <- all_variants[!duplicated(bZ)][order(bZ_aa_seq, decreasing = T),.(bZ, bZ_aa_seq)]
mw_dt[id!="WT", bZ := bz_dt[bZ!="FOS",bZ]]
mw_dt[id=="WT", bZ := "FOS"]

#Save
write.table(mw_dt, file = gsub("BindingMod", "BindingMod_bZ", modweights_file), quote = F, sep = "\t", row.names = F, na = "")


