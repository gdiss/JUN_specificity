library(parallel)
library(stringdist)
library(ShortRead)
options(stringsAsFactors = FALSE) 


load("000-data/true_barcodes.Rdata")
wt = read.delim("000-data/TableS1_wt_partners_barcodes.txt")
# this list contains the 54 wt bZIPs and some negative controls that we added. We map reads to them but they will evenutally be filtered out because a preliminary analysis showed no signal


datapath = "/tungstenfs/groups/gbioinfo/seqdata"
samples <- dir(datapath)[grep("2935F",dir(datapath))]
samples = samples[!grepl("Undetermined",samples)]
samples = samples[grep("fastq.gz",samples)]
samples_jun1 = samples[c(1,5:12,2:4)]

# we sequenced the output samples twice
# we will sum up the reads
samples <- dir(datapath)[grep("2953F",dir(datapath))]
samples = samples[!grepl("Undetermined",samples)]
samples_jun2 = samples[grep("fastq.gz",samples)]


all_samples = file.path(datapath, c(samples_jun1, samples_jun2))



################################################################################
# Process Fastq files
################################################################################

# We are quite lenient with quality filtering, because we will anyway filter reads that don't correspond to an expected barcode

reads = mclapply(all_samples,function(x){
  
  bar = readFastq(x)
  
  # filter reads with average quality < 20
  scores_bar = alphabetScore(bar)
  length_bar = width(bar)
  avQ_bar = scores_bar / length_bar
  bar = as.vector(bar@sread)
  bar = bar[avQ_bar >= 20]
  
  # count identical reads
  bar = as.data.frame(table(bar))
  bar = bar[order(-bar$Freq),]
  bar$bar = as.character(bar$bar)
  
  bar
  
}, mc.cores=18)

save(reads, file="000-data/raw_counts.Rdata")


# split into barcode regions
reads = mclapply(reads, function(x){
  
  spl = data.frame(
    substr(x$bar, 1, 20),
    substr(x$bar, 21, 26),
    substr(x$bar, 27, 50),
    x$Freq
  )
  
  # aggregate counts of the same barcodes pair (i.e. different sequences in restriction site)
  agg = aggregate(spl[,4], list(spl[,1],spl[,3]), sum)
  names(agg) = c("bcWT", "bcVar", "nReads")
  
  agg
  
},mc.cores=18)




################################################################################
# Map barcode to variants and wt
################################################################################

### WT barcodes first

min(stringdistmatrix(wt$Barcode,method="hamming"))
# the minimum distance between wt barcodes is 6. It is not possible for a read with 2 mutations or less to be amiguous, i.e. getting closer than 4bp from another barcode
# we therefore can aggregate counts from reads within a distance of 2 from the target barcodes


# unique barcode reads
wt_pool = unique(unlist(lapply(reads,function(x){
  unique(x$bcWT)
})))

# Measure Hamming distance with all 54 wild-type barcode sequences
h = do.call("cbind",mclapply(1:nrow(wt),function(x){
  stringdist(wt$Barcode[x],wt_pool,method="hamming", nthread = 1)
},mc.cores=72))

# keep Hamming distance and id to the closest of the 54 bZ + negative controls
# first get their indexes
wt_barcodes = do.call("rbind",mclapply(1:nrow(h),function(x){
  x = h[x,]
  w = which.min(x)
  c(wt$Name[w], x[w])
},mc.cores=24))
wt_barcodes = as.data.frame(cbind(wt_pool,wt_barcodes))
names(wt_barcodes) = c("bcWT","bZ","dist")
wt_barcodes$dist = as.numeric(wt_barcodes$dist)
wt_barcodes = wt_barcodes[wt_barcodes$dist <= 2,]  

rm(h)



### JUN barcodes

# Hamming distance will be measured with a Perl script. We first split the pool of unique barcode reads in 24 for parallelization
Ncore = 72
jun_pool = unique(unlist(lapply(reads,function(x){
  unique(x$bcVar)
})))
sp = split(jun_pool, cut(1:length(jun_pool), breaks=Ncore, include.lowest = T))
fnames = paste0("tmp/var_bc_jun",1:Ncore,".txt")
for(i in 1:Ncore){
  write(sp[[i]],fnames[i])
}

write(true_barcodes$bc, "tmp/true_barcodes.txt")
mclapply(1:Ncore,function(i){
  system(paste("perl 003-dist_var_barcodes.pl", fnames[i], "tmp/true_barcodes.txt"))
},mc.cores=Ncore)


jun_barcodes = do.call("rbind",mclapply(1:Ncore,function(i){
  
  read.delim(paste0(fnames,".out")[i], header=T)
  
},mc.cores=Ncore))

jun_barcodes = data.frame(bcVar = jun_pool, jun_barcodes)

# we keep only those that have a distance of 2 or less to one true barcode and the distance to the second closest is more than 2 plus the distance to the closest
jun_barcodes = jun_barcodes[jun_barcodes$dist1 <= 2 & jun_barcodes$dist2 - jun_barcodes$dist1 > 2, ]



### aggreate reads from barcode pairs with sequencing errors
# all possible pairs of JUN and wt barcodes
exp_jun = as.data.frame(as.matrix(expand.grid(var_id = c(0,unique(jun_barcodes$idx_true_bc_1)), wt_id = c(unique(wt$Name),"no match"))))
exp_jun$var_id = as.numeric(exp_jun$var_id)
exp2 = paste(exp_jun$var_id, exp_jun$wt_id)
bc_counts = do.call("cbind",mclapply(reads, function(x){
  
  # link to prey barcode
  x$id_var = jun_barcodes$idx_true_bc_1[match(x$bcVar, jun_barcodes$bcVar)]
  x$id_wt = wt_barcodes$bZ[match(x$bcWT, wt_barcodes$bcWT)]
  
  # non-recognized barcodes
  x$id_var[is.na(x$id_var)] = 0
  x$id_wt[is.na(x$id_wt)] = "no match"
  
  x = aggregate(x$nReads, list(x$id_var, x$id_wt), sum)
  names(x) = c("prey_id", "wt_id", "nReads")
  
  # order counts from all samples in the same way (defined in exp_jun and exp2)
  n = x$nReads[match(exp2, paste(x$prey_id, x$wt_id))]
  n[is.na(n)] = 0
  
  n
}, mc.cores=18))


save(bc_counts, file="000-data/bc_counts.Rdata")
save(exp_jun, file="000-data/bc_counts_annotation.Rdata") # var_id corresponds to row number in the true_barcodes object, ordered by decreasing Nbc


################################################################################
# Compute binding scores
################################################################################

### First, aggregate counts of the same amino acid variant - wt pairs

id = c("no match", true_barcodes$id)
id = id[exp_jun$var_id + 1] # add one because index 1 is now no match
jun_agg = aggregate(bc_counts, by = list(id, exp_jun$wt_id), sum)

# we then add the re-sequenced output counts
jun_agg2 = jun_agg[,1:14]
jun_agg2[,9:14] = jun_agg2[,9:14] + jun_agg[,15:20]
names(jun_agg2) = c("id","bZ","i1","i2","i3","i4","i5","i6","o1","o2","o3","o4","o5","o6")


### Format for dimsum
# For dimsum, we need to input a nt sequence, not just an id. We don't really care what this sequence is, we can just input a random nt sequence that will serve as an ID. We can then identify back the variant or partner from this sequence/ID when we import the output from dimsum

# first the wt random sequences. There are 73 wt sequences, so a 4nt long sequence will be enough
wt = rbind(wt,c("no match","","",""))
nt = c("A","T","G","C")
wt$seqID = do.call("paste0",as.data.frame(expand.grid(nt,nt,nt,nt)))[1:nrow(wt)]


# for variants we can use the actual variant sequence
var_jun = true_barcodes[!duplicated(true_barcodes$id),]

# we need to pad or strip the shorter or longer sequences
nc = nchar(var_jun$va)
var_jun$va[nc == 199] = substr(var_jun$va[nc==199],1,198)
pad = sapply(198-nc[nc<198],function(x){
  paste(rep("A",x),collapse="")
})
var_jun$va[nc < 198] = paste0(var_jun$va[nc < 198], pad)


tmp1 = var_jun$va[match(jun_agg2$id, var_jun$id)]
tmp1[is.na(tmp1)] = paste(rep("A",198),collapse="") # this is the no match
tmp2 = wt$seqID[match(jun_agg2$bZ, wt$Name)] 
jun_agg2$nt_seq = paste0(tmp1,tmp2)

out = jun_agg2[,c(15,3:14)]
write.table(out, file="000-data/dimsum_input.txt",sep="\t",row.names=F)
save(jun_agg2, file="000-data/counts_agg_by_id.Rdata")




################################################################################
# Run DimSum
################################################################################

# type ./004-dimsum.sh in terminal



################################################################################
# Run Mochi
################################################################################

load("004-dimsum_output/004-dimsum_output_fitness_replicates.RData")
all_variants$id = jun_agg2$id[match(toupper(all_variants$nt_seq), jun_agg2$nt_seq)]
all_variants$bZ = jun_agg2$bZ[match(toupper(all_variants$nt_seq), jun_agg2$nt_seq)]

# filter out variants with more than 3 mutations
n = unlist(lapply(strsplit(all_variants$id,","),length))
all_variants = all_variants[n <= 3,]

# to run Mochi, we keep only the 52 wt bZIPs and remove CEBPG and CREB1 (because they somehow dropped out) and all negative controls
all_variants = all_variants[all_variants$id != "no match" & !(all_variants$bZ %in% c("CEBPG","CREB1")),]
all_variants = all_variants[all_variants$bZ %in% wt$Name[wt$aa.sequence != ""],]
all_variants$var_aa_seq = true_barcodes$va_aa[match(all_variants$id, true_barcodes$id)]
all_variants$var_aa_seq = substr(all_variants$var_aa_seq, 1, 65)

# dummy encoding of WT bZIPs in aa sequence. FOS is encoded as repeatrs of A. Each other partner is encoded as a mutation towards C at a different position
wt_bZ = wt$Name[wt$aa.sequence != "" & !(wt$Name %in% c("CEBPG","CREB1"))]
fos_aa = paste(rep("A",51),collapse="")
wt_aa = sapply(1:51,function(i){
  substr(fos_aa,i,i) = "C"
  fos_aa
})
w = which(wt_bZ == "FOS")
wt_aa = c(wt_aa[1:(w-1)],fos_aa,wt_aa[w:51])
names(wt_aa) = wt_bZ

all_variants$bZ_aa_seq = wt_aa[all_variants$bZ]
all_variants$aa_seq = paste0(all_variants$var_aa_seq, all_variants$bZ_aa_seq)

save(all_variants, file="005-Mochi/input_Mochi.Rdata")




# Mochi was run with a development version in October 2022 that is slighty different from the released one. Re-running with the released version might lead to slightly different results, which however do not affect downstream analyses and the conclusions of the paper

# in a terminal:
# cd 005-Mochi
# ./run_mochi.py

