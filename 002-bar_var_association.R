library(stringdist)
library(stringr)
library(parallel)
library(rgl)
library(ShortRead)
options(stringsAsFactors = FALSE) 


wtJunAA = "QERIKAERKRMRNRIAASKCRKRKLERIARLEEKVKTLKAQNSELASTANMLREQVAQLKQKVMN*"
wtJunNt = "CAGGAGCGGATCAAGGCGGAGAGGAAGCGCATGAGGAACCGCATCGCTGCCTCCAAGTGCCGAAAAAGGAAGCTGGAGAGAATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTTAAACAGAAAGTCATGAACTAA"

# codon adaptation index table, used here only for the genetic code
cai = read.delim("../../data/CAI_yeast.txt",header=F)
row.names(cai) = cai$V3



################################################################################
# Functions
################################################################################

translateDNA = function(x){
  
  paste(cai[do.call("paste0",as.data.frame(t(matrix(strsplit(x,"")[[1]],nrow=3)))),"V1"],collapse="")
  
}



################################################################################
# Pre-processing
################################################################################

# read merged paired end reads
jun = readFastq("001-merged_reads/merged_jun.assembled.fastq")

# filter reads with average quality < 30
scores_jun = alphabetScore(jun)
length_jun = width(jun)
avQ_jun = scores_jun / length_jun
jun = as.vector(jun@sread)
jun = jun[avQ_jun >= 30]

# count identical reads
jun = as.data.frame(table(jun))
jun = jun[order(-jun$Freq),]
jun$jun = as.character(jun$jun)


# split the read into barcode and variant
# for the barcode, we simply take the 24 first bp. Even if there is an indel and we sequence through the spacer, that is what the barcode will look like when scoring the deepPCA
# for the 5' end of the variant, we search for a perfect match of the sequence left after recloning the intermediate plasmid.
# because if the match is not perfect due to a true mutation, the enzymes wouldn't have cut and they would have dropped from the library. Thus, any difference is for sure a sequencing error. We can discard them since the true sequence is anyway more abundant than the error
# we don't pay attention to errors or true mutations in the spacer. These have no impact at all

# look for a perfect match
chop1_jun = str_split(jun$jun,"TAGGGCTAGC")

# count the number of occurences
n = unlist(lapply(chop1_jun,length))
table(n)

# we keep only sequences were there was exactly one match
tmp = do.call("rbind", chop1_jun[n==2])
jun2 = data.frame(bar = substr(tmp[,1], 1, 24),
                  var = tmp[,2],
                  nReads = jun$Freq[n==2]
)

# the vast majority of unique sequences is 198bp long as expected
table(nchar(jun2$var))
sum(jun2$nReads) / sum(jun$Freq)
# 97.5% of the reads are kept


# aggregate counts for identical barcode-variant pairs
bv = aggregate(jun2$nReads, list(jun2$bar, jun2$var), sum)
names(bv) = c("bc","va", "Npair")


# aggregate counts for unique barcode and variant sequences
unique_bar = aggregate(bv$Npair, list(bv$bc), sum)
unique_var = aggregate(bv$Npair, list(bv$va), sum)
names(unique_bar) = c("bc", "Nbc")
names(unique_var) = c("va", "Nva")

# add information about total barcode and variant counts for each pair
bv$Nbc = unique_bar$Nbc[match(bv$bc, unique_bar$bc)]
bv$Nva = unique_var$Nva[match(bv$va, unique_var$va)]

# compute the ratio of counts for the pair to the total counts for the corresponding barcodes and variants
bv$Rbc = bv$Npair / bv$Nbc
bv$Rva = bv$Npair / bv$Nva

# re-order by decreasing counts
unique_bar = unique_bar[order(-unique_bar$Nbc),]
unique_var = unique_var[order(-unique_var$Nva),]
bv = bv[order(-bv$Npair),]


# save
save(unique_bar, file="000-data/unique_bar_jun.Rdata")
save(unique_var, file="000-data/unique_var_jun.Rdata")
save(bv, file="000-data/bv_jun.Rdata")



################################################################################
# Identification of true barcode-variant pairs (from sequencing errors)
################################################################################

# The overall strategy is described in the methods section of the manuscript

### Step1, identify unique barcode sequences

# most abundant variant associated to each unique barcode
topVar = bv[!duplicated(bv$bc),]

# filter out barcode-variant pairs sequences with less than 10 counts
hist(log10(topVar$Nbc))
hist(log10(topVar$Nbc[topVar$Nbc > 3]), breaks=40)
abline(v=1)
# 10 reads is even still a bit conservative

bv = bv[bv$Npair >= 10, ]

# recalculate total unique barcode and variant reads
unique_bar = aggregate(bv$Npair, list(bv$bc), sum)
unique_var = aggregate(bv$Npair, list(bv$va), sum)
bv$Nbc = unique_bar$x[match(bv$bc, unique_bar$Group.1)]
bv$Nva = unique_var$x[match(bv$va, unique_var$Group.1)]

# recalculate ratio of pair reads to unqiue barcode and variant reads
bv$Rbc = bv$Npair / bv$Nbc
bv$Rva = bv$Npair / bv$Nva

# most abundant variant associated to each unique barcode (many unique variants were filtered out because of less than 10 counts)
topVar = bv[!duplicated(bv$bc),]


# unique barcodes associated to each unique variant (i.e. filter out barcode sequencing errors)
# similar barcodes associated to the exact same variant are likely sequencing errors

va = split(topVar, topVar$va)

# here we put aside the variants associated to a single barcode, since there's no need to collapse
n = unlist(lapply(va,nrow))
u = do.call("rbind",va[n==1])
va = va[n>1]

# then, for every variant
true_bar = do.call("rbind",mclapply(1:length(va),function(x){
  
  # get the barcodes 
  cl = va[[x]]
  
  # measure the distance between the top one and all others
  h = stringdist(cl$bc[1], cl$bc, method="lv", nthread = 1) 
  
  # remove all that are similar, except the top one that we keep (rows are order by decreasing read counts)
  rel = cl[h <= 4,]
  out = rel[1,]
  
  # repeat this by iterating through the remaining ones
  others = cl[h > 4,]
  i = 1
  while(nrow(others)){
    
    h = stringdist(others$bc[1], others$bc, method="lv", nthread = 1)
    rel = others[h <= 4,]
    others = others[h > 4,]
    
    out = rbind(out, rel[1,])
    
    i = i+1
  }
  
  out$Ava = i # number of barcodes associated with that variant
  
  out
  
},mc.cores=24))

u$Ava = 1
true_bar = rbind(true_bar,u)

# this is now a data.frame with only unique barcodes
true_bar = true_bar[order(-true_bar$Nbc),]
row.names(true_bar) = 1:nrow(true_bar)


# We got from ~500k unique barcode sequences, to 30416 that have at least 10 reads, to 28969 by removing barcodes that were clearly errors of true barcodes because they were associated with the same variant and had a similar sequence

# Now, we also want to check whether some barcodes are similar to one another, and are associated to a variant that is also similar (but not exactly the same otherwise it would have been picked up earlier). 


### Step2, flag similar true barcodes
# true barcodes that are too similar could be confounded during deepPCA due to sequencing error
# Let's first check the distribution of distances

# we could do a stringdist matrix, but it takes too much memory
# so we're going sequencially
h = do.call("rbind",mclapply(1:nrow(true_bar), function(x){
  out = as.data.frame(table(stringdist(true_bar$bc[x], true_bar$bc, method="lv", nthread = 1))) 
  out$bc = x # numeric barcode ID, corresponds to row number. This will allow to then identify which barcode has at least one other close barcode
  out
},mc.cores=24))
h$Var1 = as.numeric(as.character(h$Var1))

# then we collapse to obtain the number of barcode pairs for each distance
tmp = aggregate(h$Freq,list(h$Var1),sum)
# we divide by two because each barcode pair is counted twice, A-B and B-A (symmetrical matrix)
# but not the first row because this corresponds to the diagonal of the matrix (i.e. distance of 0, barcode compared to itself)
tmp$x[2:nrow(tmp)] = tmp$x[2:nrow(tmp)] / 2
plot(tmp$Group.1, tmp$x, xlab="Levenstein distance", ylab="No. of barcode pairs", log="y")

# There is one barcode pair with a distance of 1, but it involves the backbone barcode (place keeper sequence in the backbone during cloning. Probably comes from undigested plasmids), so these will be filtered out anyway. Those with 3 or more differences follow the expected trend.
# The trend shows that we indeed didn't expect any barcode pair with 1 or 2 differences. This also shows that we do not expect to sample twice the exact same barcode.
# It supports that we indeed identified all true barcodes and those that remain with distances of 3 or 4 are actually expected by chance
# These should be flagged as having close neighbours We could still keep reads in deepPCA that have 2 mutations or less but don't get closer to the neighbour, but that would biase counts for these barcodes because we would filter out the fraction of their reads that get closer to their neighbour. 
# These are just 84 barcode pairs, i.e. 167 barcodes (one barcode has two different close neighbours), so we should rather filter them out (just flagged for now)

closeNeighbour = unique(h$bc[h$Var1 >= 1 & h$Var1 <= 4])
length(closeNeighbour)

true_bar$closeNeighbour = FALSE
true_bar$closeNeighbour[closeNeighbour] = TRUE


### Step3, flag promiscuous barcodes
# Next step is to flag promiscuous barcodes, i.e. barcodes that are associated to several variants.
# Variants are close in sequence by design, so two close variants might indeed get the same barcode by chance.
# Although this is not really expected given the trend in distances to happen by chance given the trend in the distribution of barcode distances plotted above
# most of the barcode-variant pairs that share the same barcodes are thus expected to be sequencing errors in the variant or chimeras due to polymerase template switching during sequencing library generation.
# Let's try to identify those and filter them out


# First, create all expected sequences from NNS mutagenesis
# these are all WT codons
exp = strsplit(wtJunNt,split="")[[1]]
exp = as.data.frame(matrix(exp,ncol=3,byrow = T))
exp = do.call("paste0",exp)

# all 32 NNS codons
possCod = do.call("paste0",expand.grid(c("A","T","G","C"),c("A","T","G","C"),c("G","C")))

# then we replace all 32 mutated positions by one of the 32 NNS codon
exp = as.vector(sapply(28:59,function(x){
  sapply(1:32,function(y){
    exp[x] = possCod[y]
    paste(exp,collapse="")
  })
}))

exp = unique(exp) # removes when the WT codon is already an NNS codon, which happens at multiple positions


true_bar$exp = true_bar$va %in% exp

# aa distance
tmp_aa = sapply(true_bar$va, translateDNA)
h = stringdist(tmp_aa, wtJunAA, method="hamming")
true_bar$nbMutAA = h




plot3d(log10(bv$Nbc), log10(bv$Nva), log10(bv$Npair), xlab="bc", ylab="va",zlab="pair", size=1, alpha=0.5, col=(bv$Nva > 3000 & bv$Nbc > 300 & bv$Npair < 100) + 1)
plot3d(log10(true_bar$Nbc), log10(true_bar$Nva), log10(true_bar$Npair), xlab="bc", ylab="va",zlab="pair", size=1, alpha=0.5)

# Pairs with low counts but high total counts of the corresponding unique barcode variant are chimeras due to the higher probability to form chimeras for barcodes and variants present at high frequency
# Let's filter them out
bv2 = bv[!(bv$Nva > 3000 & bv$Nbc > 300 & bv$Npair < 100) & bv$bc %in% true_bar$bc,]

# Pairs with a high total read count for the barcode and a low total read count for the variant are likely to be sequencing error on the variant side since it is not expected that the same barcode gets sampled twice
# We filter out pairs where a lower frequency variants shares the same barcode with the highest frequency variant of this barcode if:
# - they differ by only one mutation from the dominant variant
# - this mutation is in another codon and not a reversion towards WT
# - these variants have a total abundance much lower than the dominant one, as well as for the abundance paired with this barcode
# - if the dominant variant is WT, then we also require that the errors do not correspond to expected variants to sum their Rbc

# vector of JUN wt nucleotides
swt = strsplit(wtJunNt,"")[[1]]

bv3 = do.call("rbind", mclapply(1:nrow(true_bar),function(x){
  
  # for each unique barcode
  fbc = true_bar$bc[x]
  # take all variant sequences that share the same barcode
  same_bc = bv2[bv2$bc == fbc,]
  
  # distance between variant and WT
  h1 = stringdist(wtJunNt, true_bar$va[x],method="lv", nthread = 1)
  
  # distance towards dominant variant
  h2 = stringdist(true_bar$va[x],same_bc$va,method="lv", nthread = 1)
  
  # distance towards WT
  h3 = stringdist(wtJunNt,same_bc$va,method="lv", nthread = 1)
  
  # now we need to make sure that the variants that differs by a single base away from WT
  k = h3 == h1 + h2 & h2 == 1
  close_va = same_bc[k,] 
  
  if(nrow(close_va)){
    
    # that this mutation is in a different codon
    s_real = strsplit(true_bar$va[x],"")[[1]]
    s_real = unique(floor((which(s_real != swt) - 1)/ 3))
    cod = do.call("rbind",lapply(strsplit(close_va$va,""),function(y){
      out1 = unique(floor((which(y != swt) - 1)/ 3))
      out1 = sum(!(out1 %in% s_real)) == 1
      
      # there is a run of A that often leads to a sequencing error where on A is missing. Because of this frameshift, the amino acid sequence would be different
      # we thus want to check if by correcting this error, the sequence is now similar
      y = paste(c(y[1:64],"A",y[65:length(y)]), collapse="")
      out2 = y == true_bar$va[x]
      c(out1, out2)
    }))
    cod = cod[,1] | cod[,2]
    
    # that these variants have a much lower abundance
    ab = true_bar$Nva[x] / close_va$Nva > 100 &  true_bar$Npair[x] / close_va$Npair > 10
    
    # that these are expected or not (useful only for WT variants since the others require a second codon to be mutated and can therefore not correspond to expected variants)
    non_exp_var = !(close_va$va %in% exp)
    
    diff_va = close_va[!(cod & ab & non_exp_var),]
    rbind(same_bc[!k,], diff_va)
    
  }else{
    same_bc[!k,]
  }
  
},mc.cores=24))


# many are still left that are just above the threshold but do look like sequencing errors. However, it is better to be stringent here even if we lose some barcodes


# recalculate total unique barcode and variant reads, and ratios
unique_bar = aggregate(bv3$Npair, list(bv3$bc), sum)
unique_var = aggregate(bv3$Npair, list(bv3$va), sum)
bv3$Nbc = unique_bar$x[match(bv3$bc, unique_bar$Group.1)]
bv3$Nva = unique_var$x[match(bv3$va, unique_var$Group.1)]

# recalculate ratio of pair reads to unqiue barcode and variant reads
bv3$Rbc = bv3$Npair / bv3$Nbc
bv3$Rva = bv3$Npair / bv3$Nva

bv3$exp = bv3$va %in% exp

plot(bv3$Rbc, bv3$Nbc, xlab="Rbc", ylab="Nbc", log="y", col=bv3$exp+1, pch=".")
hist(bv3$Rbc, xlab="Rbc", main="all barcodes",  breaks=100)
hist(bv3$Rbc[bv3$exp], xlab="Rbc", main="only expected variants",  breaks=100)

nrow(bv3)
nrow(bv3[bv3$Rbc < 1,])
length(unique(bv3$bc))
length(unique(bv3$bc[bv3$Rbc < 1]))

# we can thus consider as promiscuous barcodes for which Rbc < 1, i.e. the read count for the pair is lower that the total read count for the barcode


# let's now keep only non-promiscuous barcodes, i.e. barcodes non-ambiguously associated with a single variant
true_barcodes = bv3[bv3$Rbc == 1,]
tmp = true_bar$closeNeighbour[match(true_barcodes$bc, true_bar$bc)]
true_barcodes = true_barcodes[!tmp,]


### Step 4, annotate sequence
wtJunAA2 = strsplit(wtJunAA,"")[[1]]

aa = data.frame(va = unique(true_barcodes$va), va_aa = sapply(unique(true_barcodes$va),translateDNA))
# variants with indels will anyway be ignored later

aa2 = do.call("rbind",lapply(1:nrow(aa),function(i){
  
  x = strsplit(aa$va_aa[i],"")[[1]]
  w = which(x != wtJunAA2)
  n = length(w)
  if(n == 1){
    pos = w
    wt_aa = wtJunAA2[w]
    mut_aa = x[w]
    wt_codon = substr(wtJunNt, pos*3-2, pos*3)
    mut_codon = substr(aa$va[i], pos*3-2, pos*3)
    if(w >=28 & w <= 59){
      hept = floor((w-27-1) / 7) + 1
      hept_pos = letters[1:7][(w-27-1) %% 7 + 1]
      id = paste0(hept,hept_pos,mut_aa)
    }else{
      hept = 0
      hept_pos = "z"
      id = paste0(wt_aa,pos,mut_aa)
    }
  }else{
    pos = 0
    wt_aa = paste(wtJunAA2[w],collapse=",")
    mut_aa = paste(x[w],collapse=",")
    wt_codon = ""
    mut_codon = ""
    hept = 0
    hept_pos = "z"
    id = paste(paste0(wtJunAA2[w],w,x[w]),collapse=",")
  }
  
  data.frame(nMut = n, id, pos, hept, hept_pos, wt_aa, mut_aa, wt_codon, mut_codon)
  
}))

aa2 = cbind(aa,aa2)

true_barcodes = merge(true_barcodes,aa2,by="va")
true_barcodes = true_barcodes[order(-true_barcodes$Npair),]
true_barcodes$id[true_barcodes$id == ""] = "0xx"

save(true_barcodes,file="000-data/true_barcodes.Rdata")
write.table(true_barcodes, file="000-data/TableS2_bar_var_assoc.txt",sep="\t",quote = F,row.names=F)
