################################################################################
# Packages
################################################################################

library(readxl)
library(flowCore)
library(parallel)
library(gplots)
library(ggplot2)
library(stringr)


################################################################################
# Functions
################################################################################

panel.cor.pearson <- function(x, y, digits = 2, prefix = "", cex.cor,  ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="pairwise.complete.obs", method="pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}
panel.hist40 <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks=40)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}


################################################################################
# Load data
################################################################################

load("000-data/true_barcodes.Rdata")
load("000-data/binding_scores.Rdata")
# barcode read counts, to compare bidnign variance (Fig S2)
load("000-data/bc_counts.Rdata")
load("000-data/bc_counts_annotation.Rdata")

dG = read.delim("005-Mochi_output/mochi_model_pymochi_modabs/task_1/weights/weights_Binding.txt",header=T)

# data from Diss and lehner 2018
elife = read_xlsx("000-data/elife-32472-supp1-v2.xlsx")


################################################################################
# Barcodes and variants coverage
################################################################################

# number of JUN variants obtained
length(unique(true_barcodes$id[true_barcodes$exp]))

# number of barcodes per unique amino acid variant
n = table(true_barcodes$id[true_barcodes$exp])
median(n)

# number of binding scores obtained
nrow(binding_scores[!binding_scores$filt,])

# number of unique variants
length(unique(binding_scores$id[!binding_scores$filt]))

# number of unique partners
length(unique(binding_scores$bZ[!binding_scores$filt]))

# average number of input reads
mean(c(binding_scores$count_e1_s0[!binding_scores$filt],binding_scores$count_e2_s0[!binding_scores$filt],binding_scores$count_e3_s0[!binding_scores$filt],binding_scores$count_e4_s0[!binding_scores$filt],binding_scores$count_e5_s0[!binding_scores$filt],binding_scores$count_e6_s0[!binding_scores$filt]))


################################################################################
# Fig S1: Heatmap of_absolute binding scores
################################################################################

allIDs = expand.grid(1:5,letters[1:7],unique(substr(binding_scores$id,3,3)))
allIDs = allIDs[allIDs[,3] != "x",]
allIDs = sort(paste0(allIDs[,1],allIDs[,2],allIDs[,3]))
allIDs = c("0xx",allIDs[!(substr(allIDs,1,2) %in% c("5e","5f","5g"))])

all_pairs = as.data.frame(expand.grid(bZ = unique(binding_scores$bZ), id = allIDs))
all_pairs = merge(all_pairs,binding_scores,all.x=T)

# matrix for the heatmap
m = matrix(all_pairs$fitness, ncol=52)
colnames(m) = unique(as.character(all_pairs$bZ))
m = m[,order(all_pairs$fitness[all_pairs$id == "0xx"])]
rownames(m) = allIDs

# color palette
ra = range(m,na.rm=T)
ra[1] = floor(ra[1]*100)/100
ra[2] = ceiling(ra[2]*100)/100
co = colorRampPalette(c("blue","red"))((ra[2]-ra[1])*100)

pdf("008-figures/FigS1_Heatmap_of_absolute_binding_scores.pdf",width=14)
image(1:nrow(m),1:ncol(m),m,xlab="",ylab="",axes=F, col=co, zlim=ra)
box()
axis(1, at = seq(1.5,nrow(m)+0.5,20), labels = F)
axis(1, at = seq(11,nrow(m)-10,20), labels = unique(substr(allIDs[2:length(allIDs)],1,2)), tick=F,cex.axis=0.7,line=-1)
axis(2, at = 1:52, labels=colnames(m), las=2,cex.axis=0.5)

m2 = matrix(seq(ra[1],ra[2],0.01),ncol=1)
image(seq(ra[1],ra[2],0.01),1,m2,col=co,zlim=ra,axes=F,xlab="binding score",ylab="",main="color key")
box()
axis(1)
dev.off()


################################################################################
# Fig 1d: correlation between replicate binding scores
################################################################################

corel = cor(binding_scores[!binding_scores$filt,grep("^fitness._uncorr$",names(binding_scores))])
mean(corel[upper.tri(corel)])

pdf("008-figures/Fig1d_correlation_binding_scores_rep1_and_2.pdf")
plot(binding_scores$fitness1_uncorr[!binding_scores$filt], binding_scores$fitness2_uncorr[!binding_scores$filt],pch=".",xlab="binding scores replicate 1",ylab="binding scores replicate 2")
dev.off()


################################################################################
# Fig S2: Binding scores variance as a function of average input read counts
################################################################################

# compare the binding score variance between individual barcodes and barcode reads collapsed to variants
# instead of binding scores, we compare LFC binding frequency in input and output samples

# we first need to compute these for individual barcodes
# output samples wer sequenced twice. we need to add them up
bc_counts[,7:12] = bc_counts[,7:12] + bc_counts[,13:18]
bc_counts = bc_counts[,1:12]
# filter out counts from non-expected variant-wt partner pairs
tmp = c(F,true_barcodes$exp) # non-recognized barcodes were given an index of 0
bc_counts = bc_counts[exp_jun$wt_id %in% unique(binding_scores$bZ) & tmp[exp_jun$var_id + 1],] # +1 because we just added F as index 1 for the non-recognized barcodes
exp_jun = exp_jun[exp_jun$wt_id %in% unique(binding_scores$bZ) & tmp[exp_jun$var_id + 1],]

# compute frequencies and LFCs
freq = t(t(bc_counts) / colSums(bc_counts))
lfc_bar = log10(bc_counts[,7:12] / bc_counts[,1:6])


# let's compute LFCs the same way for aggregated counts
agg_counts = binding_scores[,grep("^count",names(binding_scores))]
freq2 = t(t(agg_counts) / colSums(agg_counts))
lfc_agg = log10(agg_counts[,7:12] / agg_counts[,1:6])
var_agg = apply(lfc_agg,1,var)


# compute within-replicate variance across all barcodes mapped to the same variant-wt partner pair
pair_list = split(data.frame(bc_counts[,1:6],lfc_bar), list(exp_jun$wt_id, true_barcodes$id[exp_jun$var_id]), drop=T)
within_rep_variance = do.call("cbind",lapply(pair_list,function(x){
  
  sapply(1:6,function(i){
    
    k = is.finite(x[,i+6])
    mean_count = mean(x[k,i])
    var_lfc = var(x[k,i+6])
    c(mean_count,var_lfc)
    
  })
  
}))
within_rep_variance = within_rep_variance[,is.finite(colSums(within_rep_variance))]

# compute across-replicate variance for LFC computed from aggregated counts
var_agg = apply(lfc_agg,1,var)
mean_counts_agg = rowMeans(agg_counts)
k = is.finite(mean_counts_agg) & is.finite(var_agg)
var_agg = var_agg[k]
mean_counts_agg = mean_counts_agg[k]


# plot contours
# first, we bin the mean counts and variance for both according to the same breaks
x = seq(log10(min(c(within_rep_variance[1,],mean_counts_agg))), log10(max(c(within_rep_variance[1,],mean_counts_agg))), length.out = 20)
y = seq(log10(min(c(within_rep_variance[,2],var_agg))), log10(max(c(within_rep_variance[,2],var_agg))), length.out = 20)

bar = table(data.frame(cut(log10(within_rep_variance[1,]),breaks=x),cut(log10(within_rep_variance[2,]),breaks=y)))
agg = table(data.frame(cut(log10(mean_counts_agg),breaks=x),cut(log10(var_agg),breaks=y)))

# bin centers
xdim = (x[-1] + x[-length(x)]) / 2
ydim = (y[-1] + y[-length(y)]) / 2

pdf("008-figures/FigS2_variance_barcodes.pdf")
contour(xdim,ydim,log2(bar),xlab="average input count (log10)",ylab="LFC variance (log10)")
contour(xdim,ydim,log2(agg),add=T,col=2)
legend("topright",fill=1:2,legend=c("barcodes","aggregated"))
dev.off()


################################################################################
# Fig S3: Correlation with binding scores from an earlier study 
################################################################################

elife_singles = elife[!duplicated(elife$id2),]

binding_scores_fos = binding_scores[binding_scores$bZ == "FOS",]
binding_scores_fos$elife = elife_singles$s2_mppi[match(binding_scores_fos$id, elife_singles$id2)]

r = cor.test(binding_scores_fos$elife, binding_scores_fos$fitness, use="pairwise.complete.obs")

pdf("008-figures/FigS3_comparison_with_elife.pdf")
plot(binding_scores_fos$elife, binding_scores_fos$fitness, xlab="binding score, Diss & Lehner 2018", ylab="binding scores, this study")
text(0.5, 0, paste0("R = ", sprintf("%.2f",r$estimate)), cex=2)
dev.off()


################################################################################
# Fig 1e: Distribution binding scores
################################################################################

pdf("008-figures/Fig1e_distribution_binding_Scores.pdf")
hist(binding_scores$fitness[!binding_scores$filt],breaks=40,xlab="Binding score", ylab="Counts (x10^3)")
dev.off()


################################################################################
# Fig S4: Confirmation by splitFAST
################################################################################

# load fcs files
fcs_files = dir("000-data/FACS_data")
fcs = sapply(file.path("000-data/FACS_data",fcs_files), read.FCS, alter.names=T)
fcs <- as(fcs, "flowSet")

# format experiment names
exp_names = strsplit(fcs_files, "_")
exp_names = do.call("rbind",lapply(exp_names,function(x){
  if(x[1] == "NS"){
    strsplit(x[3],"x")[[1]]
  }else{
    x = strsplit(x[2],"-")[[1]][2]
    strsplit(x,"x")[[1]]
  }
  
}))
# some controls are present to very the gating
# non-transfected cells
# a plasmid expressing BFP only to check green channel signal from blue fluorescence
# a plasmid expressing NFAST and CFAST fragments to check for auto-complementation


pdf("007-figures/gating.pdf",width=14,height=14)
par(mfrow=c(2,2))
av_score = unlist(lapply(seq_along(fcs), function(i){
  
  x = as.data.frame(exprs(fcs[[i]]))
  
  plot(x$FSC.A, x$SSC.A, pch=".", log="xy", main = paste(exp_names[i,],collapse=" x "))
  abline(v=c(42000,260000))
  abline(h=c(25000,250000))
  x = x[x$FSC.A > 42000 & x$FSC.A < 260000 & x$SSC.A > 25000 & x$SSC.A < 250000,]
  
  plot(x$SSC.H, x$SSC.W, pch=".", log="xy")
  abline(v=c(50000,200000))
  abline(h=c(20000,85000))
  x = x[x$SSC.H > 50000 & x$SSC.H < 200000 & x$SSC.W > 20000 & x$SSC.W < 85000 ,]
  
  plot(x$Alexa.Fluor.405.A, x$GFP.A, pch=".", log="xy")
  abline(0,1)
  abline(v=500)
  x = x[!(x$GFP.A > x$Alexa.Fluor.405.A & x$Alexa.Fluor.405.A < 500),]
  
  x$ratio = x$GFP.A / x$Alexa.Fluor.405.A
  
  if(!(exp_names[i,1] %in% c("empty","CFAST"))){
    hist(x$ratio)
  }else{
    plot(1,type="n")
  }
  
  
  mean(x$ratio[is.finite(x$ratio)])
  
}))
dev.off()

exp_names[exp_names[,1] == "WT",1] = "0xx"
av_score = data.frame(id = exp_names[,1], bZ = exp_names[,2], splitFAST = av_score)

av_score2 = aggregate(av_score$splitFAST, list(av_score$id, av_score$bZ), mean)
names(av_score2) = c("id","bZ","splitFAST_score")
av_score2$splitFAST_se = aggregate(av_score$splitFAST, list(av_score$id, av_score$bZ), sd)$x / sqrt(3)

av_score2$deepPCA_score = binding_scores$fitness[match(paste(av_score2$id, av_score2$bZ), paste(binding_scores$id, binding_scores$bZ))]
av_score2$deepPCA_se = binding_scores$sigma[match(paste(av_score2$id, av_score2$bZ), paste(binding_scores$id, binding_scores$bZ))]

av_score2 = av_score2[!(av_score2$id %in% c("CFAST","empty")),]


# calculate relative binding scores
av_score3 = av_score2[av_score2$id != "0xx",]
wt = av_score2[av_score2$id == "0xx",]
av_score3$splitFAST_wt = wt$splitFAST_score[match(av_score3$bZ, wt$bZ)]
av_score3$splitFAST_se_wt = wt$splitFAST_se[match(av_score3$bZ, wt$bZ)]
av_score3$deepPCA_wt = wt$deepPCA_score[match(av_score3$bZ, wt$bZ)]
av_score3$deepPCA_se_wt = wt$deepPCA_se[match(av_score3$bZ, wt$bZ)]

av_score3$splitFAST_diff = av_score3$splitFAST_score - av_score3$splitFAST_wt
av_score3$splitFAST_diff_se = sqrt(av_score3$splitFAST_se_wt^2 + av_score3$splitFAST_se^2)
av_score3$deepPCA_diff = av_score3$deepPCA_score - av_score3$deepPCA_wt
av_score3$deepPCA_diff_se = sqrt(av_score3$deepPCA_se_wt^2 + av_score3$deepPCA_se^2)


pdf("007-figures/FigS4_confirmation_splitFAST.pdf",width=14)
par(mfrow=c(1,2))

conf.int = cbind(av_score2$deepPCA_score - av_score2$deepPCA_se*1.96, av_score2$deepPCA_score + av_score2$deepPCA_se*1.96, av_score2$splitFAST_score - av_score2$splitFAST_se*1.96, av_score2$splitFAST_score + av_score2$splitFAST_se*1.96)
r = cor.test(av_score2$deepPCA_score, av_score2$splitFAST_score)
plot(av_score2$deepPCA_score, av_score2$splitFAST_score, xlab="deepPCA binding score", ylab="splitFAST binding score", xlim=range(conf.int[,1:2]), ylim=range(conf.int[,3:4]), main = "Absolute binding scores")
text(-1.3,0.6,paste("R =", sprintf("%.3f",r$estimate), "\np =", sprintf("%.1e",r$p.value)))
segments(conf.int[,1], av_score2$splitFAST_score, conf.int[,2], av_score2$splitFAST_score)
segments(av_score2$deepPCA_score, conf.int[,3], av_score2$deepPCA_score, conf.int[,4])


conf.int = cbind(av_score3$deepPCA_diff - av_score3$deepPCA_diff_se*1.96, av_score3$deepPCA_diff + av_score3$deepPCA_diff_se*1.96, av_score3$splitFAST_diff - av_score3$splitFAST_diff_se*1.96, av_score3$splitFAST_diff + av_score3$splitFAST_diff_se*1.96)
r = cor.test(av_score3$deepPCA_diff, av_score3$splitFAST_diff)
plot(av_score3$deepPCA_diff, av_score3$splitFAST_diff, xlab="deepPCA binding score, relative to wild-type", ylab="splitFAST binding score, relative to wild-type", xlim=range(conf.int[,1:2]), ylim=range(conf.int[,3:4]), main = "Relative binding scores")
text(-1.05,0.21,paste("R =", sprintf("%.3f",r$estimate), "\np =", sprintf("%.1e",r$p.value)))
segments(conf.int[,1], av_score3$splitFAST_diff, conf.int[,2], av_score3$splitFAST_diff)
segments(av_score3$deepPCA_diff, conf.int[,3], av_score3$deepPCA_diff, conf.int[,4])
dev.off()

write.table(av_score3, file="000-data/TableS4_splitFAST_scores.txt", row.names=F,sep="\t")



################################################################################
# Fig 1f: barplot binding scores wt partners
################################################################################

wt_binding_scores = binding_scores$fitness[binding_scores$id == "0xx"]
names(wt_binding_scores) = binding_scores$bZ[binding_scores$id == "0xx"]
wt_ci = binding_scores$sigma[binding_scores$id == "0xx"]

pdf("007-figures/Fig1f_barplot_wt_partners.pdf",height=14)
barplot2(sort(wt_binding_scores)-min(wt_binding_scores)+0.05,las=2,cex.names=0.7, offset = min(wt_binding_scores)-0.05, plot.ci = T, ci.l = sort(wt_binding_scores) - wt_ci[order(wt_binding_scores)], ci.u = sort(wt_binding_scores) + wt_ci[order(wt_binding_scores)], xlim = range(wt_binding_scores), horiz=T, axes=F)
axis(3)
mtext("Binding score with wild-type JUN", line=2)
dev.off()


################################################################################
# Fig 2a: heatmap of mutational effects
# Fig S5: corresponding p-values
################################################################################

# m is the matrix of absolute binding cores plotted in Fig S1
m_rel = t(t(m) - sort(wt_binding_scores))
# remove wt
m_rel = m_rel[2:nrow(m_rel),]

ra = range(m_rel,na.rm=T)
ra[1] = floor(ra[1]*100)/100
ra[2] = ceiling(ra[2]*100)/100

co = c(colorRampPalette(c("blue","grey90"))(abs(ra[1]*100)), colorRampPalette(c("grey90","red"))(abs(ra[2]*100)))

pdf("007-figures/Fig2a_heatmap_mutational_effects.pdf",width=14)
image(1:nrow(m_rel),1:ncol(m_rel),m_rel,xlab="",ylab="",axes=F, col=co, zlim=ra)
box()
axis(1, at = seq(0.5,nrow(m_rel)+0.5,20), labels = F)
axis(1, at = seq(10,nrow(m_rel)-10,20), labels = unique(substr(allIDs[2:length(allIDs)],1,2)), tick=F,cex.axis=0.7,line=-1)
axis(2, at = 1:52, labels=colnames(m_rel), las=2,cex.axis=0.5)

m2 = matrix(seq(ra[1],ra[2],0.01),ncol=1)
image(seq(ra[1],ra[2],0.01),1,m2,col=co,zlim=ra,axes=F,xlab="binding score",ylab="",main="color key")
box()
axis(1)
dev.off()




mp = matrix(-log10(all_pairs$rel_pval), ncol=52)
colnames(mp) = unique(as.character(all_pairs$bZ))
mp = mp[,order(all_pairs$fitness[all_pairs$id == "0xx"])]
rownames(mp) = allIDs
mp = mp[2:nrow(mp),]

ra = range(mp,na.rm=T)
ra[1] = floor(ra[1]*100)/100
ra[2] = ceiling(ra[2]*100)/100

f = max(all_pairs$rel_pval[all_pairs$rel_fdr < 0.05], na.rm=T)
co = heat.colors((ra[2] - ra[1] ) * 100)
co[1:(-log10(f)*100)] = "#000000"


pdf("007-figures/FigS5_heatmap_mutational_effects_pval.pdf",width=14)
image(1:nrow(mp),1:ncol(mp),mp,xlab="",ylab="",axes=F, col=co, zlim=ra)
box()
axis(1, at = seq(0.5,nrow(mp)+0.5,20), labels = F)
axis(1, at = seq(10,nrow(mp)-10,20), labels = unique(substr(allIDs[2:length(allIDs)],1,2)), tick=F,cex.axis=0.7,line=-1)
axis(2, at = 1:52, labels=colnames(mp), las=2,cex.axis=0.5)

m2 = matrix(seq(ra[1],ra[2],0.01),ncol=1)
image(seq(ra[1],ra[2],0.01),1,m2,col=co,zlim=ra,axes=F,xlab="-log10(pval)",ylab="",main="color key")
box()
axis(1)
dev.off()


################################################################################
# Fig 2b: comparison binding scores across partners
# Fig S6: all comparisons
################################################################################

tmp = matrix(all_pairs$rel, ncol=52)
tmp[all_pairs$filt] = NA
colnames(tmp) = unique(as.character(all_pairs$bZ))
tmp = tmp[,order(all_pairs$fitness[all_pairs$id == "0xx"])]
rownames(tmp) = allIDs

pdf("007-figures/Fig2b_comparison_binding_scores_across_partners.pdf")
pairs(tmp[,c("JDP2","FOS","CREB5","ATF7")],upper.panel = NULL)
dev.off()

pdf("007-figures/FigS6_comparison_binding_scores_across_all_partners.pdf",width=28,height=28)
pairs(tmp[,ncol(tmp):1], pch=".", upper.panel = panel.cor.pearson, diag.panel = panel.hist40)
dev.off()


################################################################################
# Fig 3b: binding scores prediction vs measurements
################################################################################

pdf("007-figures/Fig3b_predicted_vs_observed_binding_scores.pdf")
plot(binding_scores$mean[!binding_scores$filt], binding_scores$fitness[!binding_scores$filt], pch=".", xlab="Predicted binding score", ylab="Observed binding score")
r = cor(binding_scores$fitness[!binding_scores$filt], binding_scores$mean[!binding_scores$filt], use="pairwise.complete.obs")^2
text(-1.5,0.5,paste0("R^2 =", sprintf("%.2f",r)), pos = 4)
dev.off()



################################################################################
# Fig 3c: Comparison with in vitro Kd
################################################################################

iv_kd = read.csv("000-data/Supplementary_Homosapiens_Kd_37degrees.csv",skip=1)
rownames(iv_kd) = iv_kd[,1]
iv_kd = iv_kd[,2:ncol(iv_kd)]
iv_kd2 = apply(iv_kd,2,as.numeric)
iv_kd3 = as.matrix(iv_kd)
iv_kd2[iv_kd3==">5000"] = 5000
iv_kd2[iv_kd3=="<1"] = 1
rownames(iv_kd2) = rownames(iv_kd)

partners = binding_scores[!duplicated(binding_scores$bZ),]
partners = partners[partners$bZ %in% colnames(iv_kd2),]
partners$Keit1 = log(iv_kd2["JUN",partners$bZ]*10^(-9)) * 1.987 * 310.15 / 1000 # divide by 1000 to have units in kcal/mol 
partners$Keit2 = log(iv_kd2[partners$bZ,"JUN"]*10^(-9)) * 1.987 * 310.15 / 1000
partners$Keit1 = partners$Keit1 - partners$Keit1[partners$bZ == "FOS"]
partners$Keit2 = partners$Keit2 - partners$Keit2[partners$bZ == "FOS"]

pdf("007-figures/Fig3c_comparison_invitro_Kds.pdf")
pairs(partners[,c("ddG_bZ","Keit1","Keit2")],upper.panel = panel.cor.pearson)
dev.off()



################################################################################
# Fig 4: Prediction global effects
################################################################################


aai = read.delim("000-data/amino_acid_properties_annotated_supplementary.txt")
aai_id = aai[,2]
names(aai_id) = aai[,1]
aai = as.matrix(aai[,3:22])
rownames(aai) = names(aai_id)
# remove incomplete features
aai_id = aai_id[!is.na(rowSums(aai))]
aai = aai[!is.na(rowSums(aai)),]

dG2 = dG[dG$Pos >= 28 & dG$Pos <= 59 & !is.na(dG$Pos),]
dG2$h = as.character(floor((dG2$Pos - 28)/ 7) + 1)
dG2$hp = letters[(dG2$Pos - 28) %% 7 + 1]
dG2$id2 = paste0(dG2$h,dG2$hp,substr(dG2$id,4,4))

# featutres for the wt residues
pr_wt = as.data.frame(t(aai[,substr(dG2$id,1,1)]))
# features for the mutant residues
pr_mut = as.data.frame(t(aai[,substr(dG2$id,4,4)]))
# change in features upon mutation
pr_diff = pr_mut - pr_wt

pr_diff = data.frame(ddG = dG2$mean , h = dG2$h, hp = dG2$hp, pr_diff)

# identify best feature
r = sapply(2:ncol(pr_diff),function(j){
  
  df2 = data.frame(ddG = pr_diff$ddG,pr_diff[,j])
  summary(lm(ddG ~ .,data=df2))$r.squared
  
})

# store best feature and remove from the feature matrix
df = data.frame(pr_diff[,which.max(r)+1])
names(df) = names(pr_diff)[which.max(r)+1]
pr_diff = pr_diff[,-(which.max(r)+1)]

# now forward search for best features
for(i in 1:4){
  
  cat("\r", i, "\t\t\t")
  
  r = sapply(2:ncol(pr_diff),function(j){
    
    df2 = data.frame(ddG = pr_diff$ddG,df,pr_diff[,j])
    
    summary(lm(ddG ~ .,data=df2))$r.squared
    
  })
  
  df = cbind(df,pr_diff[,which.max(r)+1])
  names(df)[ncol(df)] = names(pr_diff)[which.max(r)+1]
  pr_diff = pr_diff[,-(which.max(r)+1)]
  
}

r = NULL
for(i in 1:ncol(df)){
  
  df2 = data.frame(ddG = pr_diff$ddG, df[,1:i])
  r = c(r,summary(lm(ddG ~ .,data=df2))$r.squared)
  
}
names(r)[!(names(df) %in% c("h","hp"))] = aai_id[names(df)[!(names(df) %in% c("h","hp"))]]
names(r)[c(2,4)] = c("heptad position","heptad")


pdf("007-figures/Fig4a_pred_aaIndexes.pdf",width=14)
par(mar=c(5,30,4,1))
plot(r[1:5],5:1,xlim=c(0,max(r[1:5])), xlab="",ylab="",axes=F)
box()
axis(3)
mtext("Cumulative proportion of variance explained",line=3)
axis(2, at = 5:1, labels=names(r),las=2)
lines(r[1:5],5:1)
dev.off()


l = lm(ddG ~ .,data=df2)
pdf("007-figures/Fig4b_obs_vs_pred.pdf")
plot(predict(l), df2$ddG, type="n", xlab="ddG predicted from amino acid features", ylab="ddG inferred from MoCHI")
text(predict(l), df2$ddG, dG2$id2, col=as.numeric(as.factor(dG2$hp)))
dev.off()



# per hepatd position
pr_diff = pr_mut - pr_wt # reset the features matrix
pr_diff = data.frame(ddG = dG2$mean ,h = dG2$h, hp = dG2$hp, pr_diff)

h = sapply(letters[1:7],function(x){
  
  tmp = pr_diff[pr_diff$hp == x,]
  r = sapply(4:ncol(tmp),function(j){
    
    df2 = data.frame(ddG = tmp$ddG,tmp[,j])
    summary(lm(ddG ~ . ,data=df2))$r.squared
    
  })
  w = which.max(r)
  c(r[w],aai_id[w])
})


pdf("007-figures/Fig4c_amino_acid_property_per_position.pdf",height=14)
par(mar=c(25,4,4,1))
plot(1:7,h[1,7:1],axes=F,xlab="",ylab="Proportion of variance explained")
box()
axis(2)
axis(1,at = 1:7, label = paste(letters[7:1],"-", h[2,]),las=2)
dev.off()



# proportion of variance in binding explained scores per heptad position
# first split by heptad position
tmp = split(binding_scores[!binding_scores$filt,], f=substr(binding_scores$id[!binding_scores$filt],2,2))
tmp = unlist(lapply(tmp,function(x){
  
  cor(x$fitness, x$mean)^2
  
}))

pdf("007-figures/Fig4d_var_explained_by_heptad_position.pdf")
plot(1:7,tmp[letters[1:7]],xlab="Heptad positions",ylab="Proportion of variance explained",axes=F)
box()
axis(2,at=seq(0.8,1,0.05))
axis(1,at=1:7,labels=NA)
text(1:7,tmp[letters[1:7]],letters[1:7],pos=1)
dev.off()



# proportion of variance in binding scores explained per position
# first split by heptad position
tmp = split(binding_scores[!binding_scores$filt,], f=substr(binding_scores$id[!binding_scores$filt],1,2))
tmp = unlist(lapply(tmp,function(x){
  
  cor(x$fitness, x$mean)^2
  
}))
tmp = tmp[2:length(tmp)]

pdf("007-figures/Fig4e_var_explained_by_position.pdf")
plot(1:32,tmp,xlab="JUN's positions",ylab="Proportion of variance explained",axes=F,ylim=c(0.6,1))
box()
axis(2,at=seq(0.6,1,0.1))
axis(1,at=1:32,labels=NA)
text(1:32,tmp,names(tmp),pos=1)
dev.off()



################################################################################
# Fig 5a and S7: Specificity plots
################################################################################

lin = read.delim("005-Mochi_output/mochi_model_pymochi_modabs/task_1/weights/linears_weights_Binding.txt")
lin = colMeans(lin[,2:3])
ref_dG = dG$mean[dG$id == "WT"]

bl = split(binding_scores, binding_scores$bZ)

pdf("007-figures/Fig5a_specificity_plots.pdf",height=28)
par(mar=c(0,0,0,0),oma=c(4,4,4,4),mfrow=c(4,1),xpd=T)
for(i in 1:4){
  
  x = bl[[c("JDP2","FOS","CREB5","ATF7")[i]]]
  
  plot(x$ddG_mut, x$fitness, type="n", axes=F,xlab="",ylab="")
  text(x$ddG_mut, x$fitness, x$id, col=as.numeric(as.factor(substr(x$id,2,2))))
  box()
  if(i < 4){
    axis(1,labels=NA)
  }else{
    axis(1)
  }
  axis(2)
  mtext(c("JDP2","FOS","CREB5","ATF7")[i],side=4)
  
  s = seq(min(binding_scores$ddG_mut), max(binding_scores$ddG_mut), 0.01)
  model_pred = lin[1] / (1+exp(x$dG_mult[1]*(ref_dG+s+x$ddG_bZ[1]))) + lin[2]
  lines(s,model_pred)
  
  
}
mtext("Global binding energy (ddG)",side=1,outer=T,line=2)
mtext("Binding score",side=2,outer=T,line=2)
dev.off()



pdf("007-figures/FigS7_specificity_plots.pdf")
par(mar=c(0,0,0,0),oma=c(4,4,4,4),mfrow=c(4,2),xpd=T)
for(i in 1:length(bl)){
  
  x = bl[[i]]
  
  plot(x$ddG_mut, x$fitness, type="n", axes=F,xlab="",ylab="", ylim=range(binding_scores$fitness))
  text(x$ddG_mut, x$fitness, x$id, col=as.numeric(as.factor(substr(x$id,2,2))))
  box()
  if(i %% 2 == 1){
    axis(2)
  }else{
    axis(2, labels = F)
  }
  if(i %% 7 == 0 | i %% 8 == 0){
    axis(1)
  }else{
    axis(1, labels=F)
  }
  text(5.5,0.4,names(bl)[i])
  
  s = seq(min(binding_scores$ddG_mut), max(binding_scores$ddG_mut), 0.01)
  model_pred = lin[1] / (1+exp(x$dG_mult[1]*(ref_dG+s+x$ddG_bZ[1]))) + lin[2]
  lines(s,model_pred)
  
  if(i %in% seq(1,length(bl),8)){
    mtext("Global binding energy (ddG)",side=1,outer=T,line=2)
    mtext("Binding score",side=2,outer=T,line=2)
  }
  
  
}
dev.off()



################################################################################
# Fig 5c: distribution global and specific effects
################################################################################

pdf("007-figures/Fig5c_distri_glo_and_spe_effects.pdf")
par(mfrow=c(2,1))
hist(binding_scores$glob[!binding_scores$filt],breaks=40,xlab="Global effects", ylab="Frequency (x10^3)", main="", xlim=range(c(binding_scores$glob[!binding_scores$filt],binding_scores$spe[!binding_scores$filt])))
hist(binding_scores$spe[!binding_scores$filt],breaks=40,xlab="Global effects", ylab="Frequency (x10^3)", main="", xlim=range(c(binding_scores$glob[!binding_scores$filt],binding_scores$spe[!binding_scores$filt])))
dev.off()


nrow(binding_scores[binding_scores$glob_fdr < 0.05 & !is.na(binding_scores$glob_fdr),])
nrow(binding_scores[binding_scores$spe_fdr < 0.05 & !is.na(binding_scores$spe_fdr),])


################################################################################
# Fig 5de: compare mutational, global and specific effects
################################################################################

pdf("007-figures/Fig5d_mutational_vs_specific_effects.pdf")
par(mar=c(4,4,4,4))
plot(binding_scores$spe[!binding_scores$filt & binding_scores$rel_fdr >= 0.05 & !is.na(binding_scores$spe_fdr)], # here we use !is.na for spe_fdr instead of rel_fdr because it accounts for the filtering of ATF4
     binding_scores$rel[!binding_scores$filt & binding_scores$rel_fdr >= 0.05 & !is.na(binding_scores$spe_fdr)],
     xlab = "Specific effect", ylab='Mutational effect',
     pch=".",
     xlim=range(binding_scores$spe[!binding_scores$filt & !is.na(binding_scores$spe_fdr)]),
     ylim=range(binding_scores$rel[!binding_scores$filt & !is.na(binding_scores$spe_fdr)]),
)
points(binding_scores$spe[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$spe_fdr) & binding_scores$rel > 0],
       binding_scores$rel[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$spe_fdr) & binding_scores$rel > 0],
       col="red",
       pch=16,
       cex=0.5
)
points(binding_scores$spe[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$spe_fdr) & binding_scores$rel < 0],
       binding_scores$rel[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$spe_fdr) & binding_scores$rel < 0],
       col="cornflowerblue",
       pch=16,
       cex=0.5
)
abline(v=0,lty=3)
abline(h=0,lty=3)
dev.off()


pdf("007-figures/Fig5e_mutational_vs_global_effects.pdf")
par(mar=c(4,4,4,4))
plot(binding_scores$glob[!binding_scores$filt & binding_scores$rel_fdr >= 0.05 & !is.na(binding_scores$glob_fdr)], # here we use !is.na for glob_fdr instead of rel_fdr because it accounts for the filtering of ATF4
     binding_scores$rel[!binding_scores$filt & binding_scores$rel_fdr >= 0.05 & !is.na(binding_scores$glob_fdr)],
     xlab = "Global effect", ylab='Mutational effect',
     pch=".",
     xlim=range(binding_scores$glob[!binding_scores$filt & !is.na(binding_scores$glob_fdr)]),
     ylim=range(binding_scores$rel[!binding_scores$filt & !is.na(binding_scores$glob_fdr)]),
)
points(binding_scores$glob[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$glob_fdr) & binding_scores$rel > 0],
       binding_scores$rel[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$glob_fdr) & binding_scores$rel > 0],
       col="red",
       pch=16,
       cex=0.5
)
points(binding_scores$glob[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$glob_fdr) & binding_scores$rel < 0],
       binding_scores$rel[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$glob_fdr) & binding_scores$rel < 0],
       col="cornflowerblue",
       pch=16,
       cex=0.5
)
abline(v=0,lty=3)
abline(h=0,lty=3)
dev.off()




################################################################################
# Fig 6a, S8,9: heatmap of specific effects, corresponding p-values and dendrogram
################################################################################

m = matrix(all_pairs$spe[all_pairs$bZ != "ATF4" & all_pairs$id != "0xx"], ncol=51)
colnames(m) = unique(as.character(all_pairs$bZ[all_pairs$bZ != "ATF4"]))
rownames(m) = allIDs[allIDs != "0xx"]
mp = matrix(-log10(all_pairs$spe_pval[all_pairs$bZ != "ATF4" & all_pairs$id != "0xx"]), ncol=51)
m[is.na(mp)] = NA
co = cor(m,use="pairwise.complete.obs")
hc = hclust(as.dist(1-co))

m = m[,hc$order]
mp = mp[,hc$order]

ra = range(m,na.rm=T)
ra[1] = floor(ra[1]*100)/100
ra[2] = ceiling(ra[2]*100)/100

# we make the color scale symmetrical by using max(abs(ra))
color = c(colorRampPalette(c("blue","grey90"))(max(abs(ra*100))), colorRampPalette(c("grey90","red"))(max(abs(ra*100))))

pdf("007-figures/Fig6a_heatmap_specificity.pdf",width=14)
image(1:nrow(m),1:ncol(m),m,xlab="",ylab="",axes=F, col=color, zlim=c(-max(ra),max(ra)))
box()
axis(1, at = seq(0.5,nrow(m)+0.5,20), labels = F)
axis(1, at = seq(10,nrow(m)-10,20), labels = unique(substr(allIDs[allIDs != "0xx"],1,2)), tick=F,cex.axis=0.7,line=-1)
axis(2, at = 1:51, labels=colnames(m), las=2,cex.axis=0.5)

m2 = matrix(seq(ra[1],ra[2],0.01),ncol=1)
image(x=1:nrow(m2),y=1,z=m2,col=color,zlim=c(-max(ra),max(ra)),axes=F,xlab="Specific effect",ylab="",main="color key")
box()
axis(1, at = c(46, 96, 146, 196, 246)+0.5, label = c(-0.5,0,0.5,1,1.5))
dev.off()



f = max(all_pairs$spe_pval[!is.na(all_pairs$spe_fdr) & all_pairs$spe_fdr < 0.05], na.rm=T)
ra = range(mp,na.rm=T)
ra[1] = floor(ra[1]*100)/100
ra[2] = ceiling(ra[2]*100)/100

co = heat.colors((ra[2] - ra[1] ) * 100)
co[1:(-log10(f)*100)] = "#000000"

pdf("007-figures/FigS8_heatmap_pval_specificity.pdf",width=14)
image(1:nrow(mp),1:ncol(mp),mp,xlab="",ylab="",axes=F, col=co, zlim=ra)
box()
axis(1, at = seq(0.5,nrow(mp)+0.5,20), labels = F)
axis(1, at = seq(10,nrow(mp)-10,20), labels = unique(substr(allIDs[allIDs != "0xx"],1,2)), tick=F,cex.axis=0.7,line=-1)
axis(2, at = 1:51, labels=colnames(mp), las=2,cex.axis=0.5)

m2 = matrix(seq(ra[1],ra[2],0.01),ncol=1)
image(seq(ra[1],ra[2],0.01),1,m2,col=co,zlim=ra,axes=F,xlab="-log10(pval)",ylab="",main="color key")
box()
axis(1)
dev.off()


pdf("007-figures/Fig6a_dendrogram_partners.pdf")
plot(hc)
plot(as.dendrogram(hc))
dev.off()


m2 = m[rowSums(is.na(m)) < ncol(m) - 10,]
co = cor(t(m2),use="pairwise.complete.obs")
hc2 = hclust(as.dist(1-co))

pdf("007-figures/FigS9_dendro_var.pdf",width=16)
par(cex=0.2)
plot(as.dendrogram(hc2))
dev.off()

################################################################################
# FigS10: Numbers of specificity mutations and positions
################################################################################

# number of pairs with strong changes in specificity
nrow(binding_scores[
    !binding_scores$filt & 
    !is.na(binding_scores$spe_fdr) & 
    binding_scores$spe_fdr < 0.05 & 
    abs(binding_scores$spe) > 0.5,
])


# number of pairs with strong positive changes in specificity
pos_pairs = nrow(binding_scores[
    !binding_scores$filt & 
    !is.na(binding_scores$spe_fdr) & 
    binding_scores$spe_fdr < 0.05 & 
    binding_scores$spe > 0.5,
])

# number of pairs with strong negative changes in specificity
neg_pairs = nrow(binding_scores[
    !binding_scores$filt & 
    !is.na(binding_scores$spe_fdr) & 
    binding_scores$spe_fdr < 0.05 & 
    binding_scores$spe < -0.5,
])

# number of unique JUN variants with strong positive changes in specificity
pos_var = unique(binding_scores$id[
    !binding_scores$filt & 
    !is.na(binding_scores$spe_fdr) & 
    binding_scores$spe_fdr < 0.05 & 
    binding_scores$spe > 0.5
])

# number of unique JUN variants with strong negative changes in specificity
neg_var = unique(binding_scores$id[
    !binding_scores$filt & 
    !is.na(binding_scores$spe_fdr) & 
    binding_scores$spe_fdr < 0.05 & 
    binding_scores$spe < -0.5
])

dual_var = pos_var[pos_var %in% neg_var]
pos_var = pos_var[!(pos_var %in% dual_var)]
neg_var = neg_var[!(neg_var %in% dual_var)]

# number of partners affected by strictly positive specificity mutations
pos_part = unique(binding_scores$bZ[
    binding_scores$id %in% pos_var &
    !binding_scores$filt & 
    !is.na(binding_scores$spe_fdr) & 
    binding_scores$spe_fdr < 0.05 & 
    binding_scores$spe > 0.5
])

# number of partners affected by strictly negative specificity mutations
neg_part = unique(binding_scores$bZ[
    binding_scores$id %in% neg_var &
    !binding_scores$filt & 
    !is.na(binding_scores$spe_fdr) & 
    binding_scores$spe_fdr < 0.05 & 
    binding_scores$spe < -0.5
])

# number of partners affected by dual specificity mutations
dual_part = unique(binding_scores$bZ[
    binding_scores$id %in% dual_var &
    !binding_scores$filt & 
    !is.na(binding_scores$spe_fdr) & 
    binding_scores$spe_fdr < 0.05 & 
    abs(binding_scores$spe) > 0.5
])


pdf("007-figures/FigS10_quantification_specificity_mutations.pdf")
par(mfrow=c(2,2))
barplot(c(neg_pairs,pos_pairs), col = c("cornflowerblue","red"),names.arg = c("decreased","increased"),ylab="Number of variant:partner pairs",xlab="")
mtext("specific effects on binding",side=1,line=3)
barplot(c(length(neg_var),length(dual_var),length(pos_var)), col = c("cornflowerblue","purple","red"),names.arg = c("decreased","dual","increased"),ylab="Number of unique variants",xlab="")
mtext("specific effects on binding",side=1,line=3)
barplot(c(length(neg_part),length(dual_part),length(pos_part)), col = c("cornflowerblue","purple","red"),names.arg = c("decreased","dual","increased"),ylab="Number of unique partners",xlab="")
mtext("specific effects on binding",side=1,line=3)
dev.off()


################################################################################
# FigS11: Numbers of specificity mutations and positions
# while considering only partners with WT binding scores in the linear range
################################################################################

wt_linear = binding_scores$bZ[binding_scores$id == "0xx" & binding_scores$fitness > -1 & binding_scores$fitness < -0.2 & binding_scores$bZ != "ATF4"]
bind_lin = binding_scores[binding_scores$bZ %in% wt_linear,]


# number of pairs with strong positive changes in specificity
pos_pairs = nrow(bind_lin[
  !bind_lin$filt & 
    !is.na(bind_lin$spe_fdr) & 
    bind_lin$spe_fdr < 0.05 & 
    bind_lin$spe > 0.5,
])

# number of pairs with strong negative changes in specificity
neg_pairs = nrow(bind_lin[
  !bind_lin$filt & 
    !is.na(bind_lin$spe_fdr) & 
    bind_lin$spe_fdr < 0.05 & 
    bind_lin$spe < -0.5,
])

# number of unique JUN variants with strong positive changes in specificity
pos_var2 = unique(bind_lin$id[
  !bind_lin$filt & 
    !is.na(bind_lin$spe_fdr) & 
    bind_lin$spe_fdr < 0.05 & 
    bind_lin$spe > 0.5
])

# number of unique JUN variants with strong negative changes in specificity
neg_var2 = unique(bind_lin$id[
  !bind_lin$filt & 
    !is.na(bind_lin$spe_fdr) & 
    bind_lin$spe_fdr < 0.05 & 
    bind_lin$spe < -0.5
])

dual_var2 = pos_var2[pos_var2 %in% neg_var2]
pos_var2 = pos_var2[!(pos_var2 %in% dual_var2)]
neg_var2 = neg_var2[!(neg_var2 %in% dual_var2)]

# number of partners affected by strictly positive specificity mutations
pos_part = unique(bind_lin$bZ[
  bind_lin$id %in% pos_var2 &
    !bind_lin$filt & 
    !is.na(bind_lin$spe_fdr) & 
    bind_lin$spe_fdr < 0.05 & 
    bind_lin$spe > 0.5
])

# number of partners affected by strictly negative specificity mutations
neg_part = unique(bind_lin$bZ[
  bind_lin$id %in% neg_var2 &
    !bind_lin$filt & 
    !is.na(bind_lin$spe_fdr) & 
    bind_lin$spe_fdr < 0.05 & 
    bind_lin$spe < -0.5
])

# number of partners affected by dual specificity mutations
dual_part = unique(bind_lin$bZ[
  bind_lin$id %in% dual_var2 &
    !bind_lin$filt & 
    !is.na(bind_lin$spe_fdr) & 
    bind_lin$spe_fdr < 0.05 & 
    abs(bind_lin$spe) > 0.5
])


pdf("007-figures/FigS11_quantification_specificity_mutations_restricted.pdf")
par(mfrow=c(2,2))
barplot(c(neg_pairs,pos_pairs), col = c("cornflowerblue","red"),names.arg = c("decreased","increased"),ylab="Number of variant:partner pairs",xlab="")
mtext("specific effects on binding",side=1,line=3)
barplot(c(length(neg_var2),length(dual_var2),length(pos_var2)), col = c("cornflowerblue","purple","red"),names.arg = c("decreased","dual","increased"),ylab="Number of unique variants",xlab="")
mtext("specific effects on binding",side=1,line=3)
barplot(c(length(neg_part),length(dual_part),length(pos_part)), col = c("cornflowerblue","purple","red"),names.arg = c("decreased","dual","increased"),ylab="Number of unique partners",xlab="")
mtext("specific effects on binding",side=1,line=3)
dev.off()



################################################################################
# Fig6b: Determinants of specificity
################################################################################

pos_det = unique(substr(neg_var,1,2))
neg_det = unique(substr(pos_var,1,2))
dual_det = unique(substr(dual_var,1,2))
pos_det = sort(c(pos_det,dual_det))
neg_det = sort(c(neg_det,dual_det))

m = cbind(
  as.numeric(unique(substr(allIDs[allIDs != "0xx"],1,2)) %in% neg_det),
  as.numeric(unique(substr(allIDs[allIDs != "0xx"],1,2)) %in% pos_det)
)
m[m==0] = NA
m[,2] = m[,2] + 1

pdf("007-figures/Fig6b_determinants_of_specificity.pdf")
image(1:nrow(m),1:ncol(m),m,col=c("red","cornflowerblue"),axes=F, xlab="", ylab="")
box()
abline(h=1.5)
abline(v=seq(1.5,31.5,1))
axis(2,at=1:2,tick=F,labels = c("-","+"))
mtext("determinant of specificity",2,line=3)
axis(1,at=1:32,tick=F,labels=unique(substr(allIDs[allIDs != "0xx"],1,2)),cex.axis=0.5)
dev.off()


################################################################################
# Fig6c: Specific vs global effects
################################################################################

axes_lim = range(c(binding_scores$glob[!binding_scores$filt & !is.na(binding_scores$glob_fdr)],binding_scores$spe[!binding_scores$filt & !is.na(binding_scores$glob_fdr)]))

pdf("007-figures/Fig6c_specific_vs_global_effects.pdf")
par(mar=c(4,4,4,4))
plot(binding_scores$glob[!binding_scores$filt & binding_scores$rel_fdr >= 0.05 & !is.na(binding_scores$glob_fdr)], # here we use !is.na for glob_fdr instead of rel_fdr because it accounts for the filtering of ATF4
     binding_scores$spe[!binding_scores$filt & binding_scores$rel_fdr >= 0.05 & !is.na(binding_scores$glob_fdr)],
     xlab = "Global effect", ylab='Specific effect',
     pch=".",
     xlim=axes_lim,
     ylim=axes_lim,
)
points(binding_scores$glob[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$glob_fdr) & binding_scores$rel > 0],
       binding_scores$spe[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$glob_fdr) & binding_scores$rel > 0],
       col="red",
       pch=16,
       cex=0.5
)
points(binding_scores$glob[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$glob_fdr) & binding_scores$rel < 0],
       binding_scores$spe[!binding_scores$filt & binding_scores$rel_fdr < 0.05 & !is.na(binding_scores$glob_fdr) & binding_scores$rel < 0],
       col="cornflowerblue",
       pch=16,
       cex=0.5
)
abline(v=0,lty=3)
abline(h=0,lty=3)
abline(0,1,lty=3)
dev.off()


################################################################################
# Fig6d: Violin plot specific vs global effects
################################################################################

tmp = data.frame(non_spe = binding_scores$ddG_mut, spe = binding_scores$spe_fdr < 0.05)
tmp$spe2 = "n.s."
tmp$spe2[!is.na(binding_scores$spe_fdr) & binding_scores$spe_fdr < 0.05 & binding_scores$spe > 0 & binding_scores$spe <= 0.5] = "]0;0.5]"
tmp$spe2[!is.na(binding_scores$spe_fdr) & binding_scores$spe_fdr < 0.05 & binding_scores$spe > 0.5] = "> 0.5"
tmp$spe2[!is.na(binding_scores$spe_fdr) & binding_scores$spe_fdr < 0.05 & binding_scores$spe < 0 & binding_scores$spe >= -0.5] = "[-0.5;0["
tmp$spe2[!is.na(binding_scores$spe_fdr) & binding_scores$spe_fdr < 0.05 & binding_scores$spe < -0.5] = "< -0.5"
tmp$spe2 = factor(tmp$spe2, levels=c("< -0.5","[-0.5;0[","n.s.","]0;0.5]","> 0.5"))

pdf("Fig6D.pdf")
ggplot(tmp, aes(x=spe2, y=non_spe)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), col=2) +
  geom_violin(fill=NA)
dev.off()

table(tmp$spe2)


# number of significant specific and global effects
one.sample.ttest = function(m,s,mu,n){
  tstat = (m-mu)/s
  2*pt(abs(tstat), n-1, lower.tail=FALSE)
}
dG2$pval = one.sample.ttest(dG2$mean, dG2$std, 0, 6)
dG2$fdr = p.adjust(dG2$pval, method="fdr")
nrow(dG2[dG2$fdr < 0.05 & dG2$id2 %in% c(pos_var, neg_var, dual_var),])
nrow(dG2[dG2$fdr < 0.05 & ((dG2$id2 %in% pos_var & dG2$mean < 0) | (dG2$id2 %in% neg_var & dG2$mean > 0)),])
nrow(dG2[dG2$fdr < 0.05 & ((dG2$id2 %in% pos_var & dG2$mean > 0) | (dG2$id2 %in% neg_var & dG2$mean < 0)),])
nrow(dG2[dG2$fdr < 0.05 & dG2$id2 %in% pos_var & dG2$mean > 0,])
