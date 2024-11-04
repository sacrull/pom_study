# Install libraries
```R
# install.packages("ggplot2")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("Rhdf5lib")
# BiocManager::install("phyloseq")
# BiocManager::install("microbiome")
# install.packages("remotes")
# remotes::install_github("vmikk/metagMisc")
# remotes::install_github("gauravsk/ranacapa")
```
# 2. Load libraries
```R
library(ggplot2, verbose=F)
library(phyloseq, verbose=F)
library(ape, verbose=F)
library(metagMisc, verbose=F)
library(plyr, verbose=F)
library(dplyr, verbose=F)
library(vegan, verbose=F)
library(ranacapa, verbose=F)
library(microbiome, verbose=F)
library(corncob, verbose=F)
library(magrittr, verbose=F)
library(ggpubr, verbose=F)
library(ecole, verbose=F)
library(UpSetR)
library(smplot2)
setwd("~/pom_study/01-16S-analysis/02-diversity")
```
# 3. Load data
```R
set.seed(578934)
# read in ASV and taxonomy tables previously generated
seqtab <- read.table("../01-read_processing/sequence_table.merged.txt", header=T, row.names=1)
tax <- read.table("../01-read_processing/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
# read in metadata file
map <- read.table("../../meta.txt", sep="\t", header=T, row.names=1)
```
## 3.1 CHeck samples in metadata
```R
notinmeta <- setdiff(row.names(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), row.names(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw
```
## 3.2 Add in PCR blanks
```R
# get current row names
rn <- row.names(map) 
# add in rows for PCR blanks to dataframe
map[nrow(map) + seq_along(notinmeta), ] <- "blank" 
row.names(map) <- c(rn, notinmeta)
# check that it worked
notinmeta <- setdiff(row.names(seqtab), row.names(map))
print("Samples found in ASV table but not in metadata:")
notinmeta
```
## 3.3 Create phyloseq object
```R
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=F), sample_data(map), tax_table(as.matrix(tax)))
ps.dat
```
## 3.4 Remove low abundance ASVs
```R
# compute prevalence dataframe
prevdf <- apply(X=otu_table(ps.dat), MARGIN=ifelse(taxa_are_rows(ps.dat), yes=1, no=2), FUN=function(x){sum(x>0)})
# add taxa and total read counts to dataframe
prevdf <- data.frame(Prevalence=prevdf, TotalAbundance=taxa_sums(ps.dat), tax_table(ps.dat))
# which phyla are comprised as mostly low prevalence ASVs?
lowprev <- ggplot(prevdf, aes(TotalAbundance, Prevalence, nsamples(ps.dat), color="V4")) + geom_hline(yintercept=0.05, alpha=0.5, linetype=2) + geom_point(size=2, alpha=0.7) + scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + facet_wrap(~V4) + theme(legend.position="none")
pdf("./totalabund_vs_prevalence.pdf")
lowprev
dev.off()
# kept asvs must be found in at least 0.1% of all samples (~2) or have a total abundance of 1000 reads or more
ps.dat <- phyloseq_filter_prevalence(ps.dat, prev.trh=0.001, abund.trh=1000, threshold_condition="OR")
ps.dat
```
## 3.5 Remove samples with low read count
```R
# read count summary before filter
summary(rowSums(otu_table(ps.dat)))
ps.dat <- prune_samples(sample_sums(ps.dat) > 800, ps.dat)
ps.dat
# read count after filter
summary(rowSums(otu_table(ps.dat)))
```
## 3.6 Blank analysis
```R
# how many blanks are left after fitering?
blanks <- subset_samples(ps.dat, sample_type=="blank")
blanks <- prune_taxa(taxa_sums(blanks) > 1, blanks)
blanks
rownames(sample_data(blanks))
# two blanks have passed our filtering steps -- top taxa detected?
rel.abund <- transform_sample_counts(blanks, function(x) x/sum(x)) # get relative abundance
glom <- tax_glom(rel.abund, taxrank=rank_names(rel.abund)[6]) # collapse 
data <- psmelt(glom) # create dataframe from phyloseq object
data$V7 <- as.character(data$V7) # convert to character
data$V7[Abundance < 0.05] <- "< 5% abund" # rename low freq phyla
avs <- plyr::ddply(data, ~V7, function(x) c(mean=mean(x$Abundance)))
avs 
pdf("blank_proportion.pdf")
plot_bar(blanks, fill="V3")
dev.off()
system("~/.iterm2/imgcat ./blank_proportion.pdf")
#remove contaminatied samples
# before
ps.dat
# remove possible contaminated samples
remsamps <- c("DM00048V3PQ16", "DM00048V3PQ65", "DM0048V3PQ84", "DM00049V3PQ16", "DM00050V3PQ54", "DM00050V3PQ55", "DM0051V3PQ16", "DM00052V3PQ16", "DM00053V3PQ26-36", "DM00054V3PQ55", "DM00055V3PQ75", "DM00056V3PQ16", "DM00057V3PQ16", "DM00058V3PQ16", "DM00059V3PQ55", "DM00059V3PQ74", "DM00060V3PQ46", "DM00061V3PQ16", "DM00061V3PQ46", "DM00062V3PQ16", "DM00062V3PQ75")
ps.dat <- subset_samples(ps.dat, !(rownames(sample_data(ps.dat)) %in% remsamps))
# now remove blanks
ps.dat <- subset_samples(ps.dat, study_group=="pmw" | study_group=="pmo" | study_group=="pre")
# after
ps.dat
```
## 3.7 Write list of ASVs to keep and filter
```R
write.table(as.data.frame(row.names(tax_table(ps.dat))), file="asvs_to_keep", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
system("/home/$USER/miniconda3/envs/2023-USA-Nigeria/bin/seqtk subseq ../01-read_processing/rep_set.fa asvs_to_keep > rep_set.filt.fa")
system("/home/$USER/miniconda3/envs/2023-USA-Nigeria/bin/mafft --thread 8 rep_set.filt.fa > rep_set.align.fa")
# system("iqtree -s rep_set.align.fa -m MFP --seqtype DNA -T 60 1>iqtree.out 2>iqtree.err") # using modelfinder to select best phylo model
# fasttree neighbor joining tree
system("/home/$USER/miniconda3/envs/2023-USA-Nigeria/bin/fasttree -noml -nt rep_set.align.fa  > rep_set.nj.tre")
tree <- read.tree("rep_set.nj.tre")
ps.dat <- merge_phyloseq(ps.dat, tree)
ps.dat 
```
## 3.8 Save the filtered files and object
```R
write.table(as.data.frame(otu_table(ps.dat)), "sequence_table.filt.txt", sep="\t", row.names=T, col.names=T, quote=F)
# write filtered taxonomy to file
write.table(as.data.frame(tax_table(ps.dat)), "taxonomy_bac.filt.txt", sep="\t", row.names=T, col.names=T, quote=F)
# filtered metadata
write.table(as.data.frame(sample_data(ps.dat)), "map.filt.txt", sep="\t", row.names=T, col.names=T, quote=F)
     
save.image("../master_phyloseq.RData")
# or load saved image
load("../master_phyloseq.RData")
```
# 4. Dataset summary
## 4.1 Sample summary
```R
# colnames(sample_data(ps.dat))
temp <- as.data.frame(cbind(sample_data(ps.dat)$study_group, sample_data(ps.dat)$menopause,sample_data(ps.dat)$pH))
colnames(temp) <- c("study_group", "menopause", "pH")
temp <- ddply(temp, .(study_group), summarise, count = n())
colnames(temp) <- c("study_group", "count")
study_group <- factor(temp$study_group, levels=c("pre", "pmo", "pmw"))
# plot
temp <- na.omit(temp)
pdf("./sample_summary.pdf")
ggplot(temp, aes(x = study_group, y = count, fill = study_group)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme_minimal() +
  coord_flip() +
  labs(x = "Study Group", y = "Count", title = "Counts by Study Group") +
  theme(legend.position = "none")
dev.off()
system("~/.iterm2/imgcat ./sample_summary.pdf")
```
## 4.2 Top phyla in dataset
```R
options(getClass.msg=FALSE) # this is supposed to suppress class warnings but clearly does not work :|
rel.abund <- transform_sample_counts(ps.dat, function(x) x/sum(x)) # get relative abundance
glom <- tax_glom(rel.abund, taxrank=rank_names(rel.abund)[2]) # collapse 
data <- psmelt(glom) # create dataframe from phyloseq object
data$V3 <- as.character(data$V3) # convert to character
data$V3[Abundance < 0.05] <- "< 5% abund" # rename low freq phyla
avs <- plyr::ddply(data, ~V3, function(x) c(mean=mean(x$Abundance)))
avs
#                             V3         mean
# 1               Actinomycetota 1.032982e-01
# 2                    Bacillota 5.801308e-01
# 3                 Bacteroidota 3.915240e-02
# 4             Bdellovibrionota 3.106941e-05
# 5             Campylobacterota 1.237654e-02
# 6  Candidatus_Saccharibacteria 4.540141e-05
# 7               Fusobacteriota 1.050938e-02
# 8               Mycoplasmatota 3.469616e-04
# 9               Pseudomonadota 2.336323e-01
# 10               Spirochaetota 3.715550e-04
# 11                Synergistota 1.198117e-04
# 12         Terrabacteria_group 1.994871e-02
# 13     Thermodesulfobacteriota 3.691501e-05
write.table(avs, "average_abund_phyla.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
```
## 4.3 Taxonomy barchart at phylum level
```R
grouped <- data %>% group_by(study_group, V3) %>% summarize(Abundance = mean(Abundance))
# get correct order of group factors
grouped$study_group <- as.character(grouped$study_group)
grouped$study_group <- factor(grouped$study_group, levels=c("pre", "pmo", "pmw"))
pdf("./grouped_tax_barchart.pdf")
ggplot(grouped, aes(fill=V3, y=Abundance, x=study_group)) + geom_bar(position="fill", stat="identity") + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./grouped_tax_barchart.pdf")
```
## 4.4 Top genera in full dataset
```R
glom <- tax_glom(rel.abund, taxrank=rank_names(rel.abund)[7]) # collapse 
data <- psmelt(glom) # create dataframe from phyloseq object
data$V7 <- as.character(data$V7) # convert to character
data$V7[data$Abundance < 0.45] <- "< 45% abund" # rename low freq phyla
avs <- plyr::ddply(data, ~V7, function(x) c(median=median(x$Abundance)))
avs
write.table(avs, "average_abund_genera.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
#                        V7    median
# 1             < 45% abund 0.0000000
# 2             Gardnerella 0.6270431
# 3              Klebsiella 0.4793605
# 4 Lactobacillales_unknown 0.7807229
# 5           Lactobacillus 0.6631344
# 6     Limosilactobacillus 0.5240257
```
# 5. Beta diversity
```R
utiCols <- c("#8213A0", "#FA78FA", "#40A0FA")
```
## 5.1 CLR transformation
```R
ps.dat.clr <- microbiome::transform(ps.dat, transform="clr", target="OTU")
```
## 5.2 Distance based redundancy analysis
```R
ordcap <- ordinate(ps.dat.clr, "CAP", "euclidean", ~study_group)
# capscale plot by HIV status group
pdf("./bdiv_cap.uti_status.pdf")
plot_ordination(ps.dat.clr, ordcap, "samples", color="study_group") + 
    theme_minimal() + 
    scale_color_manual(values=utiCols)
dev.off()
system("~/.iterm2/imgcat ./bdiv_cap.uti_status.pdf")
permanova_pairwise(otu_table(ps.dat.clr), grp=sample_data(ps.dat.clr)$study_group, method="euclidean") # check signficane
```
## 5.3 Beta dispersion
```R
# first pull sample data from phyloseq object
metadata <- as(sample_data(ps.dat.clr), "data.frame")
# calculate aitchison distance (from CLR transformed data)
clr.dist <- dist(otu_table(ps.dat.clr), method="euclidean")
dispr <- vegan::betadisper(clr.dist, phyloseq::sample_data(ps.dat.clr)$study_group)
print("Beta disperson UTI status")
dispr
permutest(dispr)
```
# 6. Alpha diversity
```R
pdf("./adiv.study_group.pdf")
plot_richness(ps.dat, measures=c("Observed", "Shannon"), x="study_group") + 
    theme_minimal() + 
    geom_jitter(alpha=0.25) +
    geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
    geom_point() +
    stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 23,
    fill = "red"
  )
dev.off()
system("~/.iterm2/imgcat ./adiv.study_group.pdf")
```
## 6.1 Rarefaction
```R
options(warn=-1) # suppress warnings
p <- ggrare(ps.dat, step = 1000, color = "study_group", se = TRUE)
p + theme_minimal()
pdf("./rarefaction_plots.pdf")
p + theme_minimal()
dev.off()
system("~/.iterm2/imgcat ./rarefaction_plots.pdf")
options(warn=0) # back on
```