# 1. Post w/o UTi vs Post w/ UTI
```R
library(coda4microbiome)
library(tidyverse)
library(phyloseq)
library(grid)

setwd("~/pom_study/01-16S-analysis/03-balance_taxa")
load("../master_phyloseq.RData")
set.seed(12349)

# collapse data to roughly species level to minimize high sparsity
sub.dat <- subset_samples(ps.dat , study_group == "pmo" | study_group == "pmw")
glom <- tax_glom(sub.dat, "V8")
# remove any taxa with fewer than 50 counts and in at least 10% of samples post merging
# glom <- filter_taxa(glom, function(x) sum(x > 10) > (0.10*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)

#get taxonomy of ASV
taxa = as(tax_table(ps.dat), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V8)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")

#rename to have V8 level name
dat <- as.data.frame(dat)
dat <- dat %>% rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))
rownames(dat) <- dat$V8
dat <- dat[2:(length(dat)-1)] #remove last column
dat <- as.matrix(t(dat))
# merge metadata with asv table so response variable in same order
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,1:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
# response variable
dif2 <- dim(dat)[2] - 4
y <- factor(dat$study_group) #geog lcaotion

#running coda
coda <- coda_glmnet(x=x,y=y) 
sum(coda$`log-contrast coefficients`)
#positive taxa
coef<-coda$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
coda$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
coda$taxa.name[negatives[on]]
pdf("./post.uti.pdf")
coda$`signature plot`
coda$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./post.uti.pdf")
