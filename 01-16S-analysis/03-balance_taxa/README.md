# 1. Post w/ vs w/o at species level
```R
library(grid)
library(coda4microbiome)
library(tidyverse)
library(phyloseq)
library(corncob)

setwd("~/pom_study/01-16S-analysis/03-balance_taxa")
load("~/pom_study/01-16S-analysis/master_phyloseq.RData")
set.seed(12349)

ps.sub <- subset_samples(ps.dat, menopause == "post")
# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(ps.sub, "V8")
# remove any taxa with fewer than 50 counts and in at least 05% of samples post merging
glom <- filter_taxa(glom, function(x) sum(x > 50) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(glom)))
map <- sample_data(glom)

#get taxonomy of ASV
taxa = as(tax_table(ps.sub), "matrix")
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
# dif2 <- dim(dat)[2] - 4
y <- factor(dat$group) #geog lcaotion
# z <- data.frame(Tooth_Classification = as.factor(dat$Tooth_Classification)) #possible cofound


model <- coda_glmnet(x=x,y=y)

sum(model$`log-contrast coefficients`)

#positive taxa
coef<-model$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
model$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
model$taxa.name[negatives[on]]

pdf("./pm.ruti.pdf")
model$`signature plot`
model$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./pm.ruti.pdf")


#diff abundance
da_analysis <- differentialTest(formula = ~ group,
                               phi.formula = ~ group,
                               formula_null = ~ 1,
                               phi.formula_null = ~ group,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.sub,
                               fdr_cutoff = 0.05)
da_analysis
#look at sign taxa
da_analysis$significant_taxa
pdf("./diffab.post.uti.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.post.uti.pdf")
```
# 2. Post w/ vs w/o at ASV level
```R
library(grid)
library(coda4microbiome)
library(tidyverse)
library(phyloseq)
library(corncob)

setwd("~/pom_study/01-16S-analysis/03-balance_taxa")
load("~/pom_study/01-16S-analysis/master_phyloseq.RData")
set.seed(12349)

ps.sub <- subset_samples(ps.dat, menopause == "post")
# collapse data to roughly species level to minimize high sparsity
# glom <- tax_glom(ps.sub, "V8")
# remove any taxa with fewer than 50 counts and in at least 05% of samples post merging
ps.sub <- filter_taxa(ps.sub, function(x) sum(x > 50) > (0.05*length(x)), TRUE)
# pull data
dat <- t(as.data.frame(otu_table(ps.sub)))
map <- sample_data(ps.sub)

#get taxonomy of ASV
taxa = as(tax_table(ps.sub), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V8)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")

#rename to have V8 level name
dat <- as.data.frame(dat)
dat <- dat %>% rownames_to_column(var = "ASV")
dat <- left_join(dat, orderdf, by=c('ASV'='ASV'))
rownames(dat) <- paste0(dat$V8, "_", dat$ASV)
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
# dif2 <- dim(dat)[2] - 4
y <- factor(dat$group) #geog lcaotion
# z <- data.frame(Tooth_Classification = as.factor(dat$Tooth_Classification)) #possible cofound


model <- coda_glmnet(x=x,y=y)

sum(model$`log-contrast coefficients`)

#positive taxa
coef<-model$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
model$taxa.name[positives[op]]
#negative taxa
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
model$taxa.name[negatives[on]]

pdf("./post.uti.asv.pdf")
model$`signature plot`
model$`predictions plot`
dev.off()
system("~/.iterm2/imgcat ./post.uti.asv.pdf")

#diff abundance
da_analysis <- differentialTest(formula = ~ group,
                               phi.formula = ~ group,
                               formula_null = ~ 1,
                               phi.formula_null = ~ group,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.sub,
                               fdr_cutoff = 0.05)
da_analysis
#look at sign taxa
da_analysis$significant_taxa
pdf("./diffab.post.uti.asv.pdf", width = 20)
plot(da_analysis, level=c("V8"))
dev.off()
system("~/.iterm2/imgcat ./diffab.post.uti.asv.pdf")
