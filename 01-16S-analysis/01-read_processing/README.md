# 0. Combine runs
```sh
cd ~/pom_study/01-16S-analysis/data
ls run1 | sed 's/_.*//' | sort | uniq | while IFS= read -r line; do
  cat ./run1/"$line"*R1*.fastq.gz ./run2/"$line"*R1*.fastq.gz ./run3/"$line"*R1*.fastq.gz ./run4/"$line"*R1*.fastq.gz > "$line"_R1.fastq.gz;
done

ls run1 | sed 's/_.*//' | sort | uniq | while IFS= read -r line; do
  cat ./run1/"$line"*R2*.fastq.gz ./run2/"$line"*R2*.fastq.gz ./run3/"$line"*R2*.fastq.gz ./run4/"$line"*R2*.fastq.gz > "$line"_R2.fastq.gz;
done
```
# 1. DADA2 Pipeline
```R
#Load packages
library(dada2, verbose = FALSE)
library(stringr, verbose = FALSE)
library(data.table, verbose = FALSE)
library(ShortRead, verbose = FALSE)
library(Biostrings, verbose = FALSE)
library(seqinr, verbose = FALSE)
library(qualpalr, verbose = FALSE)
set.seed(12349)
#Set up paths
setwd("/home/suzanne/pom_study/01-16S-analysis/01-read_processing")
rawpath <- "/home/suzanne/pom_study/01-16S-analysis/data"
fnFs <- sort(list.files(rawpath, pattern="_R1.fastq.gz", full.names=T))
fnRs <- sort(list.files(rawpath, pattern="_R2.fastq.gz", full.names=T))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names, 50)
paste("Number of input samples: ", length(sample.names))
#Make quality plot scores
system("mkdir img") # ignore warning
fwdqual <- plotQualityProfile(fnFs[10:25])
revqual <- plotQualityProfile(fnRs[10:25])

pdf(paste("img/", "forward_quality_plot.pdf", sep=""))
fwdqual
dev.off()
pdf(paste("img/", "reverse_quality_plot.pdf", sep=""))
revqual
dev.off()
system("~/.iterm2/imgcat ./img/forward_quality_plot.pdf")
system("~/.iterm2/imgcat ./img/reverse_quality_plot.pdf")

#Filter out uncalled bases
fnFs.filtN <- file.path(rawpath, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(rawpath, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)
#Primer Removal
cutadapt <- as.character(system("which cutadapt", intern=T))
system("cutadapt --version")
path.cut <- file.path(rawpath, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc("GTGCCAGCMGCCGCGGTAA")
REV.RC <- dada2:::rc("GGACTACHVGGGTWTCTAAT")
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", "GTGCCAGCMGCCGCGGTAA", "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", "GGACTACHVGGGTWTCTAAT", "-A", FWD.RC)
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("--cores=0", R1.flags, R2.flags, "-n", 2,"-o", fnFs.cut[i], "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))
}
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE))
#Quality Filter and Trim Reads
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, minLen = c(240,175),maxN=c(0,0), maxEE=c(2,2), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$precentage_retained <- retained$reads.out/retained$reads.in * 100
retained
# Learn error rates
errF <- learnErrors(filtFs, multithread=T, random=T)
errR <- learnErrors(filtRs, multithread=T, random=T)

err.f.plt <- plotErrors(errF, nominalQ=TRUE) 
pdf(paste("./img/", "error_plot.pdf", sep=""))
err.f.plt
dev.off()
system("~/.iterm2/imgcat ./img/error_plot.pdf")
# Dereplicating
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
# Sample infrence
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, verbose = FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, verbose = FALSE)
# Merge Sample reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)
# Generate sequence tables
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Check Length distrubtion 
table(nchar(colnames(seqtab)))
# get histogram of length distribution after filter
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab))))
len.plt <- plot(x=length.histogram[,1], y=length.histogram[,2])
pdf(paste("./img/", "length_hist.pdf", sep=""))
plot(x=length.histogram[,1], y=length.histogram[,2])
dev.off()
system("~/.iterm2/imgcat ./img/length_hist.pdf")
# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
# Processing summary
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochimeras")
rownames(track) <- sample.names
track
# write to file
write.table(data.frame("row_names"=rownames(track),track),"read_retention.txt", row.names=FALSE, quote=F, sep="\t")
uniquesToFasta(seqtab.nochim, "rep_set.fa")
# fix ASV names 
system("awk '/^>/{print \">ASV\" ++i; next}{print}' < rep_set.fa > rep_set_fix.fa")
system("mv rep_set_fix.fa rep_set.fa")
# write sequence table to file, fix ASV names
my_otu_table <- t(as.data.frame(seqtab.nochim)) 
ASV.seq <- as.character(unclass(row.names(my_otu_table))) 
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') 
colnames(seqtab.nochim) <- ASV.num 
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table.merged.txt", row.names=FALSE, quote=F, sep="\t")
save.image("dada2.RData")
```
Run dada2 script
```sh
Rscript dada2.R 1> dada2.out 2> dada2.err &
disown %1
```
# 3. Assign taxonomy
```sh
kraken2 --db ~/refdb/ezbio/kraken_ezbio \
  --threads 6 \
  --use-names \
  --output rep_set.kraken.out rep_set.fa \
  --unclassified-out rep_set.unclassified.out --confidence 0.01

awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > taxids

#get lineage file
sed "s/\t//g" ~/refdb/ezbio/kraken_ezbio/taxonomy/names.dmp > rankedlineage.dmp
sort -t "|" -k 1b,1 rankedlineage.dmp > rankedlineage_sorted
cat rankedlineage_sorted | sed 's/|\{2,\}/|/g' | sed 's/|/\t/' | sed 's/|.*//'> rankedlineage_clean

grep ">" ~/refdb/ezbio/ezbiocloud_clean.fa | sed 's/_.*\t/\t/' | sed 's/>//' | sed 's/_Arch/Arch/'> ezbiocloud_id_taxonomy.txt
python3 ranked_create.py
sed -i '1i\1\troot|' rankedlineage_clean2
sed -i '1i\0\tunclassified|' rankedlineage_clean2

python3 lineages.py #output is lineage
awk -F"\t" '{print $2}' lineage  | sed 's/ /_/g' | sed 's/|/;/g' | sed 's/;$//' > taxonomy
# merge assemebleies ids and taxonomy
awk '{print $2}' rep_set.kraken.out > asvids
paste asvids taxonomy > taxonomy.txt
#filter out what was only assigned at phylum level
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python3 ~/bin/fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```

# For vince
```sh
#fixing headers
paste -d ' ' <(grep '^>' rep_set.filt.fa) <(awk '{print $8}' taxonomy_bac.txt) > combined.txt
awk '{print $1"_"$2}' combined.txt > headers.txt
paste -d "\t" <(grep '^>' rep_set.filt.fa) headers.txt | sed 's/>//g '> names.txt
seqkit replace -p "(.+)" -r '{kv}|$1' -k names.txt rep_set.filt.fa |sed 's/|.*//' > pom_species.fa
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' pom_species.fa > temp
mv temp pom_species.fa

cp sequence_table.merged.txt pom_table.txt
grep ASV pom_species.fa | sed 's/>//' | while read line; do
    modified_line=$(echo "$line" | sed 's/_.*//')
    sed -i "s/\b$modified_line\b/$line/g" pom_table.txt
done

python3 underscore.py pom_table.txt > temp
mv temp pom_table.txt
```