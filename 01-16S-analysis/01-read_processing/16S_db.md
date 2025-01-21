# 1. Format for kraken2 database
```sh
cd ~/refdb/ezbio
#fixing headers
paste -d '__' <(grep '^>' ezbiocloud_qiime_full.fasta) <(awk -F "\t" '{print $2}' ezbiocloud_id_taxonomy.txt | sed 's/ /_/g') > combined.txt
awk -F "__" '{print $1"_"$2}' combined.txt > headers.txt

paste -d "\t" <(sed 's/_.*;/_/' headers.txt | sed 's/_$//' ) headers.txt | sed 's/>//g' | sed 's/\t.*_B/__B/' | sed 's/\t.*_A/___A/'> temp
mv temp headers.txt
paste -d "\t" <(grep '^>' ezbiocloud_qiime_full.fasta) headers.txt | sed 's/>//g '> names.txt

seqkit replace -p "(.+)" -r '{kv}|$1' -k names.txt ezbiocloud_qiime_full.fasta |sed 's/|.*//' > ezbiocloud_clean.fa
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' ezbiocloud_clean.fa > temp
mv temp ezbiocloud_clean.fa
head ezbiocloud_clean.fa
sed -i 's/_$//' ezbiocloud_clean.fa
sed -i 's/__/\t/' ezbiocloud_clean.fa

sed 's/Bacteria;/\tLineage=Root,rootrank,Bacteria,domain,/g' ./ezbiocloud_clean.fa | sed 's/_Archaea;/\tLineage=Root,rootrank,Archaea,domain,/g' | sed 's/;/,phylum,/' | sed 's/;/,class,/' | sed 's/;/,order,/' | sed 's/;/,family,/' | sed 's/;/,genus,/' | sed '/^>/s/$/,species/' > ./ezbio.fa
sed -i 's/,/;/g' ezbio.fa
sed -i 's/;species;species/;species/' ezbio.fa
sed -i 's/\t\t/\t/' ezbio.fa
sed -i 's/\t_A/\tA/' ezbio.fa
```
# 2. Make kraken2 db
```sh
perl build_rdp_taxonomy.pl ezbio.fa

mkdir -p ./kraken_ezbio/library 
mv ./ezbio.fa ./kraken_ezbio/library/ezbio.fna
mkdir -p ./kraken_ezbio/taxonomy
mv names.dmp nodes.dmp ./kraken_ezbio/taxonomy
mv seqid2taxid.map ./kraken_ezbio
#building database
kraken2-build --build --db ./kraken_ezbio --threads 6
```