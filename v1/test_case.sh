# Dependencies - curl, python, taxonkit >= 0.20 (with ncbi taxdump in ~/.taxonkit), pandas, openpyxl (for excel conversion)

# first major goal is to go from accesion of the hits to the full taxonomy. First we get accessions
cut -f 2 resp1-1201107674_S85.merged.tsv | sort | uniq > accession.txt

# Now is the slowest step - getting taxid for each accesion. Usually this can be done directly from the blastdb which was
# used to create the result, but since I do not have it readily avaliable, we will go to ncbi directly
xargs -n 100 < accession.txt | while read -r batch
do
	curl -s --retry 5 -d "db=nuccore&retmode=json&id=${batch// /,}" https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi \
  | grep -oE '"(taxid|accessionversion)": *("[^"]*"|[0-9]+)' | sed -E 's/^"[a-z]+": *"?//; s/"$//'\
  | awk -v OFS='\t' 'NR%2 {t=$0; next} {print $0, t}'
	sleep 1s
done > accession2taxid.tsv


cut -f 1,2 resp1-1201107674_S85.merged.tsv | sort | uniq | wc -l
#4604
cut -f 1 resp1-1201107674_S85.merged.tsv | sort | uniq | wc -l
#4604
echo 'Every query has at most 1 sequence in hits, meaning we do not need to deal with assigning reads with multiple hits on different organisms'

# now we need to go from taxid to taxonomy. One of the tools that does that is taxonkit. If we had more than one organism
# hit per query we could use something like blobtools to classify each read
taxonkit lineage -i 2 accession2taxid.tsv > lineage.txt
taxonkit reformat -I 2 -f "{C};{a};{r};{d};{k};{K};{p};{c};{o};{f};{g};{s};{t};{S};{T}" lineage.txt > lineage_ref.txt
cut -f 2,4 lineage_ref.txt | sort | uniq | sed -e 's/\t/;/g' > lineage_ref_todf.txt

# now calculate statistics and produce a xlsx file with the results.
python3 blast_stats.py
