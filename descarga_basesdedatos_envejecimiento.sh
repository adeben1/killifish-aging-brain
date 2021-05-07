# Download databases (done at Apr 2021)

#wget https://genomics.senescence.info/genes/models_genes.zip
#wget https://genomics.senescence.info/genes/human_genes.zip
#wget https://genomics.senescence.info/cells/cellAge.zip
#wget https://genomics.senescence.info/drugs/dataset.zip
#wget https://www.dgidb.org/data/monthly_tsvs/2020-Nov/interactions.tsv 

unzip \*.zip
rm *zip

# Filter drug database to keep only the significant treatments
cut -f1,8 -d "," drugage.csv | awk -F ',' '$2 == "S"' | cut -f1 -d "," | sort -u > uniquedrugsaging.csv

