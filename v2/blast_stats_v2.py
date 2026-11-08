import pandas as pd
from Bio import Entrez
Entrez.email = "yakovkorobitcyn@gmail.com"

header = ['id', 'acc', 'pid', 'length_hit','mism', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bitsc']
table = pd.read_csv('resp1-1201107674_S85.merged.tsv', sep='\t', names=header)

accs = sorted(table['acc'].unique())

# accession to taxid, instead of taxonkit
acc2taxid = {}
for i in range(0, len(accs), 200):
    with Entrez.esummary(db="nuccore", id=",".join(accs[i:i + 200])) as h:
        for rec in Entrez.read(h):
            acc2taxid[rec["AccessionVersion"]] = int(rec["TaxId"])

# taxid to ranks
taxids = set(acc2taxid.values())
ranks = {}
with Entrez.efetch(db="taxonomy", id=",".join(map(str, taxids)), retmode="xml") as h:
    for rec in Entrez.read(h):
        r = {x["Rank"]: x["ScientificName"] for x in rec["LineageEx"]
             if x["Rank"] not in ("no rank", "clade")}
        r[rec["Rank"]] = rec["ScientificName"]
        ranks[int(rec["TaxId"])] = r

df_tax = pd.DataFrame.from_dict(ranks, orient="index").rename_axis("taxid").reset_index()
tax_names = [c for c in df_tax.columns if c not in ('no rank', 'clade')]


table['taxid'] = table['acc'].map(acc2taxid)
table = table.merge(df_tax, on='taxid', how="left")
# for all the stat calculations we will take one hit per query, the one with the highest bitscore
table_filtered = table.loc[table.groupby(["id", "acc"])["bitsc"].idxmax()]
table_filtered = table_filtered[table_filtered['acellular root'] == 'Viruses']
for tax_gr in tax_names[1:]:
    if table_filtered[tax_gr].isna().all():     # skipping columns without any data
        continue
    for stat in ['pid', 'eval', 'length_hit']:
        print(f'{stat} for {tax_gr}')
        print(table_filtered.groupby(tax_gr)[stat].median())

res_fin = {}

n_reads = table['id'].nunique()  # number of queried sequences that found at least 1 hit
n_total = table['id'].max()     # estimate of the number of total queried reads assuming continuous numbering

for rank in tax_names[1:]:
    for name, sub in table_filtered.groupby(rank):
        res_fin[name] = {f'median_{stat}': sub[stat].median() for stat in ['pid', 'eval', 'length_hit']}
        res_fin[name]['rank'] = rank
        res_fin[name]['n_reads'] = len(sub)
        res_fin[name]['fraction_of_hits'] = len(sub) / n_reads
        res_fin[name]['fraction_of_total_reads_predicted'] = len(sub) / n_total
print('Combining stats and saving them into virus_stats_v2.xlsx')
res_df = pd.DataFrame(res_fin).T.infer_objects()
print(res_df)
res_df.to_excel('virus_stats_v2.xlsx', index_label='taxon')
print('Saving taxonomy into virus_taxonomy_v2.xlsx')
res_tax_df = (table_filtered[tax_names[1:]].drop_duplicates().sort_values(tax_names[1:]).dropna(axis=1, how='all'))
res_tax_df.to_excel('virus_taxonomy_v2.xlsx', index=False)
print('Done')
