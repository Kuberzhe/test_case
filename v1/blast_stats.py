import pandas as pd

header = ['id', 'acc', 'pid', 'length_hit','mism', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bitsc']
table = pd.read_csv('resp1-1201107674_S85.merged.tsv', sep='\t', names=header)
file_taxid = open('accession2taxid.tsv', 'r').readlines()
ac2taxid = {}
tax_names = ['taxid', 'cellular root', 'acellular root',
'realm', 'domain', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
'subspecies/strain', 'subspecies', 'strain']
df_tax = pd.read_csv('lineage_ref_todf.txt', sep=';', names=tax_names)

for line in file_taxid:
    ac2taxid[line.split()[0]] = line.strip().split()[1]

table['taxid'] = [int(ac2taxid[ac]) for ac in table['acc']]
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
n_total = table['id'].max()	# estimate of the number of total queried reads assuming continuous numbering

for rank in tax_names[1:]:
    for name, sub in table_filtered.groupby(rank):
        res_fin[name] = {f'median_{stat}': sub[stat].median() for stat in ['pid', 'eval', 'length_hit']}
        res_fin[name]['rank'] = rank
        res_fin[name]['n_reads'] = len(sub)
        res_fin[name]['fraction_of_hits'] = len(sub) / n_reads
        res_fin[name]['fraction_of_total_reads_predicted'] = len(sub) / n_total
print('Combining stats and saving them into virus_stats.xlsx')
res_df = pd.DataFrame(res_fin).T.infer_objects()
print(res_df)
res_df.to_excel('virus_stats.xlsx', index_label='taxon')
print('Saving taxonomy into virus_taxonomy.xlsx')
res_tax_df = (table_filtered[tax_names[1:]].drop_duplicates().sort_values(tax_names[1:]).dropna(axis=1, how='all'))
res_tax_df.to_excel('virus_taxonomy.xlsx', index=False)
print('Done')
