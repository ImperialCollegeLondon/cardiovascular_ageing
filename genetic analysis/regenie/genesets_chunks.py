'''
usage:
python3 geneset_chunks.py [size_of_chunk] [, separated list of chroms] [out_dir] 
author: changlubio
'''

import pandas as pd
import os
import sys
import subprocess

chunk = 500
if len(sys.argv) > 1:
    chunk = int(sys.argv[1])
chroms = 'all'
if len(sys.argv) > 2:
    chroms = sys.argv[2]
out_dir = '/mnt/project/Derived/Exomes/helper_files/'


path_to_450kwes_helper_files = "/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - interim 450k release/helper_files"
set_file = os.path.join(path_to_450kwes_helper_files, 'ukb23149_450k_OQFE.sets.txt.gz')
geneset_all = pd.read_csv(set_file, sep='\t', header=None)

if chroms == 'all':
    chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X','Y']
else:
    chroms = chroms.split(',')

for chrom in chroms:
    chrom_genes = geneset_all[geneset_all[1]==chrom].sort_values(by=2)
    Nchunk, res = divmod(chrom_genes.shape[0], chunk)
    for i in range(Nchunk+1):
        outpath = 'ukb23149_450k_OQFE_c{}_n{}.sets.tsv'.format(chrom, i)
        chrom_genes.iloc[i*chunk:(i+1)*chunk,].to_csv(outpath, sep='\t', header=None, index=None)

subprocess.call('dx upload *sets.tsv --path /Derived/Exomes/helper_files/',shell=True)
