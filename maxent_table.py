# %% 
import pandas as pd
import numpy as np
import sys

sys.path.append('MaxEntScorepy')
from maxent import *
# %%

csv_file = pd.read_csv('41467_2019_12642_MOESM10_ESM.csv.gz', index_col= 'library index', compression='gzip')
# intron start - 5' - donor
matrix5 = load_matrix5()
istartmaxent5 = pd.DataFrame(index = csv_file.index, columns= range(len(csv_file.sequence.values[0])+1 )  )
for chr in istartmaxent5.index:
    for intpoint in istartmaxent5.columns:
        sequence = str(csv_file.loc[chr, 'sequence'])[int(intpoint) - 3: int(intpoint) + 6] # 9 nucleotides: 3 bases in exon and 6 bases in intron 
        if len(sequence) == 9:
            # istartmaxent5.loc[chr, intpoint] = maxent.score5(str(sequence), matrix=matrix5)
            istartmaxent5.loc[chr, intpoint] = score5(str(sequence), matrix=matrix5)
        else:
            istartmaxent5.loc[chr, intpoint] = np.NaN
# intron end - 3' - acceptor
matrix3 = load_matrix3()
iendmaxent3 = pd.DataFrame(index = csv_file.index, columns= range(len(csv_file.sequence.values[0])+1 )  )
for chr in iendmaxent3.index:
    for intpoint in iendmaxent3.columns:
        sequence = str(csv_file.loc[chr, 'sequence'])[int(intpoint) - 20: int(intpoint) + 3] # 23 nucleotides: 20 bases in the intron and 3 base in the exon
        if len(sequence) == 23:
            # iendmaxent3.loc[chr, intpoint] = maxent.score3(str(sequence), matrix=matrix3)
            iendmaxent3.loc[chr, intpoint] = score3(str(sequence), matrix=matrix3)
        else:
            iendmaxent3.loc[chr, intpoint] = np.NaN

istartmaxent5.to_csv('ce_istartmaxent5.csv.gz',header=True, index=True, sep='\t', compression='gzip')
iendmaxent3.to_csv('ce_iendmaxent3.csv.gz',header=True, index=True, sep='\t', compression='gzip')
