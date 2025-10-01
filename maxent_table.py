# %% 
import pandas as pd
import numpy as np
from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3
import argparse
# %%
# argspars arguments
parser = argparse.ArgumentParser(description='create 2 maxent score matrices: one for donors strength scores, one for acceptor strength scores')
parser.add_argument('-c', '--csv', type=str, help='path to library metadata')
parser.add_argument('-o', '--out',type=str, default='ce', help = 'output prefix')
args = parser.parse_args()

csv_file = pd.read_csv(args.csv, index_col= 'library index', sep='\t')
# intron start - 5' - donor
matrix5 = load_matrix5()
istartmaxent5 = pd.DataFrame(index = csv_file.index, columns= range(len(csv_file.sequence.values[0])+1 )  )
for chr in istartmaxent5.index:
    for intpoint in istartmaxent5.columns:
        sequence = str(csv_file.loc[chr, 'sequence'])[int(intpoint) - 3: int(intpoint) + 6] # 9 nucleotides: 3 bases in exon and 6 bases in intron 
        if len(sequence) == 9:
            istartmaxent5.loc[chr, intpoint] = maxent.score5(str(sequence), matrix=matrix5)
        else:
            istartmaxent5.loc[chr, intpoint] = np.NaN
# intron end - 3' - acceptor
matrix3 = load_matrix3()
iendmaxent3 = pd.DataFrame(index = csv_file.index, columns= range(len(csv_file.sequence.values[0])+1 )  )
for chr in iendmaxent3.index:
    for intpoint in iendmaxent3.columns:
        sequence = str(csv_file.loc[chr, 'sequence'])[int(intpoint) - 20: int(intpoint) + 3] # 23 nucleotides: 20 bases in the intron and 3 base in the exon
        if len(sequence) == 23:
            iendmaxent3.loc[chr, intpoint] = maxent.score3(str(sequence), matrix=matrix3)
        else:
            iendmaxent3.loc[chr, intpoint] = np.NaN

istartmaxent5.to_csv(header=True, index=True, sep='\t', compression='gzip')
iendmaxent3.to_csv(header=True, index=True, sep='\t', compression='gzip')