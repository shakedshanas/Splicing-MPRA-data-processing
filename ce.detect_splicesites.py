# %%
import pandas as pd
import numpy as np
import os
import argparse
from pandas.api.types import is_float_dtype # only for the row with is_float_dtype
from sys import exit
import time 
import re
# %% 
''' 
last edited on : 2 Jan 2025
14 NOV 2024 added: used super bed file - contains in addition the cigar string and the sequence
24 NOV 2024 added: improvedthe maxent scores, added underscore count and number of exons based on the underscore count
02 FEB 2024 added: changed the maxent function to be as in the ir.detect_splicesites.py script 
command
2 JAN 2025 added: filter out StartLeftBlock that are above 600. because then it shows that some middle exons are really large and that's just wrong. 
# option 1
for b in /mnt/beegfs/home/martinmikl/sshanas/ir/martin/m19/multiplxpy_subread/bedfiles/*.bed ; do input=$(echo "$b" | cut -d '/' -f11 | cut -d '.' -f1) ; sbatch -J $input -o $input.out -e $input.err --wrap "python rmdel_withbedfiles.py --chr $input" && sleep 15 ; done
# option 2 - opdated at May 29th 
for b in /lustre1/home/martinmikl/sshanas/ir/martin/ir_er/IR_Tun_6h/bedfiles/*.bed ; do input=$(echo "$b" | cut -d '/' -f11 | cut -d '.' -f1) && echo "$input"; sbatch -J $input -o $input.out -e $input.err --wrap "python rmdel_withbedfiles.py --chr $input" && sleep 15 ; done
copy is uploaded to git up on 30 Oct 2023


for b in /lustre1/home/martinmikl/sshanas/ce/AdmeraAugust24/CASLIBa/bedfiles/*.bed ; do input=$(echo "$b" | cut -d '/' -f10 | cut -d '.' -f1) ; sbatch -J $input -o $input.out -e $input.err --wrap "python splicesitedetection.ce.py --chr $input" && sleep 5 ; done
updated for 10 Oct 2024

'''
# argspars arguments
parser = argparse.ArgumentParser(description='create a df file for all observed splicing events - if using genomic data, it will be without filetring. otherwise - genomic filtering should be used. ')
parser.add_argument('--chr', type=str, help='chromosome')
parser.add_argument('-p', '--pwrdir', type=str, help='path to put all the putouts')
args = parser.parse_args()

# %%
''' upload files '''
''' in cluster '''
chr = args.chr
print('library seq: ' + str(chr))
transcript_path = args.pwrdir + '/'

if os.path.exists(transcript_path + 'bedfiles/' + 'supbedfile.' + str(chr) + '.bed'):
    bed_file = pd.read_csv(transcript_path + 'bedfiles/' + 'supbedfile.' + str(chr) + '.bed', sep='\t', header=None, 
                           names=['chr', 'StartLeftBlock', 'EndRightBlock', 'Junction_Name', 'Strand','nBlocks', 'BlockSizes', 'BlockStarts', 'cigar', 'sequence'])
t_nblocks_set = set(bed_file.nBlocks)
# print(bed_file.head())
print(t_nblocks_set)
lib = pd.read_csv('~/ce/41467_2019_12642_MOESM10_ESM.csv', index_col='library index')
imaxent5 = pd.read_csv('~/ce/ce_istartmaxent5.csv', sep='\t', index_col='library index')
imaxent3 = pd.read_csv('~/ce/ce_iendmaxent3.csv', sep='\t', index_col='library index')

imaxent5.columns = imaxent5.columns.astype(float).astype(int)
imaxent3.columns = imaxent3.columns.astype(float).astype(int)

# %%
''' functions '''

def extract_introns_lengths(df, nblocks_set):
    if max(nblocks_set) > 1:               

        BlockSizes_df = df.BlockSizes.str.split(',', expand=True).rename(columns = {i: f'BlockSize{i+1}' for i in range(df.nBlocks.max())}).fillna(0).astype(int)
        BlockStarts_df = df.BlockStarts.str.split(',', expand=True).rename(columns = {i: f'BlockStart{i+1}' for i in range(df.nBlocks.max())}).fillna(0).astype(int)
        df = BlockSizes_df.join(BlockStarts_df).join(df[['chr', 'StartLeftBlock', 'nBlocks', 'cigar', 'sumBlockSize', 'deletions']])

        df_fwd = df.iloc[::2]
        df_rev = df.iloc[1::2]
        
        
        df_fwd.reset_index(inplace=True)
        df_rev.reset_index(inplace=True)
        df_fwd = df_fwd.drop('index', axis='columns')
        df_rev = df_rev.drop('index', axis='columns')

        df_rev.columns = [str(x) + '_rev' for x in df_rev.columns]
        df_fwd.columns = [str(x) + '_fwd' for x in df_fwd.columns]
        df2 = df_fwd.join(df_rev)
        
    else: # if there is a single isoform (usually will be unspliced)
        
        df.drop(['BlockSizes', 'BlockStarts'], axis='columns', inplace=True)

        df_fwd = df.iloc[::2]
        df_rev = df.iloc[1::2]

        df_fwd.reset_index(inplace=True)
        df_rev.reset_index(inplace=True)
        df_fwd = df_fwd.drop('index', axis='columns')
        df_rev = df_rev.drop('index', axis='columns')
        
        df_rev.columns = [str(x) + '_rev' for x in df_rev.columns]
        df_fwd.columns = [str(x) + '_fwd' for x in df_fwd.columns]
        df2 = df_fwd.join(df_rev)
            
    for f in ['fwd', 'rev']:
        if (f'nBlocks_{f}' in df2.columns) and (f'StartLeftBlock_{f}' in df2.columns):
            if f'BlockSize1_{f}' in df2.columns:                                                                                 
                df2[f'IStart1_{f}'] = df2.apply(lambda x: np.where(x[f'nBlocks_{f}'] > 1, x[f'StartLeftBlock_{f}'] + x[f'BlockSize1_{f}'], 0) , axis = 1) 
            if f'BlockStart2_{f}' in df2.columns: 
                df2[f'IEnd1_{f}'] = df2.apply(lambda x:   np.where(x[f'nBlocks_{f}'] > 1, x[f'StartLeftBlock_{f}'] + x[f'BlockStart2_{f}'], 0) , axis = 1)
            if f'BlockStart3_{f}' in df2.columns:
                df2[f'IStart2_{f}'] = df2.apply(lambda x: np.where(x[f'nBlocks_{f}'] > 2, x[f'StartLeftBlock_{f}'] + x[f'BlockStart2_{f}'] + x[f'BlockSize2_{f}'], 0), axis = 1)
            if f'BlockStart3_{f}' in df2.columns:
                df2[f'IEnd2_{f}'] = df2.apply(lambda x:   np.where(x[f'nBlocks_{f}'] > 2, x[f'StartLeftBlock_{f}'] + x[f'BlockStart3_{f}'], 0), axis = 1)            
            if (f'IEnd1_{f}' in df.columns) and (f'IStart1_{f}' in df.columns): 
                df[f'Ilen1_{f}'] = df[f'IEnd1_{f}'] - df[f'IStart1_{f}']
            if (f'IEnd2_{f}' in df.columns) and (f'IStart2_{f}' in df.columns): 
                df[f'Ilen2_{f}'] = df[f'IEnd2_{f}'] - df[f'IStart2_{f}']
            if (f'IEnd3_{f}' in df.columns) and (f'IStart3_{f}' in df.columns): 
                df[f'Ilen3_{f}'] = df[f'IEnd3_{f}'] - df[f'IStart3_{f}']

    return df2

def makeisoform(pd_row):
    l_fwd = [int(pd_row[col]) for col in ['I' + point + str(n) + '_' + 'fwd' for n in range(1,max(nblocks_set_fwd)) for point in ['Start', 'End']] if col in bed_file.columns]
    l_rev = [int(pd_row[col]) for col in ['I' + point + str(n) + '_' + 'rev' for n in range(1,max(nblocks_set_rev)) for point in ['Start', 'End']] if col in bed_file.columns]
    l = list(set(l_fwd + l_rev))
    l.sort()
    return str(chr) + '_' + '_'.join([str(e) for e in l if e > 0])

# %%
bed_file['BlockSizes'] = bed_file['BlockSizes'].astype(str)
bed_file = bed_file.join(pd.DataFrame(bed_file.BlockSizes.str.split(',', expand=True).fillna(0).astype(int).sum(axis='columns'))) .rename(columns={0:'sumBlockSize'})
bed_file = bed_file.assign(deletions = bed_file.apply(lambda x: sum([int(n) for n in re.findall('[A-Z](\d+)D', x.cigar)]), axis='columns') )

bed_file = extract_introns_lengths(bed_file, t_nblocks_set)
nblocks_set_fwd = set(bed_file.nBlocks_fwd)
nblocks_set_rev = set(bed_file.nBlocks_rev)
nblocks_set_fwd =     {x for x in nblocks_set_fwd if x>0}
nblocks_set_rev =     {x for x in nblocks_set_rev if x>0}    

# %% 
# make isoform combined by paired reads & filter out isoforms with small BlockSizes
if 'Junction_Name' in bed_file.columns == False:
    bed_file[[col for col in bed_file.columns if 'cigar' not in col ]] = bed_file[[col for col in bed_file.columns if 'cigar' not in col ]].fillna(0).astype(int)
bed_file['isoform'] = bed_file.apply(lambda x: makeisoform(x), axis='columns')
bed_file = bed_file.assign(underscorecount = bed_file.isoform.str.count('_'))

''' 
filter for: 
1. small box sizes
2. many nucleotide deletions
3. StartLeftBlock below 600
'''
if os.path.exists(transcript_path + 'rawsplicedata/') == False:
    os.mkdir(transcript_path + 'rawsplicedata/')
bed_file.to_csv(transcript_path + 'rawsplicedata/'+ chr + '.csv', index=False)
bed_file = bed_file[((bed_file.sumBlockSize_rev > 110) & (bed_file.sumBlockSize_fwd > 110)) & ((bed_file.deletions_fwd < 4) & (bed_file.deletions_rev < 4)) & ((bed_file.nBlocks_fwd <= 3) & (bed_file.nBlocks_rev <= 3)) & (bed_file.StartLeftBlock_rev < 600) & (bed_file.StartLeftBlock_rev < 600) & (bed_file.underscorecount <= 4)] ###### should be added & (bed_file.StartLeftBlock_fwd < 600) (right ?) --- APRIL 29th 2025
if len(bed_file) == 0:
    print('The library sequence ' + str(chr) + ' has no isoforms. Process here is done')
    os.system(f'echo {chr} >> {transcript_path}noisoform.txt')
    exit()
# %%
# count reads
bed_file = bed_file.merge(bed_file.isoform.value_counts(), on='isoform').rename(columns={'count':'nreads_isoform'})
bed_file = bed_file.assign(nreads_chr = len(bed_file))
bed_file = bed_file.assign(chr = chr)

# %% 
# clean
bed_file.drop_duplicates('isoform', inplace=True)
bed_file = bed_file.drop(columns=['chr_' + side for side in ['fwd', 'rev'] if 'chr_' + side in bed_file.columns])
# %% 
# add library information 
bed_file['chr'] = bed_file.chr.astype(int)
bed_file = bed_file.merge(lib[['exonstart in sequence', 'exonend in sequence', 'exonlength','subset', 'gene name', 'rnareads_effective', 'first_ss', 'second_ss', 'changes']],
                          left_on='chr', right_index=True)

bed_file['exonstart in sequence'] = bed_file['exonstart in sequence'] + 273
bed_file['exonend in sequence'] = bed_file['exonend in sequence'] + 273

bed_file = bed_file.merge(lib[['maxent5', 'maxent3']], left_on='chr', right_index=True)
# %%
rename_dict = {1: 'intron1start_o', 2: 'intron1end_o', 3: 'intron2start_o', 4: 'intron2end_o', 5: 'intron3start_o', 6: 'intron3end_o', 7: 'intron4start_o', 8: 'intron4end_o'}
splicecites_df = bed_file.isoform.str.split('_', expand=True).replace('',0).fillna(0).astype(int)
print(str(chr) + ' number of junctions: ' + str(len(splicecites_df.columns) -1 ))

# splicecites_df = splicecites_df.rename(columns={col: rename_dict[col] for col in [1,2,3,4] if col in splicecites_df.columns}).drop(columns=0)
splicecites_df = splicecites_df.rename(columns={col: rename_dict[col] for col in splicecites_df.drop(columns=0).columns}).drop(columns=0)
bed_file = bed_file.join(splicecites_df)

# %%
def maxent(pd_row, s, side, dff):
    if pd_row['intron' + str(s) + side + '_o'] not in dff.columns:
        return 0, 0
    if pd_row['intron' + str(s) + side + '_o'] == dff.columns.min():
        return 0, 0
    nucleotides_locations_around = [pd_row['intron' + str(s) + side + '_o'] + n for n in range(-3,4) if pd_row['intron' + str(s) + side + '_o'] + n in dff.columns]
    if len(nucleotides_locations_around) > 1:
        max_score_local = dff.loc[int(pd_row['chr']), nucleotides_locations_around].max()       # 30 APRL 2025
        max_nuclt_local = dff.loc[int(pd_row['chr']), nucleotides_locations_around].idxmax()    # 30 APRL 2025
        return max_score_local, max_nuclt_local
    elif len(nucleotides_locations_around) == 1:
        return dff.loc[int(chr), nucleotides_locations_around], dff.loc[int(chr), nucleotides_locations_around]
    elif len(nucleotides_locations_around) == 0: 
        return 0, 0
    else: 
        return 0, 0

def makeisoform(pd_row):
    l = [pd_row[col] for col in [f'intron{n}_maxent{f}_nt' for n in range(1,4) for f in [3,5]] if col in bed_file.columns]
    l.sort()
    return str(pd_row.chr) + '_' + '_'.join([str(int(e)) for e in l if e > 0])
# %%

max_underscore_count = max(bed_file.underscorecount.unique())       # 30 APRL 2025
number_of_exons = int(np.where(max_underscore_count == 4, 3, np.where(max_underscore_count == 2, 2, 0)))
number_of_introns = number_of_exons - 1                    # 30 APRL 2025
for intrn_cnt in range(1,number_of_exons  ):               # 30 APRL 2025
    side = 'start'
    # donors
    if 'intron' + str(intrn_cnt) + side + '_o' in bed_file.columns:
        bed_file[[f'intron{intrn_cnt}_maxent5_sc', f'intron{intrn_cnt}_maxent5_nt']] = bed_file.apply(lambda x: maxent(x, s=intrn_cnt, side='start', dff=imaxent5), axis='columns').to_list()
        bed_file[[f'intron{intrn_cnt}_maxent5_sc', f'intron{intrn_cnt}_maxent5_nt']] = bed_file[[f'intron{intrn_cnt}_maxent5_sc', f'intron{intrn_cnt}_maxent5_nt']].astype(float).fillna(0).astype(int)
    side = 'end'
    # acceptors
    if 'intron' + str(intrn_cnt) + side + '_o' in bed_file.columns:

        bed_file[[f'intron{intrn_cnt}_maxent3_sc', f'intron{intrn_cnt}_maxent3_nt']] = bed_file.apply(lambda x: maxent(x, s=intrn_cnt, side='end', dff=imaxent3), axis='columns').to_list()
        bed_file[[f'intron{intrn_cnt}_maxent3_sc', f'intron{intrn_cnt}_maxent3_nt']] = bed_file[[f'intron{intrn_cnt}_maxent3_sc', f'intron{intrn_cnt}_maxent3_nt']].astype(float).fillna(0).astype(int)

bed_file = bed_file.drop(columns='isoform').assign(isoform = bed_file.drop(columns='isoform').apply(lambda x: makeisoform(x), axis='columns')   )
bed_file = bed_file.drop(columns='nreads_isoform').merge(pd.Series(bed_file.groupby('isoform').nreads_isoform.sum(), name='nreads_isoform'), left_on='isoform', right_index=True)
bed_file.drop_duplicates('isoform', inplace=True)
bed_file = bed_file.assign(underscorecount = bed_file.isoform.str.count('_'))

if ('intron1_maxent3_nt' in bed_file.columns) and  ('intron2_maxent5_nt' in bed_file.columns):
    bed_file['intron2_maxent5_nt'] = bed_file.apply(lambda x: np.where( x['intron1_maxent3_nt'] > 1100, 0, x['intron2_maxent5_nt']), axis=1)
if ('intron1_maxent3_nt' in bed_file.columns) and  ('intron2_maxent3_nt' in bed_file.columns):
    bed_file['intron2_maxent3_nt'] = bed_file.apply(lambda x: np.where( x['intron1_maxent3_nt'] > 1100, 0, x['intron2_maxent3_nt']), axis=1)

if 'isoform' in bed_file.columns:
    bed_file.drop(columns='isoform', inplace=True)
if 'underscorecount' in bed_file.columns:
    bed_file.drop(columns='underscorecount', inplace=True)
# bed_file['isoform'] = bed_file[['chr'] + bed_file.filter(like='_nt').columns.to_list()].apply(lambda x: '_'.join([str(int(e)) for e in x.values if e > 0]) , axis='columns') # 9 MAY 2025
bed_file = bed_file.assign(isoform = bed_file.apply(lambda x: makeisoform(x), axis='columns')   )
bed_file = bed_file.assign(underscorecount = bed_file.isoform.str.count('_'))
bed_file = bed_file.drop(columns=['nreads_isoform']).merge(pd.Series(bed_file.groupby('isoform').nreads_isoform.sum(), name = 'nreads_isoform'), left_on='isoform', right_index=True)
bed_file.drop_duplicates('isoform', inplace=True)
bed_file = bed_file.drop(columns=['nreads_chr']).merge(pd.Series(bed_file.groupby('chr').nreads_isoform.sum(), name = 'nreads_chr'), left_on='chr', right_index=True)
bed_file.drop(columns=['exonlength_o', 'nExons'], inplace=True, errors='ignore')

bed_file['nExons'] = bed_file.apply(lambda x: np.where(x.underscorecount==4,3,np.where(x.underscorecount==2,2,1)), axis='columns')
bed_file = bed_file.join(pd.Series(bed_file.apply(lambda x: x.underscorecount -1 , axis='columns'), name='nIntrons').astype('float') ) # 30 APRL 2025
if ('intron2_maxent5_nt' in bed_file.columns) and ('intron1_maxent3_nt' in bed_file.columns):
    bed_file = bed_file.join(pd.Series(bed_file.apply(lambda x: np.where(x.nExons==3,x.intron2_maxent5_nt - x.intron1_maxent3_nt,np.NaN), axis='columns'), name='exonlength_o').astype('float'))
# %%
# find if the splicesites are canonical or not 
const0 =  'atggtgagcaagggcgaggaggataacatggccatcatcaaggagttcatgcgcttcaaggtgcacatggagggctccgtgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgagggcacccagaccgccaagctgaaggtgaccaagggtggccccctgcccttcgcctgggacatcctgtcccctcagttcatgtacggctccaaggcctacgtgaagcaccccgccgacatccccgactacttgaagctgtccttccccgagggcttcaagtgggagcgcgtgatgaacttcgaggacggcggcgtggtgaccgtgacccaggactcctccctgcaggacggcgagttcatctacaaggtgaagctgcgcggcaccaacttcccctccgacggccccgtaatgcagaagaagaccatgggctgggaggcctcctccgagcggatgtaccccgaggacggcgccctgaagggcgagatcaagcagaggctgaagctgaaggacggcggccactacgacgctgaggtcaagaccacctacaaggccaagaagcccgtgcagctgcccggcgcctacaacgtcaacatcaagttggacatcacctcccacaacgaggactacaccatcgtggaacagtacgaacgcgccgagggccgccactccaccggcggcatggacgagctgtacaagccgGACCGt'
# len(const0) = 717
const1 =  'ACTAGTTTACGACGGGTTGGGGATGGCGTGCAGCGCAACCACGAGACGGCCTTCCAAGGTAAGGGGGTTCATTAATCGCCAAGGCCTCACTCCCTTTTTTCCATCTCTCCCCGGACTCACCCGCCAAGGGTGGGTTGGAAACCGAAACGAGTCAGTGTTGAAACGTGTCTCATCCTATTCCTGAAGCCAGAATATTCTGGCCATGAGTCATTGTTTCCGCCCATCTTGATTCTTTTGGAAATGGCAGCTCTTGTTCAAAGACCGGAAAGGGTGGGATGTCAAGACGTC'
# len(const1) = 288
const2 =  'GGCGCGCCtTGGAGTGGAAGTAGAATGAAGGATTTTTTTTAGAGAGGTGGGGATATCTAAAGGTTTTTATGACGCACGGCTGTTTGCAGGCTCTAACTAAAGGACCATTGTTTATTTGATGTTGATTTAAGTAGTGGATCCTTAGAGATAGTGGTATGGCGGTCTTGAATTGTATCAAAAATCTTGGTTTTCTCTAGGCAATTTTTTGTTCCAATTCAGTTGAATACTCTTCAGTGGATTCAAACCATGAAAAAATAAGTCACCAGGGGAGGATAGCTGAAATAATTCCTAAGGCGGTGCCTGTTTTAATGGAGAAGATATGGGGTGGAGCCTGCGTTTTAAACAAACCCAGATCTGATGCAGGATGTACTTAACTACGTTGAGAAAAACTGATCTGCGCAATTGAGGCGTTACTGAAATATTAGGTGGTGGAGATTTGAGAATAAGGGTTTTCGTCTTTTACCTCATGGGAACTCTGGAAGTCCTTTTGTTAGGATAAATCCTAATAAGACCAAGATAGTACTGTAAAATGAAGTTTAATTATCATGGGTCCCCGCTTAAGAAACTGAAGAACTTATTTTCTTTTTTTGCCCCGGGGTGAATAATAATTGGTTTACTATTGCTTTAGGGGGAAACCTTAGATATTTTAATTTACCTTCTCTCTGGATAGTAGTGTTGTAAGAGAGCAGAAACCCATACTTGAAAATGTGCTTTTCTTTTTTGTTTTCTAGGATGGGTTTGTGGAGTTCTTCCATGTAGAGGACCTAGAAGGTGGCATCAGGAATGTGCTGCTGGCTTTTGCAGGTGTTGCTGGAGTAGGAGCTGGTTTGGCATCTAGAGCTCGGACGGGTGCGCTCCATATGgtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaagcggccgcttccctttagtgagggttaatgcttcgagcagacatgataagatacattgatgagtttggacaaaccacaactagaatgcagtgaaaaaaatgctttatttgtgaaatttgtgatgctattgctttatttgtaaccattataagctgcaataaacaagttaacaacaacaattgcattcattttatgtttcaggttcagggggagatgtgggaggttttttaaagcaagtaaaacctctacaaatgtggta'
# len(const2) = 1839

ref = lib.loc[int(chr), 'sequence'][0:30] + const1 + lib.loc[int(chr), 'sequence'][45:] + const2

for s in range(1,number_of_exons):
    if f'intron{s}_maxent5_nt' in bed_file.columns:
        bed_file[f'intron{s}_maxent5_canonical'] = bed_file.apply(lambda x: 'GT' in ref[x[f'intron{s}_maxent5_nt']-2:x[f'intron{s}_maxent5_nt']+2], axis='columns')
    if f'intron{s}_maxent3_nt' in bed_file.columns:
        bed_file[f'intron{s}_maxent3_canonical'] = bed_file.apply(lambda x: 'AG' in ref[x[f'intron{s}_maxent3_nt']-2:x[f'intron{s}_maxent3_nt']+2], axis='columns')
# %%
print('end of process for library sequence: ' + chr)
if os.path.exists(transcript_path + 'csvfiles/') == False:
    os.mkdir(transcript_path + 'csvfiles/')
bed_file.to_csv(transcript_path + 'csvfiles/'+ chr + '.csv', index=False)
# if os.path.exists(transcript_path + 'csvfiles/') == False:
#     os.mkdir(transcript_path + 'csvfiles/')
# bed_file.to_csv('pm425a'+ chr + '.csv', index=False)