'''
For each read:

1. git the location of the barcode: by finding the constant region (I use only the last 10 nt of the constant region because we don't allow mismatches in that area) 
2. extract the barcode using the location
3. using the barcode, check which library variant it corresponds to (bcindexdict[barcode]) 
4. match the library variant to the readID from the FASTQ files.

'''

# %%
# load modules
import pandas as pd
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import pickle
import os
# %%
parser = argparse.ArgumentParser(description='Demultiplex paired reads to library IDs')
parser.add_argument('-r1', '--read1', type=str, help='read1 gzipped fastq path')
parser.add_argument('-r2', '--read2',type=str, help = 'read2 gzipped fastq path')
parser.add_argument('-o', '--outdir',type=str, default= '~/ce', help = 'output path directory (do not include "/" at the end of the string)')

args = parser.parse_args()
# %%
# load files
csvfile = pd.read_csv('41467_2019_12642_MOESM10_ESM.csv.gz', index_col='library index', compression='gzip')
csvfile['barcode'] = csvfile.sequence.apply(lambda x: x[18:18+12])

libbcdict=csvfile.barcode.to_dict() # keys: library variant IDs ; vals: barcode sequence

# %%
# prepare variables 
# count the number of reads in the file
readcount=0
barcodelength = 12    
bcindexdict={barcodeseq:libID for libID, barcodeseq in libbcdict.items()} # keys: barcode sequence ; vals: library variant IDs

readcounts={idx:0 for idx in libbcdict.keys()} # keys: library IDs ; vals: read counts per index
readids={idx:[] for idx in libbcdict.keys()} # keys: library IDs ; vals: list of read IDs

unmappedreadcount= 0
unmappedreadseqs=[]
mappedreadcount=0
mappedreadseqs=[]
differentreadbcs = dict()
ce_bc = 'CGGTATGCGC'
# %%
# open with gzip the gz fastq files
with gzip.open(args.read1, 'rt') as r1, gzip.open (args.read2, 'rt') as r2:
    
    # 1. sort read IDs according to library barcodes. 
    
    # add values to barcode_sort library 
    for title, seq, qual in FastqGeneralIterator(r1): # find the read IDs per variant for read 1
        readcount +=1
        
        # 1.1. find location of last 10 nucleotides of constant region
        
        primerstart= seq.find(ce_bc) # find location of last 10 nucleotides of constant region
        
        # 1.2. if not found
        if primerstart == -1:
            unmappedreadcount+=1
            unmappedreadseqs.append([title, seq, qual])

        else: # 1.3. if found
            
            # seqeunce after the constant region, till the end of the barcode (the whole read barcode)
            readbc = seq[primerstart+10:primerstart+10+barcodelength]
            differentreadbcs[title.split(' ')[0]] = readbc
            mappedreadcount+=1
            mappedreadseqs.append([title, seq, qual])
            
            if readbc in bcindexdict.keys(): # if the barcode is in the barcode keys
                readcounts[bcindexdict[readbc]]+=1 # add one to readcounts to the ID
                readids[bcindexdict[readbc]].append(title.split(' ')[0]) # add the read ID to the dicionary 
                    
    # find the read IDs per variant for read 2                                    
    for title, seq, qual in FastqGeneralIterator(r2):
        primerstart= seq.find(ce_bc) # find location of last 10 nucleotides of constant region
        
        if primerstart == -1:
            unmappedreadcount+=1
            unmappedreadseqs.append([title, seq, qual])
        else:
            readbc=seq[primerstart+10:primerstart+10+barcodelength] # seqeunce after the constant region, till the end of the barcode (the whole read barcode)

            if readbc in bcindexdict.keys(): # if the barcode is in the barcode keys
                readcounts[bcindexdict[readbc]]+=1 # add one to readcounts to the ID
                readids[bcindexdict[readbc]].append(title.split(' ')[0]) # add the read ID to the dicionary (keys are library ids)
            else:
                unmappedreadcount+=1
                unmappedreadseqs.append(seq)

if args.outdir[-1] != '/':
    path_to_output = args.outdir + '/'
else:
    path_to_output = args.outdir 

if not os.path.isdir(path_to_output):
    os.mkdir(path_to_output)

f = open(path_to_output + "readids.pkl","wb")
# write the python object (dict) to pickle file
pickle.dump(readids,f)

f.close()
