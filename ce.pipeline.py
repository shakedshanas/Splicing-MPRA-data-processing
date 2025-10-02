''' 


'''
# %%
import pickle
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import os
import time 
from sys import exit
import pandas as pd
# %%
parser = argparse.ArgumentParser(description='run pipeline for each chromosome', epilog = 'I wish you all an enjoyable mapping & processing')
parser.add_argument('-r1', '--read1', type=str, help='read 1')
parser.add_argument('-r2', '--read2', type=str, help='read 2')
parser.add_argument('--libseq', type=str, help='index')
parser.add_argument('-p', '--pwrdir', type=str, help='path to put all the putouts')

args = parser.parse_args()
start = time.time() # start timer
# %%
if args.pwrdir[-1] != '/':
    path_to_output = args.pwrdir + '/'
else:
    path_to_output = args.pwrdir 

p = open(path_to_output + 'readids.pkl', 'rb')
readid = pickle.load(p) # dict keys: library variant ID ; dict values: readIDs
csv_file = pd.read_csv('41467_2019_12642_MOESM10_ESM.csv.gz', compression='gzip')
# %%

chr = int(csv_file.loc[int(args.libseq), 'library index'])
print('>>> library sequence: ' + str(chr) + ' <<<')
if len(readid[chr]) == 0:
    print(f'chromosome {chr} has no reads in this experiment.')
    os.system(f'echo {chr} >> {path_to_output}lost_sequences.txt')
    exit()

print(f'make read list for {chr}')
file = open(path_to_output + f'{chr}.lst','w')
for id in readid[chr]:
	file.write(id+"\n")
file.close()

'''
to create the fastq files out of the list I extraced, use the commands:
seqtk subseq ir_test_read1.fastq.gz 34078.lst > 34078.r1.lst
seqtk subseq ir_test_read2.fastq.gz 34078.lst > 34078.r2.lst
'''
print('make fastq')
while True: 
    if os.path.exists(path_to_output + f'{chr}.lst')==True:
        break
os.system(f'seqtk subseq {args.read1} {path_to_output}{chr}.lst > {path_to_output}{chr}.r1.fastq ; seqtk subseq {args.read2} {path_to_output}{chr}.lst > {path_to_output}{chr}.r2.fastq')


'''create fasta file for the library variant'''
### constant reagions
tot_seq = 'atggtgagcaagggcgaggaggataacatggccatcatcaaggagttcatgcgcttcaaggtgcacatggagggctccgtgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgagggcacccagaccgccaagctgaaggtgaccaagggtggccccctgcccttcgcctgggacatcctgtcccctcagttcatgtacggctccaaggcctacgtgaagcaccccgccgacatccccgactacttgaagctgtccttccccgagggcttcaagtgggagcgcgtgatgaacttcgaggacggcggcgtggtgaccgtgacccaggactcctccctgcaggacggcgagttcatctacaaggtgaagctgcgcggcaccaacttcccctccgacggccccgtaatgcagaagaagaccatgggctgggaggcctcctccgagcggatgtaccccgaggacggcgccctgaagggcgagatcaagcagaggctgaagctgaaggacggcggccactacgacgctgaggtcaagaccacctacaaggccaagaagcccgtgcagctgcccggcgcctacaacgtcaacatcaagttggacatcacctcccacaacgaggactacaccatcgtggaacagtacgaacgcgccgagggccgccactccaccggcggcatggacgagctgtacaagccgGACCGtNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACTAGTTTACGACGGGTTGGGGATGGCGTGCAGCGCAACCACGAGACGGCCTTCCAAGGTAAGGGGGTTCATTAATCGCCAAGGCCTCACTCCCTTTTTTCCATCTCTCCCCGGACTCACCCGCCAAGGGTGGGTTGGAAACCGAAACGAGTCAGTGTTGAAACGTGTCTCATCCTATTCCTGAAGCCAGAATATTCTGGCCATGAGTCATTGTTTCCGCCCATCTTGATTCTTTTGGAAATGGCAGCTCTTGTTCAAAGACCGGAAAGGGTGGGATGTCAAGACGTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGCGCGCCtTGGAGTGGAAGTAGAATGAAGGATTTTTTTTAGAGAGGTGGGGATATCTAAAGGTTTTTATGACGCACGGCTGTTTGCAGGCTCTAACTAAAGGACCATTGTTTATTTGATGTTGATTTAAGTAGTGGATCCTTAGAGATAGTGGTATGGCGGTCTTGAATTGTATCAAAAATCTTGGTTTTCTCTAGGCAATTTTTTGTTCCAATTCAGTTGAATACTCTTCAGTGGATTCAAACCATGAAAAAATAAGTCACCAGGGGAGGATAGCTGAAATAATTCCTAAGGCGGTGCCTGTTTTAATGGAGAAGATATGGGGTGGAGCCTGCGTTTTAAACAAACCCAGATCTGATGCAGGATGTACTTAACTACGTTGAGAAAAACTGATCTGCGCAATTGAGGCGTTACTGAAATATTAGGTGGTGGAGATTTGAGAATAAGGGTTTTCGTCTTTTACCTCATGGGAACTCTGGAAGTCCTTTTGTTAGGATAAATCCTAATAAGACCAAGATAGTACTGTAAAATGAAGTTTAATTATCATGGGTCCCCGCTTAAGAAACTGAAGAACTTATTTTCTTTTTTTGCCCCGGGGTGAATAATAATTGGTTTACTATTGCTTTAGGGGGAAACCTTAGATATTTTAATTTACCTTCTCTCTGGATAGTAGTGTTGTAAGAGAGCAGAAACCCATACTTGAAAATGTGCTTTTCTTTTTTGTTTTCTAGGATGGGTTTGTGGAGTTCTTCCATGTAGAGGACCTAGAAGGTGGCATCAGGAATGTGCTGCTGGCTTTTGCAGGTGTTGCTGGAGTAGGAGCTGGTTTGGCATCTAGAGCTCGGACGGGTGCGCTCCATATGgtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaagcggccgcttccctttagtgagggttaatgcttcgagcagacatgataagatacattgatgagtttggacaaaccacaactagaatgcagtgaaaaaaatgctttatttgtgaaatttgtgatgctattgctttatttgtaaccattataagctgcaataaacaagttaacaacaacaattgcattcattttatgtttcaggttcagggggagatgtgggaggttttttaaagcaagtaaaacctctacaaatgtggta'
# len(tot_seq) = 3038
const0 =  'atggtgagcaagggcgaggaggataacatggccatcatcaaggagttcatgcgcttcaaggtgcacatggagggctccgtgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgagggcacccagaccgccaagctgaaggtgaccaagggtggccccctgcccttcgcctgggacatcctgtcccctcagttcatgtacggctccaaggcctacgtgaagcaccccgccgacatccccgactacttgaagctgtccttccccgagggcttcaagtgggagcgcgtgatgaacttcgaggacggcggcgtggtgaccgtgacccaggactcctccctgcaggacggcgagttcatctacaaggtgaagctgcgcggcaccaacttcccctccgacggccccgtaatgcagaagaagaccatgggctgggaggcctcctccgagcggatgtaccccgaggacggcgccctgaagggcgagatcaagcagaggctgaagctgaaggacggcggccactacgacgctgaggtcaagaccacctacaaggccaagaagcccgtgcagctgcccggcgcctacaacgtcaacatcaagttggacatcacctcccacaacgaggactacaccatcgtggaacagtacgaacgcgccgagggccgccactccaccggcggcatggacgagctgtacaagccgGACCGt'
# len(const0) = 717
const1 =  'ACTAGTTTACGACGGGTTGGGGATGGCGTGCAGCGCAACCACGAGACGGCCTTCCAAGGTAAGGGGGTTCATTAATCGCCAAGGCCTCACTCCCTTTTTTCCATCTCTCCCCGGACTCACCCGCCAAGGGTGGGTTGGAAACCGAAACGAGTCAGTGTTGAAACGTGTCTCATCCTATTCCTGAAGCCAGAATATTCTGGCCATGAGTCATTGTTTCCGCCCATCTTGATTCTTTTGGAAATGGCAGCTCTTGTTCAAAGACCGGAAAGGGTGGGATGTCAAGACGTC'
# len(const1) = 288
const2 =  'GGCGCGCCtTGGAGTGGAAGTAGAATGAAGGATTTTTTTTAGAGAGGTGGGGATATCTAAAGGTTTTTATGACGCACGGCTGTTTGCAGGCTCTAACTAAAGGACCATTGTTTATTTGATGTTGATTTAAGTAGTGGATCCTTAGAGATAGTGGTATGGCGGTCTTGAATTGTATCAAAAATCTTGGTTTTCTCTAGGCAATTTTTTGTTCCAATTCAGTTGAATACTCTTCAGTGGATTCAAACCATGAAAAAATAAGTCACCAGGGGAGGATAGCTGAAATAATTCCTAAGGCGGTGCCTGTTTTAATGGAGAAGATATGGGGTGGAGCCTGCGTTTTAAACAAACCCAGATCTGATGCAGGATGTACTTAACTACGTTGAGAAAAACTGATCTGCGCAATTGAGGCGTTACTGAAATATTAGGTGGTGGAGATTTGAGAATAAGGGTTTTCGTCTTTTACCTCATGGGAACTCTGGAAGTCCTTTTGTTAGGATAAATCCTAATAAGACCAAGATAGTACTGTAAAATGAAGTTTAATTATCATGGGTCCCCGCTTAAGAAACTGAAGAACTTATTTTCTTTTTTTGCCCCGGGGTGAATAATAATTGGTTTACTATTGCTTTAGGGGGAAACCTTAGATATTTTAATTTACCTTCTCTCTGGATAGTAGTGTTGTAAGAGAGCAGAAACCCATACTTGAAAATGTGCTTTTCTTTTTTGTTTTCTAGGATGGGTTTGTGGAGTTCTTCCATGTAGAGGACCTAGAAGGTGGCATCAGGAATGTGCTGCTGGCTTTTGaataaCGCGCCTTAGCTCGGACGGGTGCGCTC'
# len(const2) = 1839
# %%
csv_file.set_index('library index', inplace=True)
print('make fasta')
with open(path_to_output + f"{chr}.fasta", "w") as output_handle:
    record = SeqIO.SeqRecord(Seq(csv_file.loc[chr, 'sequence'][0:30] + const1 + csv_file.loc[chr, 'sequence'][45:] + const2),
                                id=str(chr),
                                name='MPRA_experiment',
                                description='cassette_exons')
    SeqIO.write(record, output_handle, "fasta")        
# %%
print('make genome')
if os.path.exists(path_to_output + f'{chr}.GenomeTempDir')==True:
	os.rmdir(path_to_output + f'{chr}.GenomeTempDir')

while True:
    ''' if the files exists, continue to STAR generate genome '''
    if os.path.exists(path_to_output + f'{chr}.fasta')==True:
        break
os.system(f'STAR --runThreadN 8 --runMode genomeGenerate --genomeFastaFiles {path_to_output}{chr}.fasta --genomeDir {path_to_output}{chr}.genome --genomeSAindexNbases 4 --outTmpDir {path_to_output}{chr}.GenomeTempDir')
# %% 
# map the files to the chr
print('map reads')
while True:
    if os.path.exists(path_to_output + f'{chr}.genome/Genome') & os.path.exists(path_to_output + f'{chr}.genome/SA') & os.path.exists(path_to_output + f'{chr}.genome/SAindex') & os.path.exists(path_to_output + f'{chr}.genome/chrLength.txt') & os.path.exists(path_to_output + f'{chr}.genome/chrName.txt') & os.path.exists(path_to_output + f'{chr}.genome/genomeParameters.txt') & os.path.exists(path_to_output + f'{chr}.r1.fastq') & os.path.exists(path_to_output + f'{chr}.r2.fastq'):
        break

if os.path.exists(path_to_output + f'{chr}.MappingTmpDir')==True:
    os.rmdir(path_to_output + f'{chr}.MappingTmpDir')

os.system(f'STAR --runThreadN 8 --genomeDir {path_to_output}{chr}.genome --readFilesIn {path_to_output}{chr}.r1.fastq {path_to_output}{chr}.r2.fastq --outFileNamePrefix {path_to_output}{chr}. --outSAMtype BAM Unsorted --outTmpDir {path_to_output}{chr}.MappingTmpDir')
if os.path.exists(path_to_output + 'otheroutputs/') == False:
    os.mkdir(path_to_output + 'otheroutputs/')

# clean 
os.system(f"mv {path_to_output}{chr}.lst {path_to_output}otheroutputs/")

while True:
    if os.path.exists(path_to_output + f'{chr}.Aligned.out.bam') & os.path.exists(path_to_output + f'{chr}.Log.final.out') & os.path.exists(path_to_output + f'{chr}.Log.out') & os.path.exists(path_to_output + f'{chr}.Log.progress.out')  & os.path.exists(path_to_output + f'{chr}.SJ.out.tab'):
        break

print('pair bam file')
while True:
    if os.path.exists(path_to_output + f'{chr}.Aligned.out.bam') ==True:
        break
if os.path.exists(path_to_output + '/bamfiles/') == False:
    os.mkdir(path_to_output + '/bamfiles/')
os.system(f'samtools view -bf 0x2 {path_to_output}{chr}.Aligned.out.bam > {path_to_output}bamfiles/{chr}.paired.bam')
end = time.time()
os.system(f'echo "chr {chr} took {(end - start)/60} minutes to process bam file"')
os.system(f"echo -e {chr},$(samtools view -F 4 -c {path_to_output}bamfiles/{chr}.paired.bam),$(samtools view -f 4 -c {path_to_output}bamfiles/{chr}.paired.bam),$(samtools view -F 4 {path_to_output}bamfiles/{chr}.paired.bam | awk '$6 ~ /N/' | wc -l),$(samtools view -F 4 {path_to_output}bamfiles/{chr}.paired.bam | awk '$6 !~ /N/' | wc -l),$(samtools view -F 4 {path_to_output}bamfiles/{chr}.paired.bam | awk '$6 ~ /S/' | wc -l),$(samtools view -F 4 {path_to_output}bamfiles/{chr}.paired.bam | awk '$6 !~ /S/' | wc -l),demux_star >> {path_to_output}mapsummary.csv")

# clean
os.system(f"mv {path_to_output}{chr}.Aligned.out.bam {path_to_output}otheroutputs/")
os.system(f"rm {path_to_output}{chr}.fasta {path_to_output}{chr}.Log.*.out {path_to_output}{chr}.tab {path_to_output}{chr}.*.fastq") # clean
os.system(f"rm {path_to_output}{chr}.SJ.out.tab") # clean
os.system(f"rm -r {path_to_output}{chr}.genome") # clean

print('make bed file')
if os.path.exists(path_to_output + '/bedfiles') == False:
    os.mkdir(path_to_output + '/bedfiles/')
while True:
    if os.path.exists(path_to_output + f'bamfiles/{chr}.paired.bam') ==True:
        break
print('make bed file part 1')
os.system(f'bamToBed -split -bed12 -i {path_to_output}bamfiles/{chr}.paired.bam | cut -f1,2,3,4,6,10,11,12 > {path_to_output}bedfiles/{chr}.p1.bed')
print('make file for cigar string and sequence')
os.system(f'samtools view {path_to_output}bamfiles/{chr}.paired.bam | cut -f 6,10 > {path_to_output}bedfiles/{chr}.p2.bed')


print('combine the two together')
while True:
    if os.path.exists(path_to_output + f'bedfiles/{chr}.p1.bed') & os.path.exists(path_to_output + f'bedfiles/{chr}.p2.bed'):
        break
os.system(f"paste -d '	' {path_to_output}bedfiles/{chr}.p1.bed {path_to_output}bedfiles/{chr}.p2.bed > {path_to_output}bedfiles/supbedfile.{chr}.bed")
os.system(f"rm {path_to_output}{chr}.tab {path_to_output}{chr}.Log.out") # clean

print('process bed file')
while True:
    if os.path.exists(path_to_output + f'bedfiles/supbedfile.{chr}.bed') == True:
        break
os.system(f'python ce.detect_splicesites.py --chr {chr} --pwrdir ' + args.pwrdir)