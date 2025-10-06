# Description
Mapping and upstream analysis pipeline for Illumina reads from splicing MPRA experiments. The aim is to detect splicing isoforms in an unbiased manner.
A computational pipeline was developed to describe the entire spectrum of splice site choices from MPRA paired-end RNA sequencing data experiments.  
The pipeline performs barcode-based demultiplexing, splice-aware read mapping, and quantitative splice site analysis across 9,608 library variants simultaneously.
This is done using a computational workflow comprising of 3 stages:  
1. Demultiplexing and barcode assignment using `1.create_bcread_dict.py`  
2. Read mapping and splice site detection using `ce.pipeline.py` and `ce.detect_splicesites.py`, for simultanously processing, this is run using `2.ce.activate_star_pipeline.sh`  
3. Data integration of all the library variants using `3.ce.concat.py`  

### Computational Requirements

The pipeline was executed on a high-performance computing cluster with:
- Operating System: Red Hat Enterprise Linux 9 operating system distribution
- SLURM job scheduling for parallel processing (v23.02.4)
- Array jobs processing up to 9,608 library variants simultaneously
- 8 CPU cores per job
- Temporary directory management for intermediate files


# Setup pipeline

## Download repository, activate conda environment, download csv and fastq files
1. Install the repository
2. Create a new working conda environemnt from the `.yml` file: 
```bash
conda env create -f spl_mpra_map.yml
```
3. Download cassette exon FASTQ files 


# Setup splice strenght score matrices using [maxentpy](https://github.com/kepbod/maxentpy)
Pre-computes MaxEnt splice site strength scores for all possible nucleotide positions across all library variant sequences.
The script imports the maxentpy module from a local directory (`MaxEntScorepy`), which provides functions for calculating splice site scores using the Maximum Entropy algorithm. The module includes `load_matrix5()` and `load_matrix3()` functions for loading donor and acceptor site position weight matrices, and `score5()` and `score3()` functions for calculating scores.

**Output files:**  
1. MaxEnt donor splice site strength score: `ce_istartmaxent5.csv.gz`  
2. MaxEnt acceptor splice site strength score: `ce_iendmaxent3.csv.gz`  


## Instructions
Run the Following:
```bash
python maxent_table.py
```


# stage 1: Barcode Demultiplexing and Read Assignment

## Input Data Processing

Paired-end FASTQ files (read1 and read2) from multiplexed sequencing experiments were processed using a custom Python script (`1.create_bcread_dict.py`) built on BioPython.  

## Barcode Identification and Mapping

Read assignment to library variants employed a two-step barcode identification process:

1. **Constant Region Location**: Each read was scanned for the 10-nucleotide constant region sequence `CGGTATGCGC` using Python's string `find()` method. This sequence represents the terminal portion of the constant region, chosen to minimize false positives while maintaining mapping specificity.

2. **Barcode Extraction**: Upon a constant region detection, the following 12 nucleotides were extracted as the putative barcode sequence.

3. **Library Assignment**: Extracted barcodes were matched against the pre-constructed barcode dictionary using exact string matching. Successfully matched reads were assigned to their corresponding library variants, with read identifiers stored in variant-specific lists.

Read assignment results were serialized using Python's pickle module for efficient processing reads from every library variants simultaneously.

## Instructions

1. Edit the file `1.activate_1.create_bcread_dic.sh`:
    a. Change the path and FASTQ files for `read1` and `read2`
    b. Name an output directory to be created. 
2. Run the file

Or directly run:
```bash
sbatch -J <output_dir_name> -o <output_dir_name>.out -e <output_dir_name>.err --wrap "python 1.create_bcread_dict.py -r1 <read1>.fastq.gz -r2 <read2>.fastq.gz -o <output_dir_name>"
```


# Stage 2: Mapping Reads corresponding to library every library variant

For per variant mapping using STAR mapper. Mapping will be proceed after the following processes are done: 
1. Pair-read FASTQ of a single library variant
2. STAR Reference Genome files for a single library variant


#### Synthetic Reference Generation and Genome indexing with STAR

For each library variant, a synthetic reference sequence was constructed.  
Individual STAR genome indices were generated for each library variant using the following parameters:
- `--genomeSAindexNbases 4`: Optimized for short reference sequences
- `--runThreadN 8`: Parallel processing on 8 CPU cores
- `--runMode genomeGenerate`: Genome generation mode

#### Read Subset Extraction

Library-specific reads were extracted from the original FASTQ files using seqtk with read identifier lists generated during demultiplexing:
```bash
seqtk subseq [input.fastq.gz] [variant.lst] > [variant.fastq]
```

#### Splice-Aware Read Mapping

STAR alignment was performed with parameters optimized for splice junction detection:
- `--outSAMtype BAM Unsorted`: Direct BAM output without coordinate sorting
- `--runThreadN 8`: Multi-threaded processing
- Paired-end mode with R1 and R2 FASTQ inputs

#### Post-Mapping Quality Control

Aligned reads were filtered to retain only properly paired reads using samtools:
```bash
samtools view -bf 0x2 input.bam > output.paired.bam
```

Mapping statistics were collected including:
- Total aligned reads
- Unmapped reads
- Reads containing splice junctions (CIGAR 'N' operations)
- Reads without splice junctions
- Reads with soft-clipping (CIGAR 'S' operations)


## Instructions

1. Edit the file `2.ce.activate_star_pipeline.sh`:
    a. Within the python run command line change:
         i. The path and FASTQ files for `read1` and `read2` (to be identical as the former sh script)
         ii. The path assigned for `--pwrdir` means path workking directory, change it to the output directory.
    b. Change the following parameters:
   ```bash
   --job-name= <change_job_name>
   --output=<change_path>/.%A_%a.output
   --partition=<change_slurm_partition_to_run_1_day>
   ```
2. Run the file

These instructions will run the script `ce.pipeline.py` as a job array of ~9600 jobs (one job per library variant).



# Stage 3: Splice Site Detection and Quantification

#### BAM to BED Conversion

Aligned reads were converted to BED12 format using bedtools with preservation of splice junction information:
```bash
bamToBed -split -bed12 -i input.bam
```
Additional read sequence and CIGAR string information were extracted using SAMtools and merged with BED coordinates to create comprehensive splice junction records.

#### Paired-Read Integration

Since the analysis used paired-end sequencing, forward and reverse read information was integrated by:
1. Separating even-indexed (forward) and odd-indexed (reverse) records
2. Joining paired records based on read identifiers
3. Combining block coordinate information from both reads

Up to two introns per transcript isoform were supported, accommodating sequences with up to three exons.

#### Quality Filtering

Multiple quality control filters were applied:

1. **Read Length Filters**: Both forward and reverse reads required >110 nucleotides of aligned sequence
2. **Deletion Filters**: Reads with >4 nucleotides of deletions (CIGAR 'D' operations) were excluded
3. **Complexity Filters**: Maximum of 3 exons per isoform were allowed
4. **Position Filters**: Reads mapping with start positions >600 nucleotides were excluded to eliminate mis-mapping artifacts

#### MaxEnt Splice Site Scoring

Splice site strength was evaluated using MaxEnt algorithm scores. For each detected splice site, the algorithm searched a 7-nucleotide window (±3 nucleotides downstream and upstread apart) around the detected position to identify the optimal splice site and assign the maximum score within this region.

#### Isoform Quantification

Unique isoforms were defined by their complete set of splice junction coordinates. Read counts were aggregated per isoform, and the following metrics were calculated:
- Reads per isoform
- Total reads per library variant  
- Relative isoform abundance
- Number of exons per isoform
- Exon length measurements for multi-exonic isoforms


## Instructions
No running instructions since the script `ce.pipeline.py` will run the script `ce.detect_splicesites.py` automatically.


# Stage 4: Data Integration and Output

#### Result Compilation

Individual library variant results were concatenated using pandas DataFrame operations. The final dataset included:
- Library variant identifiers and metadata
- Complete splice junction coordinate sets
- MaxEnt scores for all detected splice sites
- Canonical/non-canonical classifications
- Read count quantifications
- Isoform structural annotations

#### Output Format

Results were exported as comma-separated value (CSV) files containing one row per unique isoform with columns for:
- Library identification and metadata
- Splice site coordinates and scores
- Quantification metrics
- Quality control flags


## Instructions  

While running the command:
```bash
python 3.ce.concat.py -p <path_of_inputs_directory> -n <prefix_for_output_filename>
```

the `-p` flag input string should be identical to the `<output_dir_name>` when the script `1.create_bcread_dict.py` was run.

## Output file  
On the same directory the scripts are found, the integrated splice isoform descriptions are shown in a CSV file.

