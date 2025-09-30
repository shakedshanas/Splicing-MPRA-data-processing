# Description
A computational pipeline was developed to describe the entire spectrum of splice site choices from MPRA paired-end RNA sequencing data experiments.  
The pipeline performs barcode-based demultiplexing, splice-aware read mapping, and quantitative splice site analysis across 9,608 library variants simultaneously.
This is done using a computational workflow comprising of 3 stages:  
1. Demultiplexing and barcode assignment using `1.create_bcread_dict.py`  
2. Read mapping and splice site detection using `ce.pipeline.py` and `ce.detect_splicesites.py`, for simultanously processing, this is run using `2.ce.activate_star_pipeline.sh`  
3. Data integration of all the library variants using `3.ce.concat.py`  

### Computational Requirements

The pipeline was executed on a high-performance computing cluster with:
- SLURM job scheduling for parallel processing
- Array jobs processing up to 9,608 library variants simultaneously
- 8 CPU cores per job
- Temporary directory management for intermediate files


# Download repository, activate conda environment, download csv and fastq files
1. Install the repository
2. Create a new working conda environemnt from the `.yml` file: 
```bash
conda env create -f spl_mpra_map.yml
```
3. download CSV files used in the workflow throught <<link>>  :  
   a. **41467_2019_12642_MOESM10_ESM.csv.gz** (library variants metadata)  
   b. **ce_istartmaxent5.csv.gz** (library variante donor maxent scores)  
   c. **ce_iendmaxent3.csv.gz** (library variante acceptor maxent scores)  
```bash
wget <link>41467_2019_12642_MOESM10_ESM.csv.gz  .
wget <link>ce_istartmaxent5.csv.gz  .
wget <link>ce_iendmaxent3.csv.gz  .
```
4. Download cassette exon FASTQ files 




# stage 1: Barcode Demultiplexing and Read Assignment

## Input Data Processing

Paired-end FASTQ files (R1 and R2) from multiplexed sequencing experiments were processed using a custom Python script (`1.create_bcread_dict.py`) built on BioPython.  
The pipeline accepts gzip-compressed FASTQ files and processes them iteratively to minimize memory usage.

## Library Reference Construction

Library variant information was loaded from a CSV reference file containing:
- Library index identifiers  
- Complete sequence constructs  
- 12-nucleotide barcode sequences extracted from positions 18-30 of each construct  
- Associated metadata including gene names, subset classifications, and expected splice sites  

## Barcode Identification and Mapping

Read assignment to library variants employed a two-step barcode identification process:

1. **Constant Region Location**: Each read was scanned for the 10-nucleotide constant region sequence `CGGTATGCGC` using Python's string `find()` method. This sequence represents the terminal portion of the constant region, chosen to minimize false positives while maintaining mapping specificity.

2. **Barcode Extraction**: Upon constant region detection, the subsequent 12 nucleotides were extracted as the putative barcode sequence. The barcode coordinates were calculated as:
   ```python
   barcode_start = constant_region_position + 10
   barcode_end = barcode_start + 12
   ```

3. **Library Assignment**: Extracted barcodes were matched against the pre-constructed barcode dictionary using exact string matching. Successfully matched reads were assigned to their corresponding library variants, with read identifiers stored in variant-specific lists.

#### Quality Control and Statistics

The demultiplexing process tracked several quality metrics:
- Total read count processed
- Successfully mapped reads per library variant
- Unmapped reads due to absent constant regions
- Unmapped reads with unrecognized barcodes

Read assignment results were serialized using Python's pickle module for efficient downstream processing.

## Stage 2: Reference Construction and Read Mapping

#### Synthetic Reference Generation

For each library variant, a synthetic reference sequence was constructed by concatenating:

1. **Variable Region 1**: Nucleotides 0-30 from the library sequence
2. **Constant Region 1**: 288-nucleotide linker sequence containing regulatory elements
3. **Variable Region 2**: Nucleotides 45 to end of library sequence  
4. **Constant Region 2**: 1,839-nucleotide 3' region containing reporter genes and regulatory sequences

Two versions of Constant Region 2 were maintained depending on whether CASLIB primers were used in library construction, differing in a 25-nucleotide sequence modification.

#### STAR Genome Indexing

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

## Stage 3: Splice Site Detection and Quantification

#### BAM to BED Conversion

Aligned reads were converted to BED12 format using bedtools with preservation of splice junction information:
```bash
bamToBed -split -bed12 -i input.bam
```

Additional read sequence and CIGAR string information were extracted using samtools and merged with BED coordinates to create comprehensive splice junction records.

#### Paired-Read Integration

Since the analysis used paired-end sequencing, forward and reverse read information was integrated by:
1. Separating even-indexed (forward) and odd-indexed (reverse) records
2. Joining paired records based on read identifiers
3. Combining block coordinate information from both reads

#### Splice Junction Coordinate Extraction

For multi-block alignments (indicating splice junctions), intron coordinates were calculated using block size and start position information:

- **Intron Start Coordinates**: `StartLeftBlock + BlockSize[i]`
- **Intron End Coordinates**: `StartLeftBlock + BlockStart[i+1]`

Up to two introns per transcript isoform were supported, accommodating sequences with up to three exons.

#### Quality Filtering

Multiple quality control filters were applied:

1. **Read Length Filters**: Both forward and reverse reads required >110 nucleotides of aligned sequence
2. **Deletion Filters**: Reads with >4 nucleotides of deletions (CIGAR 'D' operations) were excluded
3. **Complexity Filters**: Maximum of 3 exons per isoform were allowed
4. **Position Filters**: Reads mapping with start positions >600 nucleotides were excluded to eliminate mis-mapping artifacts
5. **Isoform Complexity**: Maximum of 4 splice junctions per unique isoform

#### MaxEnt Splice Site Scoring

Splice site strength was evaluated using MaxEnt algorithm scores:

1. **Donor Sites (5' splice sites)**: Scored using position weight matrices for GT dinucleotides and surrounding sequence context
2. **Acceptor Sites (3' splice sites)**: Scored using matrices for AG dinucleotides and branch point regions

For each detected splice site, the algorithm searched a 7-nucleotide window (±3 nucleotides) around the detected position to identify the optimal splice site and assign the maximum score within this region.

#### Canonical Splice Site Classification

Splice sites were classified as canonical or non-canonical based on dinucleotide composition:
- **Canonical donor sites**: GT dinucleotide at positions +1/+2 of introns
- **Canonical acceptor sites**: AG dinucleotide at positions -2/-1 of introns

#### Isoform Quantification

Unique isoforms were defined by their complete set of splice junction coordinates. Read counts were aggregated per isoform, and the following metrics were calculated:
- Reads per isoform
- Total reads per library variant  
- Relative isoform abundance
- Number of exons per isoform
- Exon length measurements for multi-exonic isoforms

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

