# description
A computational pipeline was developed to describe the entire spectrum of splice site choices from MPRA paired-end RNA sequencing data experiments.  
The pipeline performs barcode-based demultiplexing, splice-aware read mapping, and quantitative splice site analysis across 9,608 library variants simultaneously.

# step 0: install scripts, anaconda environemnt and csv files
1. intall the repository
2. create a new working conda environemnt from the `.yml` file: 
```conda env create -f spl_mpra_map.yml```
3. install CSV files used in the workflow throught <<link>>  :
   a. **41467_2019_12642_MOESM10_ESM.csv.gz** (library variants metadata)
   b. **ce_istartmaxent5.csv.gz** (library variante donor maxent scores)
   c. **ce_iendmaxent3.csv.gz** (library variante acceptor maxent scores)

``` wget <link>41467_2019_12642_MOESM10_ESM.csv.gz  . ```
``` wget <link>ce_istartmaxent5.csv.gz  . ```
``` wget <link>ce_iendmaxent3.csv.gz  . ```
