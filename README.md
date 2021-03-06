# Regulatory sequence library

Transcriptional levels of regulatory sequences can be measured in high-throughput using Massively Parallel Reporter Assays (MPRA). Materials and code used to process and analyze raw sequencing data from regulatory sequence MPRA experiments are provided here. This pipeline was utilized in the following manuscripts: "Systematic dissection of σ70 sequence diversity and function in bacteria" and "High-throughput transcriptional measurements of regulatory sequences derived from bacterial biosynthetic gene clusters".

Manuscripts can be accessed here: [add link] [add link]  
Raw data for the manuscript can be accessed here: E-MTAB-9111 [σ70] [add link], E-MTAB-9223 [BGC] [add link]


# Dependencies
- Python 3.X
    - biopython
    - pandas
    - numpy 
    - xlrd 
- bbMerge

# 1. Raw file processing and other prerequisites
A single sample should consist of four raw sequencing files: DNA read1, DNA read2, RNA read1, and RNA read2. Download the files in the demo directory. To demo this pipeline, download these samples from the E-MTAB-9111 repository `176_DNA_1_S9_R1.fastq.gz, 176_DNA_1_S9_R2.fastq.gz, 176_RNA_1_1_S21_R1.fastq.gz, 176_RNA_1_1_S21_R2.fastq.gz`.

Create the following output directories:
```
0_bccounts/
1_lowq/
2_missingadapter/
3_frag/
4_badbc/
05_bc_counts/
5_goodbc_badalign/
6_goodbc_goodalign/
7_goodbc_perfectalign/
8_goodbc_goodalign_bccounts/
9_goodbc_perfectalign_bccounts/
10_log_files/
merged/
raw/
unmerged/
```

Library sequences for the demo is provided from `3_lib_final.csv` file. Using that file as a reference template, other custom libraries can be analyzed.
Move the raw sequencing files (fastq.gz files) to the raw directory. In the `mk_mergefile.py` file, edit the di (directory) variable as appropriate and run the python script to generate a bash script, `merge.sh` that will merge all the files in the `raw/` directory.  


# 2. Parse sequence counts
Run the python script `parse.py` (alternatively, jupyter notebook of the same name is also provided). This iterates through merged DNA and RNA files, performs necessary quality and length threshold checks, and counts abundances of each regulatory sequence in the library. Variables in the `parse.py` script can be modified to adapt the script to other libraries of varying adapter sequences, lengths, etc. 

# 3. Quantify transcriptional activities
Run the python script `analyze.py` (alternatively, jupyter notebook of the same name is also provided). This uses the counts data to compute relative transcriptional activity (Tx values) for each regulatory sequence. Final output is written to `masterdf.csv`
