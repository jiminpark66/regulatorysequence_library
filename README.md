# Regulatory Sequence Library

Code used to quantify transcriptional activities values from a multiplexed reporter assay referenced in manuscript X. Currently, a simple demo is provided and the repository will be updated later for more general applications. 

Manuscript reference X.
Raw data can be accessed at X.

# Dependencies and required softwares
bbMerge
pandas  
numpy  
xlrd  

# Demo
For the purposes of the demo, download a single sample from the repository accession. A single sample should consist of four raw sequencing files; DNA read1, DNA read2, RNA read1, RNA read2.  
Next, download demo directory, contains empty output file directories, library reference file, etc.  
Download these files:  

# Prep sequencing files
- Create a new directory "raw" and populate with raw sequencing files. Also create "unmerged" and "merged" files for output.
- Pair-end merge sequencing files using SeqPrep using default settings.
- To easily prep multiple samples at the same time, use the mk_mergefile.py provided. Just update the mk_mergefile script with the appropriate parameters
- Run the resulting bashscript merge.sh

# Parse merged sequencing files
- Create these following folders [0_bccounts
1_lowq
2_missingadapter
3_frag
4_badbc
05_bc_counts
5_goodbc_badalign
6_goodbc_goodalign
7_goodbc_perfectalign
8_goodbc_goodalign_bccounts
9_goodbc_perfectalign_bccounts
10_log_files]
- Library reference file is provided as X (3_lib_final.xlsx)
- Run the python file parse.py (or the jupyter notebook of the same name) which iterates through merged DNA and RNA files, checks sequence read qualities and writes counts of each regulatory sequence elements


# Quantify transcriptional activities
- Run the python file analyze.py (or the jupyter notebook of the same name)
- Uses counts of each elements from each sample to compute relative transcriptional activity (tx values) 
- Final output is written to masterdf.csv
