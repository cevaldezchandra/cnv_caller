# cnv_caller

This R script takes the outputs from the program EXCAVATOR (Magi, A., et al., EXCAVATOR: detecting copy number variants from whole-exome sequencing data. Genome Biol, 2013. 14(10): p. R120) and produces a gene level amplification summary for the genes of interest. 

# Setting input paths for EXCAVATOR calculations

### Input should look like:
Rscript CNV_v7.R /path/to/sample.bam targetname /path/to/output

### Takes in 3 commandline prompts: 
1. Path to .bam file
2. Target name used to name target initialization folder and output folder
3. Path to desired output directory - WITHOUT folder name, the script will
   generate individual folders using the targetname

### Input files and data needed (bullets 2-4 are found in main EXCAVATOR folder)
1. Sample.RCNorm_txt (output from EXCAVATOR)
2. 51_total_exons.csv (51 unpaired samples for pooled background)
3. gene_names.txt
4. exon_names.txt
