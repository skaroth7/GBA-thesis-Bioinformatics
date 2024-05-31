#Create conda container
- conda create -n ONT_amplicon_pipeline

#Activate conda environment
- conda activate ONT_amplicon_pipeline

#requirements 
Install requirements
- conda install bcftools
- conda install vcftools
- conda install porechop
- conda install minimap2
- conda install nanopolish
- conda install whatshap
- conda install snakemake
- conda install qualimap

#first run prepare_fastq_data.py
- python3 prepare_fastq_data.py "/sequenceout_directory/"

#modify the config.json to specify your required sequencing files and number of cores per sample. You can run multiple samples using multiple threads additionally, by specifying the number of cores when running the snakefile with the instructions described below (INT = number of samples to run together).

run the snakemake pipeline
- snakemake "snakefile" --cores INT
