Create conda container
- conda create -n ONT_amplicon_pipeline

Activate conda environment
- conda activate ONT_amplicon_pipeline

Requirements:
Install requirements
- conda install bcftools
- conda install vcftools
- conda install porechop
- conda install minimap2
- conda install nanopolish
- conda install whatshap
- conda install snakemake
- conda install qualimap

First run prepare_fastq_data.py
- place the prepare_fastq_data.py file into the base directory containing all of the sequencing output folders
- python3 prepare_fastq_data.py "/sequenceout_directory/"

Modify the config.json to specify your required sequencing files and number of cores per sample. Additionally, you can run multiple samples using multiple threads by specifying the number of cores when running the snakefile with the instructions described below (INT = number of samples to run together).

Run the snakemake pipeline
- snakemake "snakefile" --cores INT
