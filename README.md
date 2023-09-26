#Runs on the Conda package, so first you need to install conda.

#On conda then create a new Enviornment:
#conda create -n new_enviornment

#then activate this enviornment:
#conda activate new_enviornment

#The follow packages are required:
#minimap2
#samtools
#qualimap
#samtools
#nanopolish- nanopolish needs to be installed manually. Then Specify directory in the config file
#whatshap
#bcftools
#snakemake
#conda install -c conda-forge bzip2
#conda install -c anaconda hdf5

#Use
#conda install -n new_envionrnment [package_name] (conda install -c bioconda minimap2 may be needed))

#Place the Snakefile in the Enviornment along with the config.json file then:


modify the config file to specificy the base directory, the reference sequence, the fast5 file, the sequencing summary file and the nanopolish reference name and coordinates: 

E.G
{
    "data": "/WORKSPACE/Oscar/snake_make_workflows/example",
    "reference": "/WORKSPACE/Oscar/snake_make_workflows/example/reference_sequence/GBA.fasta",
    "fast5": "/WORKSPACE/Oscar/snake_make_workflows/example/fast5_pass",
    "nanopolish_reference": GBA_Reference:0-6945
    "sequencing_summary": "/WORKSPACE/Oscar/snake_make_workflows/example/sequencing_summary.txt"
}


now enter in the console:
snakemake --cores [int]
