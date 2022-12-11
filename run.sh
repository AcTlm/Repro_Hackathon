conda activate snakemake
snakemake -s wkf_download-data_8SRR_index.txt --cores 8 --use-singularity
snakemake -s wkf_map-analyze_8SRR.txt --cores 2 --use-singularity