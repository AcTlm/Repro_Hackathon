git clone https://github.com/AcTlm/Repro_Hackathon 

apt-get install -y wget ###  telecharger wget et install pour tout le workflow
apt-get install  -y gunzip 

samples = ["SRR628582","SRR628583",...]

# Fichiers en sortie du pipeline 
rule all:
    input:
        #### a ecrire 
	
# Creer une liste des SRR pour tout le worlflow 
rule all:
	input:
		expand("fastqc/{sample}_fastqc.html", sample=samples) ###
		default_target: True ### cette règle ets la 

### telechargement des transcriptomes des individus malades
rule DownloadSRR_fastq: 
    input: 
       
    output: "{sample}.fastq"
        
    singularity:
        "docker://AcTlm/docker-sratools"     ????

    shell: 
	docker build . -t dockersra -f  Dockerfile.sra-tools
	docker run -it -v $PWD/ data:rw ncbi/sra-tools prefetch -O /data {samples}"
	docker run -it -v $PWD/ data:rw ncbi/sra-tools fasterq-dump -O /data /data/{samples}/{samples}.sra"
# verifier si ils sont dezipé avec prefetch et fasterq dump 

### verifier la qualité de ces sequences 
rule quality_analysis:
	input: "{samples}.fastq"
	
	output: 
		"{samples}_fastqc.html"               ### "SRR628582_fastqc.zip"
	
	singularity:
        "docker://AcTlm/Dockerfile_fastqc"
	
	shell:
	docker build. -t dockerfastqc -f  Dockerfile.fastqc
	docker run -it -v $PWD/ data 
	"fastqc {input} > {output} "

### télécharger le génome humain
rule Download_ref_human:
    input:
       
    output:
		ref_human.fa

    singularity:
        
    shell: 
	wget -O chromosome.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz
	gunzip -c chromosome.fa.gz > {output}
	rm chromosome.fa.gz


rule DownloadGFF:
    input:
       
    output: 
	Homo_sapiens.GRCh38.101.chr.gtf

    singularity:
        
    shell: 
        wget -O $PWD/Homo_sapiens.GRCh38.101.chr.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
	  gunzip -c Homo_sapiens.GRCh38.101.chr.gtf.gz > {output}
	  rm Homo_sapiens.GRCh38.101.chr.gtf.gz

# Index du génome de référence - génome humain avant l'alignement des reads (séquences)
rule index_REFhumain:
    input: ref_human.fa
       
    output: saindex.sra ## extension???

    shell:
	docker build. -t dockerSTAR -f  Dockerfile.STAR
	docker run -it -v $PWD/data:rw /STAR --runThreadN <nb cpus> --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ref_human.fa

rule mapping:
    input: 
	sample}.fastq"
        saindex.sratools
	Homo_sapiens.GRCh38.101.chr.gtf
    output:
  		mapping/{sample}.sor.bam # quelle est l'extension 

    shell: ### docker run de STAR     

-----------------------------------------------


