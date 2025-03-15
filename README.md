# DNRBG-An-integrative-pipeline-for-assembly-of-bacterial-genomes
An integrative pipeline for assembling bacterial genomes.

This pipeline is used for de-novo reference-guided assembly, which operates in three phases.

1) Pre-processing and Contig Generation: The process begins with raw reads in FASTQ format. We first assess the quality of these reads, removing any low-quality sequences and Illumina adapters. The remaining short reads are assembled into contigs.

2) Selection of best reference: The best reference genome is identified, and the percentage alignment calculation is performed.

3) Reference-based scaffolding: Reference-based scaffolding is conducted.

*Installation*

This pipeline is available as a LINUX executable.**

*Dependencies*

This pipeline requires various dependencies.

fastQC Available at https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

multiQC Available at https://multiqc.info/docs/getting_started/installation/

Trimmomatic Available at https://github.com/usadellab/Trimmomatic

FLASH Available at https://ccb.jhu.edu/software/FLASH/

Unicycler Available at https://github.com/rrwick/Unicycler

QUAST Available at https://github.com/ablab/quast

Plentyofbugs Available at https://github.com/nickp60/plentyofbugs

Bowtie2 Available at https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

AlignGraph Available at https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

BUSCO Available at https://busco.ezlab.org/

*USAGE*

./dnrgb.sh read1.fastq read2.fastq 
where fastq_1 and fastq_2 are raw sequences, paired-end read one and two, respectively.

Docker To install docker you can see https://docs.docker.com/engine/install/ubuntu/. To use docker please use pull command docker pull somil131/dnrgb_genome_assembly:latest.

*Docker Command*

please make sure that your reference genome directory is named as reference_genome.

You can use the command as mentioned below:-

docker run -it -v $(pwd):/data -v /Path/to/reference_genome:/data/reference_genome somil131/dnrgb_genome_assembly:latest /data/path/to/Your_files_1.fastq /data/path/to/Your_files_2.fastq
