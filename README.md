# DNRBG: An Integrative Pipeline for Bacterial Genome Assembly

DNRBG is an integrative pipeline for de novo and reference-guided bacterial genome assembly. The pipeline operates in three phases.

1. **Pre-processing and contig generation:** Raw paired-end FASTQ reads are assessed for quality, and low-quality sequences and Illumina adapters are removed. The processed reads are then assembled into contigs.
2. **Best-reference selection:** The best candidate reference genome is identified, and the percentage of alignment is calculated.
3. **Reference-based scaffolding:** Reference-based scaffolding is performed to generate a reference-guided assembly.

## Dependencies

The pipeline uses the following tools:

- **FastQC:** https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **MultiQC:** https://multiqc.info/docs/getting_started/installation/
- **Trimmomatic:** https://github.com/usadellab/Trimmomatic
- **FLASH:** https://ccb.jhu.edu/software/FLASH/
- **Unicycler:** https://github.com/rrwick/Unicycler
- **QUAST:** https://github.com/ablab/quast
- **PlentyofBugs:** https://github.com/nickp60/plentyofbugs
- **Bowtie2:** https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
- **AlignGraph:** https://github.com/baoe/AlignGraph
- **BUSCO:** https://busco.ezlab.org/

Additional dependencies and executables may be required depending on the selected assembly configuration.

## Usage on Linux

The pipeline can be run using the `dnrbg.sh` shell script.

```bash
./dnrbg.sh \
  -1 /path/to/read1.fastq \
  -2 /path/to/read2.fastq \
  -g /path/to/reference_genome
```

Arguments:

- `-1`: Path to the forward FASTQ file.
- `-2`: Path to the reverse FASTQ file.
- `-g`: Directory containing candidate reference genome FASTA files.

Additional parameters can be specified as needed. Run the following command to view the available options:

```bash
./dnrbg.sh -h
```

## Running with Docker

### Install Docker

Install Docker by following the official instructions:

https://docs.docker.com/engine/install/ubuntu/

### Pull the Docker image

```bash
docker pull mkulab/dnrbg_latest:v1.0
```

### Docker command

Make sure that your candidate reference genome directory is named `reference_genome`. The directory should contain the reference FASTA files used for candidate-reference selection.

Run the following command from the directory containing your paired-end reads:

```bash
docker run -it \
  -v "$(pwd):/data" \
  -v /path/to/reference_genome:/data/reference_genome \
  mkulab/dnrbg_latest:v1.0 \
  -1 /data/read1.fastq \
  -2 /data/read2.fastq \
  -g /data/reference_genome
```

Replace `/path/to/reference_genome` with the actual path to the reference genome directory on your host system. Replace `read1.fastq` and `read2.fastq` with the actual FASTQ filenames located in the current working directory.

Additional options can be supplied using the flags supported by `dnrbg.sh`.

## Running with Nextflow

### Install Nextflow

Install Nextflow using the official documentation:

https://www.nextflow.io/docs/latest/install.html

**Required version:** `25.04.2.5947`

Check the installed version:

```bash
nextflow -version
```

### Run the pipeline

Run `dnrbg.nf` using the following command. Replace the placeholder paths with the actual paths to the reads, adapter file, assembler, reference genomes, and padding script.

```bash
nextflow run dnrbg.nf \
  --forward_read /path/to/forward.fastq \
  --reverse_read /path/to/reverse.fastq \
  --thread 4 \
  --phred 33 \
  --PATH_TO_ADAPTER_CONTAM_FILE /path/to/adapter/file \
  --leading 3 \
  --trailing 3 \
  --slidingwindow 4:15 \
  --minlength 36 \
  --max_overlap 150 \
  --assembler /path/to/skesa \
  --reference_genome /path/to/reference_genome \
  --pad_read_path /path/to/pad_reads.py \
  --distancelow 100 \
  --distancehigh 1000 \
  --outputDir1 fastqc_out \
  --outputDir2 multiqc_out \
  --outputDir3 trimmomatic_out \
  --outputDir4 flash_out \
  --outputDir5 unicycler_out \
  --outputDir6 quast_out \
  --outputDir7 plentyofbugs_out \
  --outputDir8 bowtie2_out \
  --outputDir9 reference_based_assembly \
  --outputDir10 busco
```

### Parameters and output directories

Provide the correct paths for all required input files and executables. Adjust parameters such as thread count, Phred score, trimming settings, FLASH overlap, assembler, and alignment distances according to your requirements.

The output directory parameters `--outputDir1` through `--outputDir10` correspond to the pipeline's individual output stages. Keep their names unchanged unless you also update the relevant workflow configuration.

Ensure that `dnrbg.nf` and `pad_reads.py` are available at their specified paths and that the required software dependencies are installed or otherwise accessible to the workflow.
