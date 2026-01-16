#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters read1 and read2
params.forward_read = null
params.reverse_read = null

// Define parameters for trimmomatic
params.thread = 4
params.phred = 33
params.PATH_TO_ADAPTER_CONTAM_FILE = 'Path/to/your/adapter/file'
params.leading = 3
params.trailing = 3
params.slidingwindow = '4:15'
params.minlength = 36

// Define parameters for flash
params.max_overlap = 150

// Define parameters for plentyofbugs
params.assembler = 'path/to/skesa'
params.reference_genome = 'path/to/your/reference/genome/dir'

// Define input parameter path for pad reads
params.pad_read_path = 'path/to/your/pad_reads.py'

// Define input parameters for AlignGraph
params.distancelow = 100
params.distancehigh = 1000

// Define input parameters for quast2 (add this with your other params)
params.minContigLength = 500

// Define output directories
params.outputDir1 = "fastqc_out"
params.outputDir2 = "multiqc_out"
params.outputDir3 = "trimmomatic_out"
params.outputDir4 = "flash_out"
params.outputDir5 = "unicycler_out"
params.outputDir6 = "quast_out"
params.outputDir7 = "plentyofbugs_out"
params.outputDir8 = "bowtie2_out"
params.outputDir9 = "reference_based_assembly"
params.outputDir10 = "busco_out"

// Define the first process (FastQC)
process fastqc {
    input:
    path forward_read
    path reverse_read

    output:
    path "${params.outputDir1}", emit: fastqc_out

    script:
    """
    mkdir -p ${params.outputDir1}
    fastqc -o ${params.outputDir1} -f fastq "$forward_read" "$reverse_read"
    """
}

// Define the second process (MultiQC)
process multiqc {
    input:
    path fastqc_out

    output:
    path "${params.outputDir2}", emit: multiqc_out

    script:
    """
    mkdir -p ${params.outputDir2}
    multiqc ${fastqc_out} -o ${params.outputDir2}
    """
}

// Define the third process (Trimmomatic)
process trimmomatic {
    input:
    path forward_read
    path reverse_read
    val thread
    val phred
    path PATH_TO_ADAPTER_CONTAM_FILE
    val leading
    val trailing
    val slidingwindow
    val minlength

    output:
    path "${params.outputDir3}", emit: trimmomatic_out

    script:
    """
    mkdir -p ${params.outputDir3}
    trimmomatic PE -threads ${thread} -phred${phred} "$forward_read" "$reverse_read" \
        ${params.outputDir3}/output_1P.fq ${params.outputDir3}/output_1U.fq \
        ${params.outputDir3}/output_2P.fq ${params.outputDir3}/output_2U.fq \
        ILLUMINACLIP:${PATH_TO_ADAPTER_CONTAM_FILE}:2:30:10 LEADING:${leading} \
        TRAILING:${trailing} SLIDINGWINDOW:${slidingwindow} MINLEN:${minlength}
    """
}

// Define the fourth process (FLASH)
process flash {
    input:
    path trimmomatic_out
    val max_overlap

    output:
    path "${params.outputDir4}", emit: flash_out

    script:
    """
    mkdir -p ${params.outputDir4}
    flash --max-overlap ${max_overlap} \
        ${trimmomatic_out}/output_1P.fq ${trimmomatic_out}/output_2P.fq \
        -d ${params.outputDir4}
    """
}

// Define the fifth process (Unicycler)
process unicycler {
    input:
    path trimmomatic_out
    path flash_out

    output:
    path "${params.outputDir5}/assembly.fasta", emit: assembly_file

    script:
    """
    mkdir -p ${params.outputDir5}
    unicycler -1 ${trimmomatic_out}/output_1P.fq -2 ${trimmomatic_out}/output_2P.fq \
        -s ${flash_out}/out.extendedFrags.fastq -o ${params.outputDir5}
    """
}

// Define the sixth process (Quast)
process quast {
    input:
    path assembly_file

    output:
    path "${params.outputDir6}", emit: quast_out

    script:
    """
    mkdir -p ${params.outputDir6}
    quast.py -o ${params.outputDir6} ${assembly_file}
    """
}

// Define plentyofbugs
process plentyofbugs {
    input:
    path quast_out
    path assembler
    path reference_genome
    path trimmomatic_out

    output:
    path "${params.outputDir7}/best_reference", emit: best_reference
    path "${params.outputDir7}/assembly/contigs.fasta", emit: contigs_file

    script:
    """
    plentyofbugs --assembler ${assembler} \
        -f ${trimmomatic_out}/output_1P.fq \
        -r ${trimmomatic_out}/output_2P.fq \
        -g ${reference_genome} \
        -o ${params.outputDir7}
    """
}

// Define process bowtie2_build
process bowtie2_build {
    input:
    path best_reference
    path reference_genome

    output:
    path "${params.outputDir8}/reference_index*", emit: bowtie2_index

    script:
    """
    mkdir -p ${params.outputDir8}
    best_genome=\$(awk '{gsub(/.*\\//, ""); print \$1}' "${best_reference}")
    genome_file="${reference_genome}/\${best_genome}"
    bowtie2-build -f "\${genome_file}" "${params.outputDir8}/reference_index"
    """
}

// Define process bowtie2
process bowtie2 {
    input:
    path bowtie2_index
    path trimmomatic_out

    output:
    path "${params.outputDir8}/out.sam", emit: sam_file

    script:
    """
    mkdir -p ${params.outputDir8}
    indexPath=\$(find -L ${bowtie2_index} -type f -name '*.rev.1.bt2' | sed 's/\\.rev.1.bt2\$//')
    echo "Using index path: \${indexPath}"
    bowtie2 -x \${indexPath} \
            -1 ${trimmomatic_out}/output_1P.fq \
            -2 ${trimmomatic_out}/output_2P.fq \
            -S "${params.outputDir8}/out.sam"
    """
}

// Define the process for seqtk
process seqtk {
    input:
    path sam_file
    path trimmomatic_out

    output:
    path "${params.outputDir9}/paired_forward_1.fa", emit: forward_fa
    path "${params.outputDir9}/paired_forward_2.fa", emit: reverse_fa

    script:
    """
    mkdir -p ${params.outputDir9}
    seqtk seq -A ${trimmomatic_out}/output_1P.fq > ${params.outputDir9}/paired_forward_1.fa
    seqtk seq -A ${trimmomatic_out}/output_2P.fq > ${params.outputDir9}/paired_forward_2.fa
    """
}

// Define process for pad read
process pad_read {
    input:
    path forward_fa
    path reverse_fa
    path pad_read_path

    output:
    path "${params.outputDir9}/padded_out1.fa", emit: padded_file1
    path "${params.outputDir9}/padded_out2.fa", emit: padded_file2

    script:
    """
    mkdir -p ${params.outputDir9}
    python ${pad_read_path} ${forward_fa} ${params.outputDir9}/padded_out1.fa 150
    python ${pad_read_path} ${reverse_fa} ${params.outputDir9}/padded_out2.fa 150
    """
}

// Define the process for AlignGraph
process AlignGraph {
    input:
    path best_reference
    path reference_genome
    path padded_file1
    path padded_file2
    path assembly_file
    val distancelow
    val distancehigh

    output:
    path "${params.outputDir9}/sample_extendedcontig.fasta", emit: extended_contigs
    path "${params.outputDir9}/sample_remainingcontig.fasta", emit: remaining_contigs

    script:
    """
    # Create output directory
    mkdir -p ${params.outputDir9}

    # Resolve the genome file path properly
    best_genome=\$(awk '{gsub(/.*\\//, ""); print \$1}' "${best_reference}")
    genome_path="\$(readlink -f "${reference_genome}/\${best_genome}")"

    if [ ! -f "\${genome_path}" ]; then
        echo "ERROR: Genome file not found at \${genome_path}"
        echo "Tried to resolve: ${reference_genome}/\${best_genome}"
        echo "Directory contents:"
        ls -l "${reference_genome}"
        exit 1
    fi

    echo "Using genome file: \${genome_path}"
    AlignGraph \\
        --read1 "${padded_file1}" \\
        --read2 "${padded_file2}" \\
        --contig "${assembly_file}" \\
        --genome "\${genome_path}" \\
        --distanceLow ${distancelow} \\
        --distanceHigh ${distancehigh} \\
        --extendedContig "${params.outputDir9}/sample_extendedcontig.fasta" \\
        --remainingContig "${params.outputDir9}/sample_remainingcontig.fasta"
    """
}

// Define the process for quast2 (final Quast)
process quast2 {
    input:
    path extended_contigs
    path remaining_contigs
    val min_length

    output:
    path "${params.outputDir9}/quast_output", emit: quast2_out

    script:
    """
     mkdir -p ${params.outputDir9}/quast_output

    # Check if extended_contigs exists AND has at least one contig
    if [[ -s "${extended_contigs}" ]] && grep -q ">" "${extended_contigs}"; then
        quast.py -o ${params.outputDir9}/quast_output --min-contig ${min_length} "${extended_contigs}"
    else
        echo "WARNING: Using remaining_contigs (extended_contigs was empty/invalid)"
        quast.py -o ${params.outputDir9}/quast_output --min-contig ${min_length} "${remaining_contigs}"
    fi
    """
}

process busco {
    // Simple process with direct BUSCO execution
    cpus 4  // Sets default CPU count

    input:
    path extended_contigs
    path remaining_contigs

    output:
    path "busco_output/*", emit: busco_results
    path "busco_output/short_summary.*.txt", emit: summary

    script:
    """
    # Select input file (extended contigs preferred)
    if [[ -s "${extended_contigs}" ]] && grep -q ">" "${extended_contigs}"; then
        INPUT="${extended_contigs}"
        echo "Using extended contigs as input"
    else
        INPUT="${remaining_contigs}"
        echo "Using remaining contigs as input"
        [[ -s "\$INPUT" ]] || { echo "ERROR: No valid input files"; exit 1; }
    fi

    # Run BUSCO
    busco \\
        -i "\$INPUT" \\
        -o busco_output \\
        -l bacteria_odb10 \\
        -m genome \\
        -c ${task.cpus} \\
        --force
    """
}

workflow {
    // Run fastqc and get the output
    fastqc_ch = fastqc(params.forward_read, params.reverse_read)

    // Run MultiQC on the FastQC output directory
    multiqc_ch = multiqc(fastqc_ch.fastqc_out)

    // Run trimmomatic
    trimmomatic_ch = trimmomatic(
        params.forward_read,
        params.reverse_read,
        params.thread,
        params.phred,
        params.PATH_TO_ADAPTER_CONTAM_FILE,
        params.leading,
        params.trailing,
        params.slidingwindow,
        params.minlength
    )

    // Run flash
    flash_ch = flash(trimmomatic_ch.trimmomatic_out, params.max_overlap)

    // Run unicycler
    unicycler_ch = unicycler(trimmomatic_ch.trimmomatic_out, flash_ch.flash_out)

    // Run quast
    quast_ch = quast(unicycler_ch.assembly_file)

    // Run plentyofbugs
    plentyofbugs_ch = plentyofbugs(
        quast_ch.quast_out,
        params.assembler,
        params.reference_genome,
        trimmomatic_ch.trimmomatic_out
    )

    // Run bowtie2_build process
    bowtie2_build_ch = bowtie2_build(
        plentyofbugs_ch.best_reference,
        params.reference_genome
    )

    // Run bowtie2 process
    bowtie2_ch = bowtie2(
        bowtie2_build_ch.bowtie2_index,
        trimmomatic_ch.trimmomatic_out
    )

    // Run seqtk
    seqtk_ch = seqtk(bowtie2_ch.sam_file, trimmomatic_ch.trimmomatic_out)

    // Run the pad_read process
    pad_read_ch = pad_read(
        seqtk_ch.forward_fa,
        seqtk_ch.reverse_fa,
        params.pad_read_path
    )

    // Run AlignGraph
    aligngraph_ch = AlignGraph(
        plentyofbugs_ch.best_reference,
        params.reference_genome,
        pad_read_ch.padded_file1,
        pad_read_ch.padded_file2,
        unicycler_ch.assembly_file,
        params.distancelow,
        params.distancehigh
    )

    // Run quast2 for AlignGraph
    quast2_ch = quast2(aligngraph_ch.extended_contigs, aligngraph_ch.remaining_contigs, params.minContigLength)

    // Run busco
    busco_ch = busco(aligngraph_ch.extended_contigs, aligngraph_ch.remaining_contigs)
}
(base) manish_kumar@DESKTOP-G29M48F:~/scripts/DNRBG$ ls -lhrt
total 60K
-rwxr-xr-x 1 manish_kumar manish_kumar  12K Jan 16 06:32 dnrgb_new_deep.nf
-rwxr-xr-x 1 manish_kumar manish_kumar  26K Jan 16 06:33 dnrgb.sh
-rwxr-xr-x 1 manish_kumar manish_kumar  701 Jan 16 06:33 pad_reads_1.py
-rwxr-xr-x 1 manish_kumar manish_kumar  205 Jan 16 06:33 env.yml
-rw-r--r-- 1 manish_kumar manish_kumar 9.8K Jan 16 06:39 Dockerfile
(base) manish_kumar@DESKTOP-G29M48F:~/scripts/DNRBG$ cat dnrgb.sh
#!/bin/bash
###############################################
#colour coding format for echo
green='\033[0;32m'  #colour coding for green
red='\033[0;31m'    #colour coding for red
yellow='\033[1;33m' #colour coding for yellow
reset='\033[0m'     #reset colour
##############################################
# Help
##############################################
Help()
{
# Colors for formatting
green='\033[0;32m'
yellow='\033[1;33m'
reset='\033[0m'
# Display Help
echo -e "${green}     This is an integrative genome assembly pipeline.${reset}"
echo -e "${yellow}    Usage:${reset}"
echo -e "${yellow}    ./genome_assembly.sh read_1 read_2${reset}"
echo -e "${yellow}    where:${reset}"
echo -e "${yellow}    read_1 - Path to the first input read file (fastq format)${reset}"
echo -e "${yellow}    read_2 - Path to the second input read file (fastq format)${reset}"
echo -e "${green}     Please ensure that the input files are in fastq format.${reset}"
}
###############################################
#########pre procressing step one ##############
#get working directory
current_dir=$(pwd)
#input files that is raw reads in fastq format
forward_read=$1
reverse_read=$2
#check if the arguments are given
if [ "$#" -eq 2 ]; then # check for two arguments
echo -e "${green}Input provided....${reset}"
else
echo -e "${red}Proper input is not provided please provide raw reads...${reset}"
Help # Call the help function to display usage
exit 1
fi
######################################################
#Read Quality Assessment using FastQC
if command -v fastqc &>/dev/null; then # check for the fastqc programme
echo -e "${green}fastqc programme in installed...${reset}" # if found run fastqc here
mkdir fastqc_out && # make fastqc output folder
export LC_ALL=C &&
echo -e "${green}executing fastqc...${reset}" &&
fastqc -o fastqc_out -f fastq "$forward_read" "$reverse_read" || exit 1 # running fastqc .....
else
# if not found ask for path and run fastqc here
echo -e "${red}fastqc programme is not found please check this programme is eighter installed or present in your path...${reset}"
echo -e "${yellow}please provide absolute path to the executable....${reset}"
read -rp "Pleae provide absolue path to fastqc executable here : " fastqc
echo -e "${green}Making fastqc output folder...${reset}"
export LC_ALL=C &&
mkdir fastqc_out && # make fastqc output folder
echo -e "${green}executing fastqc...${reset}"
"$fastqc" -o fastqc_out -f fastq "$forward_read" "$reverse_read" || exit 1
wait for command to run using wait command
fastqc_pid=$!
wait $fastqc_pid
fi || exit
echo -e "${green} fastqc run sucessfully..."
#############################################################
# use multiqc to compile results
mkdir multiqc_out && #create a directory named multiqc_out
# navigate into the directory named multiqc_out
if command -v multiqc &>/dev/null; then # check for multiqc programme
echo -e "${green} multiqc programme is installed...${reset}" # if found run multiqc here
echo -e "${green} executing multiqc...${reset}"
multiqc  "$current_dir/fastqc_out" -o "$current_dir/multiqc_out" || exit 1
else # if not found ask for path and run multiqc here
echo -e "${red} multiqc programme is not found please check this programme is eighter installed or present in your path...${reset}"
echo -e "${yellow} please provide absolute path to the multiqc executable....${reset}"
read -rp "Pleae provide absolute path to multiqc executable here :" multiqc
echo -e "${green} executing multiqc...${reset}"
"$multiqc" -o "$current_dir/fastqc_out" ../multiqc_out || exit 1
cd "$current_dir" &&
# Wait for command to run using wait command
multiqc_pid=$!
wait "$multiqc_pid"
fi || exit
echo -e "${green} multiqc run sucessfully...${reset}"
###############################################################
# Use trimmomatic to trim unpaired reads and remove adapter sequences
# Check for java
if command -v java &>/dev/null; then
echo -e "${green} java installed...${reset}"
else # if not installed ask for installation
echo -e "${red}Please install java and check if installed already check if it is present in your path..."
fi || exit
# Make directory for trimmomatic
echo -e "${green}makeing output directory for trimmomatic...${reset}"
mkdir trimmomatic_out &&
# Check for trimmomatic
# Taking parameters using if else statement
# Navigate to trimmomatic_out directory and captute its path
if command -v trimmomatic &>/dev/null; then # check for trimmomatic
echo -e "${green} Trimmomatic programme is installed...${reset}"
echo -e "${yellow} Please provide the relevent parameters for the tools...${reset}"
read -rp "Please provide number of threads to run the tool : " thread
read -rp "Please provide phread score : " phread
read -rp "Please provide path to adapter file:" PATH_TO_ADAPTER_CONTAM_FILE
read -rp "Please provide parameter LEADING: " leading
read -rp "Please provide parameter TRAILING: " trailing
read -rp "Please provide parameter SLIDINGWINDOW: " slidingwindow
read -rp "Please provide parameter MINLEN: " minlength
echo -e "${green} executing Trimmomatic...${reset}"
trimmomatic PE -threads "$thread" -phred"$phread" "$forward_read" "$reverse_read" output_1P.fq output_1U.fq output_2P.fq output_2U.fq ILLUMINACLIP:"$PATH_TO_ADAPTER_CONTAM_FILE":2:30:10 LEADING:"$leading" TRAILING:"$trailing" SLIDINGWINDOW:"$slidingwindow" MINLEN:"$minlength" || exit 1  #To remove low quality base pairs and adapter sequences
else # ask path to the trimmomatic programme path
echo -e "${red} please install trimmomatic and check if already installed check if it is present in your path...${reset}"
echo -e "${yellow} Please provide the absolute path to the trimmomatic executable here : ${reset}"
read -rp "Please provide your trimmomatic path here : " trimmomatic
echo -e "${yellow} Please provide the relevent parameters for the tools...${reset}"
read -rp "Please provide number of threads to run the tool : " thread
read -rp "Please provide phread score : " phread
read -rp "Please provide path to adapter file :" PATH_TO_ADAPTER_CONTAM_FILE
read -rp "Please provide parameter LEADING: " leading
read -rp "Please provide parameter TRAILING: " trailing
read -rp "Please provide parameter SLIDINGWINDOW: " slidingwindow
read -rp "Please provide parameter MINLEN: " minlength
echo -e "${green}executing Trimmomatic...${reset}"
java -jar "$trimmomatic" PE -threads "$thread" -phred"$phread" "$forward_read" "$reverse_read" output_1P.fq output_1U.fq output_2P.fq output_2U.fq ILLUMINACLIP:"$PATH_TO_ADAPTER_CONTAM_FILE":2:30:10 LEADING:"$leading" TRAILING:"$trailing" SLIDINGWINDOW:"$slidingwindow" MINLEN:"$minlength" || exit 1   #To remove low quality base pairs and adapter sequences
fi || exit
# Wait for the command to be compleated
trimmomatic_pid=$!
wait $trimmomatic_pid
# Move all output  to trimmomatic_out directory
mv output_1P.fq output_1U.fq output_2P.fq output_2U.fq trimmomatic_out/ &&
# Capturing path of trimmomatic_out directory
cd trimmomatic_out && # navigate to trimmomatic directory
trimmomatic_path=$(pwd) && # capturing directory path of trimmed reads
cd "$current_dir" || exit 1
#########################################################################################################
# Merging of Overlapping Paired-end Reads using FLASH
# Checking if flash is installed
# Using if else statement to check and run command
# Make directory for flash output
mkdir flash_out && #make flash directory
cd flash_out &&    #navigate to flash directory
if command -v flash &>/dev/null; then
echo -e "${green}flash is installed...${reset}"
echo -e "${green}executing flash...${reset}"
read -rp "Please provide input for parameter maximum overlap: " max_overlap
flash --max-overlap "$max_overlap" "$trimmomatic_path/output_1P.fq" "$trimmomatic_path/output_2P.fq" > output -d "$PWD" || exit 1 #flash will generate notCombined and extended fastq files
else # ask for the absolute path to the flash command
echo -e "${red} flash is not installed if you think it is installed please check if it is in your path...${reset}"
echo -e "${yellow}Please provide absolute path to the flash command...${reset}"
read -rp "Please provide absolute path to the flash command here : " flash
read -rp "Please provide input for parameter maximum overlap: " max_overlap
echo -e "${green}executing flash...${reset}"
"$flash" --max-overlap "$max_overlap" "$trimmomatic_path/output_1P.fq" "$trimmomatic_path/output_2P.fq" > output -d "$PWD" || exit 1 #flash will generate notCombined and extended fastq file
fi || exit
# Wait for the command to compleated
flash_pid=$!
wait $flash_pid
# Capture flash path
flash_path=$(pwd) &&  #use this in unicycler command to get extended frags
# Change to the old directory
cd "$current_dir" &&
#############################################################################################################
# De novo Genome Assembly using Unicycler
# Check if unicycler is installed
# Using if else statement to check and run command
# Make directory for unicycler output
# Make directory for unicycler output
mkdir unicycler_out && cd unicycler_out || exit
# Check if both unicycler and spades are installed
if command -v spades &>/dev/null && command -v unicycler &>/dev/null; then
echo -e "${green}Both spades and unicycler are installed.${reset}"
               unicycler  -1 "$trimmomatic_path/output_1P.fq" \
                          -2 "$trimmomatic_path/output_2P.fq" -s "$flash_path/out.extendedFrags.fastq" \
                          -o assembly  || exit 1
                                     fi || exit
# Check if neither unicycler nor spades are installed
if ! command -v spades &>/dev/null && ! command -v unicycler &>/dev/null; then
echo -e "${red}Both unicycler and spades are not installed.${reset}"
echo -e "${yellow}Please provide absolute paths to unicycler and spades.${reset}"
read -rp "Please provide absolute path to unicycler: " unicycler
echo -e "${green}Running unicycler.${reset}"
        "$unicycler" -1 "$trimmomatic_path/output_1P.fq" \
                    -2 "$trimmomatic_path/output_2P.fq" -s "$flash_path/out.extendedFrags.fastq" \
                    -o assembly || exit 1
fi || exit
cd assembly &&
assemblypath=$(pwd)
# Wait for unicycler command to complete
unicycler_pid=$!
wait $unicycler_pid
# Navigate to old directory
cd "$current_dir" || exit
#######################################################################################
# Assessment of Assembly Quality using QUAST
# Make quast output directory
mkdir quast_out &&
cd quast_out &&
# Check if quast is installed
if command -v quast &>/dev/null; then
echo -e "${green}quast is installed...${reset}"
echo -e "${green}Running quast...${reset}"
quast -o quast_output "$assemblypath/assembly.fasta" || exit 1
else # take absolute path to quast
echo -e "${yellow}Please provide absolute path to quast.py and make sure python is installed...${reset}"
read -rp "Please provide absolute path to quast here : " quast
echo -e "${green}Running quast...${reset}"
python "$quast" -o quast_output "$assemblypath/assembly.fasta" || exit 1
fi || exit
# Wait for quast command to complete
quast_pid=$!
wait $quast_pid
# Print completion message
echo -e "${green} Phase 1 complete!${reset}"
##################################################################################################
# Phase 2
##################################################################################################
cd "$current_dir" || exit
cd trimmomatic_out || exit # navigate to trimmomatic_out directory
trimmomatic_output_path=$(pwd) # capturing path
cd "$current_dir" || exit
# ask for path to reference genome
echo -e "${yellow} Please provide path to reference genome...${reset}"
read -rp "Please provide your path here : " reference_genome
#check for required programme
# List of required programs
required_programs=("mash" "skesa" "spades" "seqtk" "unicycler" "plentyofbugs")
# Check if each required program is installed
missing_programs=()
for program in "${required_programs[@]}"; do
if ! command -v "$program" &> /dev/null; then
missing_programs+=("$program")
fi
done
# If there are missing programs, install them
if [ ${#missing_programs[@]} -gt 0 ]; then
echo -e "${red}The following programs are missing: ${missing_programs[*]}${reset}"
echo -e "${missing_programs[@]}"
echo -e "${yellow}Please install the following programmes.${reset}"
echo -e "${missing_programs[@]}"
######################################################################
elif [ ${#missing_programs[@]} == 0 ]; then
######################################################################
echo -e "${yellow}Please provide name of assembler skesa or spades.${reset}"
read -rp "Please provide your assembler name here: " assembler
echo -e "${green}Running plenty of bugs.${reset}"
plentyofbugs --assembler "$assembler" \
             -f "${trimmomatic_output_path}/output_1P.fq" \
             -r "${trimmomatic_output_path}/output_2P.fq" \
             -g "$reference_genome" \
             -o plentyofbugs_out
######################################################################
elif  ! command -v plentyofbugs &>/dev/null && ! command -v spades && ! command -v skesa && command -v seqtk && command -v unicycler && command -v seqtk && command -v mash; then
######################################################################
# ask for absolute path to unicycler and skesa or spades
echo -e "${yellow}Please provide absolute path to plentyofbbugs.${reset}"
read -rp "${yellow}Please provide absolute path to plentyofbugs here : " plentyofbugs
echo -e "${yellow}Please provide absolute path to assembler skesa or spades.${reset}"
read -rp "Please provide your assembler path here : " assembler
"$plentyofbugs" --assembler "$assembler" \
                -f "${trimmomatic_output_path}/output_1P.fq" \
                -r "${trimmomatic_output_path}/output_2P.fq" \
                -g "$reference_genome" \
                -o plentyofbugs_out
fi || exit
#######make function###########################
function activate_conda() {
# If statement to check if the environment exists
if conda env list | grep -q "install_env"; then
echo -e "${green} Enviornment 'install_env' already exists. Skipping environment creation...${reset}"
# Initiate conda environment
echo -e "${yellow} Initiating conda environment...${reset}" &&
source activate install_env &&
# Install plentyofbugs command not found using pip
if ! command -v plentyofbugs; then
echo -e "${red}plentyofbugs not found.${reset}"
echo -e "${green}trying installing plentyofbugs using pip.${reset}"
pip install plentyofbugs
fi
# List of required programs
required_programs=("mash" "skesa" "spades" "seqtk" "unicycler" "plentyofbugs")
# Check if each required program is installed
missing_programs=()
for program in "${required_programs[@]}"; do
if ! command -v "$program" &> /dev/null; then
missing_programs+=("$program")
fi
done
# If there are missing programs, install them
if [ ${#missing_programs[@]} -gt 0 ]; then
echo -e "${yellow}The following programs are missing and will be installed: ${missing_programs[*]}${reset}"
echo "${missing_programs[@]}"
# Install missing programs
for program in "${missing_programs[@]}"; do
echo -e "${yellow}Installing $program...${reset}"
conda install --yes "bioconda::$program"
done
fi
# run plenty of bugs
echo -e "${green}Running plentyofbugs.${reset}"
echo -e "${yellow}Please provide name of assembler skesa or spades.${reset}"
read -rp "Please provide your assembler name here: " assembler
echo -e "${green}Running plenty of bugs.${reset}"
plentyofbugs --assembler "$assembler" \
             -f "${trimmomatic_output_path}/output_1P.fq" \
             -r "${trimmomatic_output_path}/output_2P.fq" \
             -g "$reference_genome" \
             -o plentyofbugs_out
#conda deactivate
conda deactivate
else
# Ask if the following command will be installed or not
echo -p "${yellow} Would you like to install the following in conda environment? Please provide an answer with yes or no below !.${reset}"
read -rp "Please provide your answer here : " answer
# Use if else statement to evaluate the command line input
if [ "$answer" == yes ]; then
# Create the conda environment
echo -e "${yellow} Creating conda environment...${reset}"
conda create --name install_env &&
# Initiate conda environment
echo -e "${yellow} Initiating conda environment...${reset}"
source activate install_env &&
# Install plentyofbugs command not found using pip
if ! command -v plentyofbugs; then
echo -e "${red}plentyofbugs not found.${reset}"
echo -e "${green}trying installing plentyofbugs using pip.${reset}"
pip install plentyofbugs
fi
# List of required programs
required_programs=("mash" "skesa" "spades" "seqtk" "unicycler" "plentyofbugs")
# Check if each required program is installed
missing_programs=()
for program in "${required_programs[@]}"; do
if ! command -v "$program" &> /dev/null; then
missing_programs+=("$program")
fi
done
# If there are missing programs, install them
if [ ${#missing_programs[@]} -gt 0 ]; then
echo -e "${yellow}The following programs are missing and will be installed: ${missing_programs[*]}${reset}"
echo -e "${yellow}The following programs are missing and will be installed:${reset}"
echo "${missing_programs[@]}"
# Install missing programs
for program in "${missing_programs[@]}"; do
echo -e "${yellow}Installing $program..."
conda install --yes "bioconda::$program"
done
echo -e "${green}Installation of required programs complete."
# Run Plentyofbugs
plentyofbugs --assembler "$assembler" \
             -f "${trimmomatic_output_path}/output_1P.fq" \
             -r "${trimmomatic_output_path}/output_2P.fq" \
             -g "$reference_genome" \
            -o plentyofbugs_out &&  #Selection of best genome based upon Mash Value
conda deactivate
else
echo -e "${green}All required programs are already installed. Skipping installation."
fi
echo -e "${green} Installation complete...${reset}"
# Run command plenty of bugs
echo -e "${yellow} Please select assembler skesa or spades...${reset}"
read -rp "Please provide your assembler skesa/spades here : " assembler
plentyofbugs --assembler "$assembler" \
             -f "${trimmomatic_output_path}/output_1P.fq" \
             -r "${trimmomatic_output_path}/output_2P.fq" \
             -g "$reference_genome" \
             -o plentyofbugs_out &&  #Selection of best genome based upon Mash Value
conda deactivate
else
echo -e "${red} ERROR: Please check all the tools are properly installed to run the program."
fi
fi
 }
# Check if all required commands are missing
if ! command -v plentyofbugs &>/dev/null && \
   ! command -v spades &>/dev/null && \
   ! command -v skesa &>/dev/null && \
   ! command -v seqtk &>/dev/null && \
   ! command -v unicycler &>/dev/null && \
   ! command -v mash &>/dev/null; then
     echo -e "${red}There are some programs missing and will be installed and run in a separate conda environment.${reset}"
     echo -e "${yellow}Would you like to install them in a separate conda environment? If the problem persists, you can manually install them in the 'install_env' environment and then rerun the assembler.${reset}"
     read -rp "${yellow}Please provide your answer here (yes/no): " yesno
if [ "$yesno" == "yes" ]; then
activate_conda
else
echo -e "${red}Please install the programs manually!"
fi
else
echo -e "${green}status ...................ok!.${reset}"
fi
############################################################################################
########################Percentage alignment calculation using Bowtie2
# Change to current directory
cd "$current_dir" || exit
# Taking path to reference genome
cd plentyofbugs_out &&
reference_path=$(cat best_reference | awk '{print $1}') &&
# Change path to current directory
cd "$current_dir" || exit
# Make bowtie2 directory
mkdir bowtie2_out || exit
# Change directory to bowtie2
cd bowtie2_out  || exit
# Do percentage alignment calculations
if command -v bowtie2  &>/dev/null; then
echo -e "${green}Bowtie2 is installed${reset}"
echo -e "${green}Executing Bowtie2${reset}"
echo -e "${green}Indexing reference genome${reset}"
# Indexing reference genomes
read -rp "Please provide full path to plentyofbugs best reference here : " reference_path_plentyofbugs
bowtie2-build -f "$reference_path_plentyofbugs"  reference_index &&
bowtie2build_pid=$!
wait $bowtie2build_pid
# Change to current directory
cd "$current_dir" || exit
# Change path to trimmomatic_out
cd trimmomatic_out &&
trimmomatic_output_path=$(pwd) &&
# Change to current directory
cd "$current_dir" || exit
# Change path to bowtie2_out
cd bowtie2_out || exit
bowtie2 -x  reference_index -1 "$trimmomatic_output_path/output_1P.fq" -2 "$trimmomatic_output_path/output_2P.fq" -S out.sam &&
bowtie2_pid=$!
wait $bowtie2_pid
else
echo "${red}Error: Bowtie2 is not installed. Please install Bowtie2 and try again.${reset}"
exit 1
# Change to current directory
cd "$current_dir" || exit
fi || exit
#########################################################################
# Check if all the command exist
#check for seqtk if not found store its path in a variable
if command -v seqtk &>/dev/null; then
echo -e "${green} seqtk is installed...${reset}"
else
# ask for path
echo "${red} Please install seqtk...${reset}"
fi || exit
#check for AlignGraph if not found store its path in a variable
if command -v AlignGraph &>/dev/null; then
echo -e "${green} AlignGraph is installed...${reset}"
else
# ask to AlignGraph
echo -e "${red} Please install AlignGraph...${reset}"
fi || exit
####################################################################################
# Store the original directory
cd "$current_dir" || exit
##################################################################################
# Reference Based Assembly
cd plentyofbugs_out || exit
reference_genome=$(cat best_reference | awk '{print $1}')
cd "$current_dir" || exit
cd "$current_dir"/plentyofbugs_out/assembly/ || exit
contigpath=$(pwd)
cd "$current_dir" || exit
cd trimmomatic_out &&
trimmed_read=$(pwd)
cd "$current_dir" || exit
# make directory for reference based assembly
mkdir reference_based_assembly &&
cd reference_based_assembly || exit
reference_based_assembly=$(pwd)
if command -v seqtk &>/dev/null && command -v AlignGraph &>/dev/null; then
echo -e "${green}Running seqtk to convert fastq files to fasta files${reset}"
seqtk seq -A "$trimmed_read/output_1P.fq" > output_1P.fa &&
seqtk seq -A "$trimmed_read/output_2P.fq" > output_2P.fa &&
seqtk_pid=$!
wait $seqtk_pid
#######################################################
read -rp "please provide absolute pad_reads_1.py path": pad
read -rp "please provide your raw read length input for parameter: " Pad_read
python "$pad"  output_1P.fa padded_out1.fa "$Pad_read"
python "$pad"  output_2P.fa padded_out2.fa "$Pad_read"
#######################################################
echo -e "${yellow} Please provide required input parameters for AlignGraph${reset}"
read -rp "Please provide your input for parameter --distanceLow : "  distancelow
read -rp "please provide your input for parameter --distancehigh : " distancehigh
#read -rp "Please provide unicycler contig with absolute path  from unicycler --contig : " contigs
read -rp "Please provide absolute path of best reference from plentyofbugs : "  best_reference_path
echo -e "${green} Running AlignGraph${reset}"
AlignGraph --read1 padded_out1.fa --read2 padded_out2.fa --contig "$assemblypath/assembly.fasta" --genome "$best_reference_path" --distanceLow "$distancelow" --distanceHigh "$distancehigh" --extendedContig "$reference_based_assembly/sample_extendedcontig.fasta" --remainingContig "$reference_based_assembly/sample_remainingcontig.fasta"
else
echo "${red}please install AlignGraph or seqtk if already installed please make sure they are in you current path${reset}"
fi
cd "$current_dir" || exit
###################################################
#use quast
#make quast output dir
mkdir quast2_out &&
cd quast2_out &&

# Check if quast is installed
if command -v quast &>/dev/null; then
    echo -e "${green}quast is installed...${reset}"
    echo -e "${green}Running quast...${reset}"

    file="${reference_based_assembly}/sample_extendedcontig.fasta"

    # Check if file does not exist or is empty
    if [[ ! -e "$file" || ! -s "$file" ]]; then
        echo -e "${red}Warning: sample_extendedcontig.fasta is empty. sample_remainingcontig.fasta will be used${reset}"
        echo -e "${green}Running quast...${reset}"
        read -rp "Please provide minimum contig length for running quast here: " contig_length_quast

        quast -o quast_output --min-contig "$contig_length_quast" "${reference_based_assembly}/sample_remainingcontig.fasta" || exit 1
    else
        read -rp "Please provide minimum contig length for running quast here: " contig_length_quast
        echo -e "${green}Running quast...${reset}"
        quast -o quast_output --min-contig "$contig_length_quast" "${reference_based_assembly}/sample_extendedcontig.fasta" || exit 1
    fi
else
    echo -e "${red}quast is not installed. Please install it first.${reset}"
    exit 1
fi
cd "$current_dir" || exit
###################################################
#####################################################
# Change to reference-based assembly directory and get its absolute path
cd reference_based_assembly || exit
reference_based_assembly_path=$(pwd) || exit
cd "$current_dir" || exit

# Create and enter BUSCO output directory
mkdir -p busco || exit
cd busco || exit

# Check if BUSCO is installed
if command -v busco &>/dev/null; then
    echo -e "${green}BUSCO is installed.${reset}"

    # Check if sample_extendedcontig.fasta exists and is not empty
    busco_input="${reference_based_assembly_path}/sample_extendedcontig.fasta"
    if [[ ! -e "$busco_input" || ! -s "$busco_input" ]]; then
        echo -e "${red}Error: sample_extendedcontig.fasta not found or empty in ${reference_based_assembly_path}.${reset}"
        exit 1
    fi

    # Ask user for lineage and run BUSCO
    read -rp "Specify the name of the BUSCO lineage to be used: " busco_lineage
    echo -e "${green}Running BUSCO on sample_extendedcontig.fasta...${reset}"

    busco -i "$busco_input" -m genome -l "$busco_lineage" -o busco_output || exit 1

else
    echo -e "${red}BUSCO is not installed. Please install it before proceeding.${reset}"
    exit 1
fi

# Return to the original directory
cd "$current_dir" || exit
#####################################################################
