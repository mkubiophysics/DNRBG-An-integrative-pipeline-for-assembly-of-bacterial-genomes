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
params.outputDir10 = "busco"

// Define the first process (FastQC)
process fastqc {

    input:
    path forward_read
    path reverse_read

    output:
    path "${params.outputDir1}"

    script:
    """
    mkdir -p fastqc_out
    fastqc -o fastqc_out -f fastq "$forward_read" "$reverse_read"
    """
}

// Define the second process (MultiQC)
process multiqc {

    input:
    path fastqc_out

    output:
    path "${params.outputDir2}"

    script:
    """
    mkdir -p multiqc_out
    multiqc "$fastqc_out" -o ${params.outputDir2}
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
    path "${params.outputDir3}"

    script:
    """
    mkdir -p ${params.outputDir3}
    trimmomatic PE -threads "$thread" -phred"$phred" "$forward_read" "$reverse_read" \
        "${params.outputDir3}/output_1P.fq" "${params.outputDir3}/output_1U.fq" \
        "${params.outputDir3}/output_2P.fq" "${params.outputDir3}/output_2U.fq" \
        ILLUMINACLIP:"$PATH_TO_ADAPTER_CONTAM_FILE":2:30:10 LEADING:"$leading" \
        TRAILING:"$trailing" SLIDINGWINDOW:"$slidingwindow" MINLEN:"$minlength"
    """
}

// Define the fourth process (FLASH)
process flash {

    input:
    path trimmomatic_out
    val max_overlap

    output:
    path "${params.outputDir4}"

    script:
    """
    mkdir -p flash_out
    flash --max-overlap "$max_overlap" \
        "${params.outputDir3}/output_1P.fq" "${params.outputDir3}/output_2P.fq" \
        -d ${params.outputDir4}
    """
}

// Define the fifth process (Unicycler)
process unicycler {

    input:
    path trimmomatic_out
    path flash_out

    output:
    path "${params.outputDir5}"

    script:
    """
    mkdir -p unicycler_out
    unicycler -1 "${params.outputDir3}/output_1P.fq" -2 "${params.outputDir3}/output_2P.fq" \
        -s "${params.outputDir4}/out.extendedFrags.fastq" -o "${params.outputDir5}/assembly"
    """
}

// Define the sixth process (Quast)
process quast {

    input:
    path unicycler_out

    output:
    path "${params.outputDir6}"

    script:
    """
    mkdir -p quast_out
    quast.py -o "${params.outputDir6}/quast_output" "${params.outputDir5}/assembly/assembly.fasta"
    """
}

// Define plentyofbugs
process plentyofbugs {
    input:
    path quast_result
    path assembler
    path reference_genome
    path trimmomatic_out

    output:
    path "${params.outputDir7}/plentyofbugs_out/*", emit: plentyofbugs_result
    path "${params.outputDir7}/plentyofbugs_out/assembly/contigs.fasta" , emit: contigs_file

    script:
    """
    mkdir -p "${params.outputDir7}"
    plentyofbugs --assembler "$assembler" -f "${trimmomatic_out}/output_1P.fq" -r "${trimmomatic_out}/output_2P.fq" -g "$reference_genome" -o "${params.outputDir7}/plentyofbugs_out"
    """
}

// Define process bowtie2_build
process bowtie2_build {
    input:
    path plentyofbugs_result
    path reference_genome

    output:
    path "${params.outputDir8}/bowtie2_index/*", emit: bowtie2_index

    script:
    """
    mkdir -p "${params.outputDir8}/bowtie2_index"
    best_genome=\$(awk '{gsub(/.*\\//, ""); print \$1}' "best_reference")
    genome_file="${params.reference_genome}/\${best_genome}"
    bowtie2-build -f "\$genome_file" "${params.outputDir8}/bowtie2_index/reference_index"
    """
}

// Define process bowtie2
process bowtie2 {
    input:
    path bowtie2_index
    path trimmomatic_out

    output:
    path "${params.outputDir8}/out.sam"

    script:
    """
    mkdir -p bowtie2_out
    indexPath=\$(find -L ${bowtie2_index} -type f -name '*.rev.1.bt2' | sed 's/\\.rev.1.bt2\$//')
    echo "Using index path: \${indexPath}"
    bowtie2 -x \${indexPath} \
            -1 ${trimmomatic_out}/output_1P.fq \
            -2 ${trimmomatic_out}/output_2P.fq \
            -S "bowtie2_out/out.sam"
    """
}

// Define the process for seqtk
process seqtk {
    input:
    path bowtie2_out
    path trimmomatic_out

    output:
    path "${params.outputDir9}/*", emit: seqtk_out

    script:
    """
    mkdir -p "${params.outputDir9}"
    seqtk seq -A  "${trimmomatic_out}/output_1P.fq" > "${params.outputDir9}/paired_forward_1"
    seqtk seq -A  "${trimmomatic_out}/output_2P.fq" > "${params.outputDir9}/paired_forward_2"
    """
}

//Define process for pad read
process pad_read {
    input:
    tuple path(forward_read), path(reverse_read)
    path pad_read_path

    output:
    path "${params.outputDir9}/padded_out1.fa", emit: padded_file1
    path "${params.outputDir9}/padded_out2.fa", emit: padded_file2

    script:
    """
    mkdir -p "${params.outputDir9}"
    python "${params.pad_read_path}" "$forward_read" "${params.outputDir9}/padded_out1.fa" 150
    python "${params.pad_read_path}" "$reverse_read" "${params.outputDir9}/padded_out2.fa" 150
    """
}

//Define process for AlignGraph
process AlignGraph {
    input:
    path plentyofbugs_results
    path reference_genome
    path padded_file1
    path padded_file2
    path contigs_file
    val distancelow
    val distancehigh

    output:
    path "reference_based_assembly/sample_extendedcontig.fasta", emit: extended_contigs
    path "reference_based_assembly/sample_remainingcontig.fasta", emit: remaining_contigs

    script:
    """
    mkdir -p "reference_based_assembly"
    echo "Contig file path: $contigs_file"
    best_genome=\$(awk '{gsub(/.*\\//, ""); print \$1}' "best_reference")
    genome_file="${params.reference_genome}/\${best_genome}" 
    AlignGraph --read1 "$padded_file1" --read2 "$padded_file2" --contig "$contigs_file" --genome "\$genome_file"  --distanceLow "$distancelow"  --distanceHigh "$distancehigh"  --extendedContig "reference_based_assembly/sample_extendedcontig.fasta" --remainingContig "reference_based_assembly/sample_remainingcontig.fasta"
    """
}


// Define the process for quast2 (final Quast)
process quast2 {

    input:
    path extended_contigs

    output:
    path "${params.outputDir9}/quast_output"

    script:
    """
    mkdir -p ${params.outputDir9}/quast_output
    quast.py -o "${params.outputDir9}/quast_output" "$extended_contigs"
    """
}

// Define the process for BUSCO
process busco {

    input:
    path AlignGraph_results
    path extended_contigs
    path quast2_results

    output:
    path "${params.outputDir10}/busco_output"

    script:
    """
    mkdir -p ${params.outputDir10}
    busco -i "$extended_contigs" -o "${params.outputDir10}/busco_output" -l bacteria_odb10  --mode genome
    """
}

 workflow {

    // Run fastqc  and get the output
       def fastqc_result = fastqc(params.forward_read, params.reverse_read)

    // Run MultiQC on the FastQC output directory
       multiqc(fastqc_result) // Pass the FastQC output directory to MultiQC

    // Run trimmomatic on multiqc output dir
      def trimmomatic_results = trimmomatic(params.forward_read, params.reverse_read, params.thread, params.phred, params.PATH_TO_ADAPTER_CONTAM_FILE, params.leading, params.trailing, params.slidingwindow, params.minlength)

   // Run flash on trimmomatic dir
      def flash_results = flash(trimmomatic_results, params.max_overlap)

   // Run unicycler
      def unicycler_results = unicycler(trimmomatic_results, flash_results)

   // Run quast
      def quast_result = quast(unicycler_results)

   // Run plentyofbugs
      def (plentyofbugs_results, contigs_file) = plentyofbugs(quast_result, params.assembler, params.reference_genome, trimmomatic_results)

   // Run bowtie2_build process
      def bowtie2_build_results = bowtie2_build(plentyofbugs_results, params.reference_genome)

   // Run bowtie2 process using the emitted bowtie2_index
      def bowtie2_align_result = bowtie2(bowtie2_build_results.bowtie2_index, trimmomatic_results)

   // Run seqtk
      def seqtk_results = seqtk(bowtie2_align_result, trimmomatic_results)

      def forward_read = file("${seqtk_results.seqtk_out}/paired_forward_1")
      def reverse_read = file("${seqtk_results.seqtk_out}/paired_forward_2")

   // Run the pad_read process
      def pad_read_results =  pad_read(seqtk_results, params.pad_read_path)

   // Run AlignGraph
     def AlignGraph_results = AlignGraph(plentyofbugs_results, params.reference_genome, pad_read_results.padded_file1, pad_read_results.padded_file2, contigs_file, params.distancelow, params.distancehigh)

  // Run quast2 for AlignGraph
     def quast2_results = quast2(AlignGraph_results.extended_contigs)

  //Run busco
    busco(AlignGraph_results, quast2_results)
}
