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
    rm -rf ${params.outputDir7}
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
    path output:
    path"bowtie2_out/reference_index*", emit: bowtie2_index 

    script:
    """
    mkdir -p ${params.outputDir8}
    best_genome=\$(awk '{gsub(/.*\\//, ""); print \$1}' "best_reference")
    genome_file="${params.reference_genome}/\${best_genome}" 
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
    best_genome=\$(awk '{gsub(/.*\\//, ""); print \$1}' "best_reference")
    genome_file="${params.reference_genome}/\${best_genome}" 
    AlignGraph --read1 ${padded_file1} --read2 ${padded_file2} \
        --contig ${assembly_file} \
        --genome "\${genome_file}" \
        --distanceLow ${distancelow} \
        --distanceHigh ${distancehigh} \
        --extendedContig ${params.outputDir9}/sample_extendedcontig.fasta \
        --remainingContig ${params.outputDir9}/sample_remainingcontig.fasta
    """
}

// Define the process for quast2 (final Quast)
process quast2 {
    input:
    path extended_contigs

    output:
    path "${params.outputDir9}/quast_output", emit: quast2_out

    script:
    """
    mkdir -p ${params.outputDir9}/quast_output
    quast.py -o ${params.outputDir9}/quast_output ${extended_contigs}
    """
}

// Define the process for BUSCO
process busco {
    input:
    path extended_contigs

    output:
    path "${params.outputDir10}", emit: busco_out

    script:
    """
    mkdir -p ${params.outputDir10}
    busco -i ${extended_contigs} -o ${params.outputDir10} -l bacteria_odb10 --mode genome
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
    quast2_ch = quast2(aligngraph_ch.extended_contigs)
    
    // Run busco
    busco_ch = busco(aligngraph_ch.extended_contigs)
}
