#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ============================================================
// Parameters
// ============================================================
params.forward_read  = null
params.reverse_read  = null
params.thread        = 4
params.phred         = 33
params.PATH_TO_ADAPTER_CONTAM_FILE = null
params.leading       = 3
params.trailing      = 3
params.slidingwindow = '4:15'
params.minlength     = 36
params.max_overlap   = 10
params.assembler     = 'skesa'
params.reference_genome = null
params.pad_read_path = null
params.pad_read_len  = 150
params.distancelow   = 200
params.distancehigh  = 1000
params.busco_lineage = 'bacteria_odb10'
params.align_threshold = 75
params.outdir        = "results"

// ============================================================
// PHASE 1
// ============================================================

process FASTQC {
    publishDir "${params.outdir}/fastqc_out", mode: 'copy'
    input:
    path fwd
    path rev
    output:
    path "fastqc_out", emit: out_dir
    script:
    """
    mkdir -p fastqc_out
    fastqc -o fastqc_out -f fastq -t ${params.thread} ${fwd} ${rev}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc_out", mode: 'copy'
    input:
    path fastqc_out
    output:
    path "multiqc_out", emit: out_dir
    script:
    """
    mkdir -p multiqc_out
    multiqc ${fastqc_out} -o multiqc_out
    """
}

process TRIMMOMATIC {
    publishDir "${params.outdir}/trimmomatic_out", mode: 'copy'
    input:
    path fwd
    path rev
    val  adapter_path
    output:
    path "trimmomatic_out", emit: trim_out
    script:
    def adapter_clip = adapter_path ? "ILLUMINACLIP:${adapter_path}:2:30:10" : ""
    """
    mkdir -p trimmomatic_out
    trimmomatic PE -threads ${params.thread} -phred${params.phred} \\
        ${fwd} ${rev} \\
        trimmomatic_out/output_1P.fq trimmomatic_out/output_1U.fq \\
        trimmomatic_out/output_2P.fq trimmomatic_out/output_2U.fq \\
        ${adapter_clip} \\
        LEADING:${params.leading} \\
        TRAILING:${params.trailing} \\
        SLIDINGWINDOW:${params.slidingwindow} \\
        MINLEN:${params.minlength}
    """
}

process FLASH {
    publishDir "${params.outdir}/flash_out", mode: 'copy'
    input:
    path trim_out
    output:
    path "flash_out", emit: flash_out
    script:
    """
    mkdir -p flash_out
    flash --max-overlap ${params.max_overlap} --threads ${params.thread} \\
        --output-prefix out \\
        --output-directory flash_out \\
        ${trim_out}/output_1P.fq \\
        ${trim_out}/output_2P.fq
    """
}

process UNICYCLER {
    publishDir "${params.outdir}/unicycler_out", mode: 'copy'
    input:
    path flash_out
    output:
    path "unicycler_out", emit: unicycler_out
    script:
    """
    mkdir -p unicycler_out
    unicycler \\
        -1 ${flash_out}/out.notCombined_1.fastq \\
        -2 ${flash_out}/out.notCombined_2.fastq \\
        -s ${flash_out}/out.extendedFrags.fastq \\
        -t ${params.thread} \\
        -o unicycler_out/assembly
    """
}

process QUAST_DENOVO {
    publishDir "${params.outdir}/quast_out", mode: 'copy'
    input:
    path unicycler_out
    output:
    path "quast_out", emit: quast_out
    script:
    """
    quast -o quast_out -t ${params.thread} --min-contig 0 \\
        ${unicycler_out}/assembly/assembly.fasta
    """
}

// ============================================================
// PHASE 2
// ============================================================

process PLENTYOFBUGS {
    publishDir "${params.outdir}/plentyofbugs_out", mode: 'copy'
    input:
    path trim_out
    val  assembler
    val  reference_dir
    output:
    path "plentyofbugs_out", emit: pob_out
    script:
    """
    plentyofbugs --assembler ${assembler} \\
        -f ${trim_out}/output_1P.fq \\
        -r ${trim_out}/output_2P.fq \\
        -g ${reference_dir} \\
        -o plentyofbugs_out
    """
}

process BOWTIE2_BUILD {
    publishDir "${params.outdir}/bowtie2_out", mode: 'copy'
    input:
    path pob_out
    val  reference_dir
    output:
    path "bt2_index",           emit: index
    path "best_reference.txt",  emit: best_ref
    script:
    """
    mkdir -p bt2_index
    best_ref=\$(awk 'NR==1 {print \$1}' ${pob_out}/best_reference)
    if [ ! -f "\$best_ref" ]; then
        best_ref="${reference_dir}/\$(basename "\$best_ref")"
    fi
    echo "\$best_ref" > best_reference.txt
    bowtie2-build -f "\$best_ref" bt2_index/reference_index
    """
}

process BOWTIE2_ALIGN {
    publishDir "${params.outdir}/bowtie2_out", mode: 'copy'
    input:
    path bt2_index
    path trim_out
    output:
    path "out.sam",            emit: sam
    path "alignment_rate.txt", emit: rate
    script:
    """
    bowtie2 -p ${params.thread} \\
        -x ${bt2_index}/reference_index \\
        -1 ${trim_out}/output_1P.fq \\
        -2 ${trim_out}/output_2P.fq \\
        -S out.sam

    total=\$(samtools view -S -c out.sam)
    aligned=\$(samtools view -S -F 4 -c out.sam)
    rate=\$(awk -v a="\$aligned" -v t="\$total" \\
        'BEGIN {if(t>0) printf "%.2f", (a/t)*100; else print "0.00"}')
    {
        echo "total=\$total"
        echo "aligned=\$aligned"
        echo "rate=\$rate"
    } > alignment_rate.txt
    echo "Alignment rate: \$rate% (threshold: ${params.align_threshold}%)"
    """
}

process CHECK_ALIGNMENT {
    publishDir "${params.outdir}/bowtie2_out", mode: 'copy'
    input:
    path rate_file
    output:
    path "alignment_check.txt", emit: check
    script:
    """
    rate=\$(grep '^rate=' ${rate_file} | cut -d= -f2)
    if awk -v r="\$rate" -v t="${params.align_threshold}" 'BEGIN {exit !(r >= t)}'; then
        echo "PASS: rate \$rate% >= threshold ${params.align_threshold}%" > alignment_check.txt
    else
        echo "WARN: rate \$rate% < threshold ${params.align_threshold}% - continuing anyway" > alignment_check.txt
    fi
    cat alignment_check.txt
    """
}

process SEQTK {
    publishDir "${params.outdir}/reference_based_assembly", mode: 'copy'
    input:
    path trim_out
    output:
    path "seqtk_out", emit: seqtk_out
    script:
    """
    mkdir -p seqtk_out
    seqtk seq -A ${trim_out}/output_1P.fq > seqtk_out/output_1P.fa
    seqtk seq -A ${trim_out}/output_2P.fq > seqtk_out/output_2P.fa
    """
}

process PAD_READ {
    publishDir "${params.outdir}/reference_based_assembly", mode: 'copy'
    input:
    path seqtk_out
    path pad_script
    output:
    path "padded_out1.fa", emit: padded1
    path "padded_out2.fa", emit: padded2
    script:
    """
    python ${pad_script} ${seqtk_out}/output_1P.fa padded_out1.fa ${params.pad_read_len}
    python ${pad_script} ${seqtk_out}/output_2P.fa padded_out2.fa ${params.pad_read_len}
    """
}

process ALIGNGRAPH {
    publishDir "${params.outdir}/reference_based_assembly", mode: 'copy'
    input:
    path pob_out
    path unicycler_out
    path padded1
    path padded2
    val  reference_dir
    output:
    path "sample_extendedcontig.fasta",  emit: extended
    path "sample_remainingcontig.fasta", emit: remaining
    script:
    """
    best_ref=\$(awk 'NR==1 {print \$1}' ${pob_out}/best_reference)
    if [ ! -f "\$best_ref" ]; then
        best_ref="${reference_dir}/\$(basename "\$best_ref")"
    fi
    AlignGraph \\
        --read1 ${padded1} \\
        --read2 ${padded2} \\
        --contig ${unicycler_out}/assembly/assembly.fasta \\
        --genome "\$best_ref" \\
        --distanceLow ${params.distancelow} \\
        --distanceHigh ${params.distancehigh} \\
        --extendedContig sample_extendedcontig.fasta \\
        --remainingContig sample_remainingcontig.fasta
    """
}

// ============================================================
// MERGE: extended + remaining → one complete assembly file
// ============================================================
process MERGE_ALIGNGRAPH {
    publishDir "${params.outdir}/reference_based_assembly", mode: 'copy'
    input:
    path extended
    path remaining
    path unicycler_out
    output:
    path "complete_assembly.fasta", emit: complete
    path "merge_report.txt",        emit: report
    script:
    """
    # Merge extended + remaining into one complete assembly
    if [ -s ${extended} ] && [ -s ${remaining} ]; then
        cat ${extended} ${remaining} > complete_assembly.fasta
        echo "STATUS=OK"                                > merge_report.txt
        echo "MODE=extended+remaining"                  >> merge_report.txt
        echo "Extended contigs: \$(grep -c '^>' ${extended})"   >> merge_report.txt
        echo "Remaining contigs: \$(grep -c '^>' ${remaining})"  >> merge_report.txt
    elif [ -s ${extended} ]; then
        cp ${extended} complete_assembly.fasta
        echo "STATUS=EXTENDED_ONLY"                     > merge_report.txt
        echo "MODE=extended_only"                       >> merge_report.txt
    elif [ -s ${remaining} ]; then
        cp ${remaining} complete_assembly.fasta
        echo "STATUS=REMAINING_ONLY"                    > merge_report.txt
        echo "MODE=remaining_only"                      >> merge_report.txt
    else
        cp ${unicycler_out}/assembly/assembly.fasta complete_assembly.fasta
        echo "STATUS=DENOVO_FALLBACK"                   > merge_report.txt
        echo "MODE=de_novo_fallback"                    >> merge_report.txt
    fi

    echo "Total contigs: \$(grep -c '^>' complete_assembly.fasta)" >> merge_report.txt
    echo "Total bp: \$(awk '/^>/{next}{s+=length}END{print s}' complete_assembly.fasta)" >> merge_report.txt
    cat merge_report.txt
    """
}

process QUAST_REFGUIDED {
    publishDir "${params.outdir}/quast2_out", mode: 'copy'
    input:
    path complete_asm
    output:
    path "quast_out", emit: quast_out
    script:
    """
    quast -o quast_out -t ${params.thread} --min-contig 0 ${complete_asm}
    """
}

process BUSCO {
    publishDir "${params.outdir}/busco_out", mode: 'copy'
    input:
    path complete_asm
    output:
    path "busco_out", emit: busco_out
    script:
    """
    mkdir -p busco_out
    busco -i ${complete_asm} -m genome -l ${params.busco_lineage} \\
        -c ${params.thread} -o busco_run --out_path busco_out
    """
}

// ============================================================
// WORKFLOW
// ============================================================

workflow {
    if (!params.forward_read)     error("--forward_read is required")
    if (!params.reverse_read)     error("--reverse_read is required")
    if (!params.reference_genome) error("--reference_genome is required")
    if (!params.pad_read_path)    error("--pad_read_path is required")

    fwd        = file(params.forward_read,  checkIfExists: true)
    rev        = file(params.reverse_read,  checkIfExists: true)
    pad_script = file(params.pad_read_path, checkIfExists: true)

    // Phase 1
    FASTQC(fwd, rev)
    MULTIQC(FASTQC.out.out_dir)
    TRIMMOMATIC(fwd, rev, params.PATH_TO_ADAPTER_CONTAM_FILE)
    FLASH(TRIMMOMATIC.out.trim_out)
    UNICYCLER(FLASH.out.flash_out)
    QUAST_DENOVO(UNICYCLER.out.unicycler_out)

    // Phase 2
    PLENTYOFBUGS(TRIMMOMATIC.out.trim_out, params.assembler, params.reference_genome)
    BOWTIE2_BUILD(PLENTYOFBUGS.out.pob_out, params.reference_genome)
    BOWTIE2_ALIGN(BOWTIE2_BUILD.out.index, TRIMMOMATIC.out.trim_out)
    CHECK_ALIGNMENT(BOWTIE2_ALIGN.out.rate)
    SEQTK(TRIMMOMATIC.out.trim_out)
    PAD_READ(SEQTK.out.seqtk_out, pad_script)
    ALIGNGRAPH(
        PLENTYOFBUGS.out.pob_out,
        UNICYCLER.out.unicycler_out,
        PAD_READ.out.padded1,
        PAD_READ.out.padded2,
        params.reference_genome
    )

    // Merge extended + remaining → complete assembly
    MERGE_ALIGNGRAPH(
        ALIGNGRAPH.out.extended,
        ALIGNGRAPH.out.remaining,
        UNICYCLER.out.unicycler_out
    )

    // QUAST + BUSCO on merged complete assembly
    QUAST_REFGUIDED(MERGE_ALIGNGRAPH.out.complete)
    BUSCO(MERGE_ALIGNGRAPH.out.complete)
}
