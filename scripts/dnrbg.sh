#!/usr/bin/env bash
#
# dnrbg_1.sh — integrative genome assembly pipeline (with resume + fallback)
#
# Phase 1: FastQC -> MultiQC -> Trimmomatic -> FLASH -> Unicycler -> QUAST
# Phase 2: plentyofbugs -> Bowtie2 -> AlignGraph -> QUAST -> BUSCO
#
# Resume:
#   ./dnrbg_1.sh --resume <WORK_DIR> [other flags]
#
# AlignGraph fallback:
#   If sample_extendedcontig.fasta is empty (0 bytes), the pipeline
#   automatically uses sample_remainingcontig.fasta for QUAST and BUSCO.
#
set -euo pipefail

###############################################################################
# Colours
###############################################################################
if [[ -t 1 ]]; then
    green=$'\033[0;32m'
    red=$'\033[0;31m'
    yellow=$'\033[1;33m'
    blue=$'\033[0;34m'
    reset=$'\033[0m'
else
    green=""; red=""; yellow=""; blue=""; reset=""
fi

log()   { echo -e "${green}[INFO]${reset}  $*"; }
warn()  { echo -e "${yellow}[WARN]${reset}  $*" >&2; }
error() { echo -e "${red}[ERROR]${reset} $*" >&2; }
die()   { error "$*"; exit 1; }

###############################################################################
# Defaults
###############################################################################
THREADS=4
PHRED=33
ADAPTER_FILE=""
LEADING=3
TRAILING=3
SLIDINGWINDOW="4:15"
MINLEN=36

MAX_OVERLAP=10
ASSEMBLER="skesa"
BUSCO_LINEAGE="bacteria_odb10"
PAD_SCRIPT=""
PAD_READ_LEN=150
DIST_LOW=200
DIST_HIGH=1000
REFERENCE_DIR=""
SKIP_PHASE2=0
ALIGN_THRESHOLD=75
FRESH=0
RESUME_DIR=""

FORWARD_READ=""
REVERSE_READ=""

# Set inside main() AFTER parse_args so --resume takes effect
WORK_DIR=""
CHECKPOINT_DIR=""
LOG_FILE=""
REPORT_FILE=""
CONFIG_FILE=""
CANDIDATE_REFS=()

# Set during execution
ALIGNGRAPH_STATUS="not run"
ALIGNGRAPH_USED_FILE=""

###############################################################################
# Help
###############################################################################
usage() {
    cat <<EOF
${green}Integrative genome assembly pipeline${reset}

${yellow}Usage:${reset}
  $(basename "$0") -1 <R1.fastq> -2 <R2.fastq> -g <REF_DIR> [options]
  $(basename "$0") --resume <WORK_DIR> [options]
  $(basename "$0") --fresh -1 <R1.fastq> -2 <R2.fastq> [options]

${yellow}Resume:${reset}
  --resume DIR   Resume a previous run (skip completed steps)
  --fresh        Force a new run

${yellow}Required:${reset}
  -1 FILE        Path to forward reads (FASTQ)
  -2 FILE        Path to reverse reads (FASTQ)
  -g DIR         Directory of candidate reference FASTA files

${yellow}Trimmomatic:${reset}
  -T INT         Threads                       (default: $THREADS)
  -P INT         Phred (33|64)                 (default: $PHRED)
  -A FILE        Adapter FASTA
  -L INT         LEADING                       (default: $LEADING)
  -G INT         TRAILING                      (default: $TRAILING)
  -W STR         SLIDINGWINDOW                 (default: $SLIDINGWINDOW)
  -M INT         MINLEN                        (default: $MINLEN)

${yellow}Other tools:${reset}
  -o INT         FLASH max-overlap             (default: $MAX_OVERLAP)
  -X STR         Assembler (skesa|spades)      (default: $ASSEMBLER)
  -b STR         BUSCO lineage                 (default: $BUSCO_LINEAGE)
  -C FILE        Path to pad_reads_1.py
  -R INT         Read length for padding       (default: $PAD_READ_LEN)
  -d INT         AlignGraph --distanceLow      (default: $DIST_LOW)
  -D INT         AlignGraph --distanceHigh     (default: $DIST_HIGH)
  -m INT         Alignment-rate threshold %    (default: $ALIGN_THRESHOLD)
  -S             Skip Phase 2
  -h             Show this help
EOF
}

###############################################################################
# Argument parsing
###############################################################################
parse_args() {
    local args=("$@")
    local i=0

    while [[ $i -lt ${#args[@]} ]]; do
        case "${args[$i]}" in
            --resume)
                RESUME_DIR="${args[$((i+1))]:-}"
                [[ -n "$RESUME_DIR" ]] || die "--resume requires a directory."
                [[ -d "$RESUME_DIR" ]] || die "Resume dir not found: $RESUME_DIR"
                i=$((i+2))
                ;;
            --fresh)
                FRESH=1
                i=$((i+1))
                ;;
            *) break ;;
        esac
    done

    set -- "${args[@]:$i}"
    OPTIND=1
    while getopts ":1:2:T:P:A:L:G:W:M:o:X:b:C:R:d:D:g:m:Sh" opt; do
        case "$opt" in
            1) FORWARD_READ=$OPTARG ;;
            2) REVERSE_READ=$OPTARG ;;
            T) THREADS=$OPTARG ;;
            P) PHRED=$OPTARG ;;
            A) ADAPTER_FILE=$OPTARG ;;
            L) LEADING=$OPTARG ;;
            G) TRAILING=$OPTARG ;;
            W) SLIDINGWINDOW=$OPTARG ;;
            M) MINLEN=$OPTARG ;;
            o) MAX_OVERLAP=$OPTARG ;;
            X) ASSEMBLER=$OPTARG ;;
            b) BUSCO_LINEAGE=$OPTARG ;;
            C) PAD_SCRIPT=$OPTARG ;;
            R) PAD_READ_LEN=$OPTARG ;;
            d) DIST_LOW=$OPTARG ;;
            D) DIST_HIGH=$OPTARG ;;
            g) REFERENCE_DIR=$OPTARG ;;
            m) ALIGN_THRESHOLD=$OPTARG ;;
            S) SKIP_PHASE2=1 ;;
            h) usage; exit 0 ;;
            :) die "Option -$OPTARG requires an argument." ;;
            \?) die "Unknown option: -$OPTARG" ;;
        esac
    done
}

###############################################################################
# Directory / checkpoint helpers
###############################################################################
ensure_dir() {
    local d=$1
    [[ -d "$d" ]] || mkdir -p "$d"
    echo "$d"
}

step_done() {
    [[ -f "$CHECKPOINT_DIR/$1.done" ]]
}

mark_done() {
    ensure_dir "$CHECKPOINT_DIR" >/dev/null
    date > "$CHECKPOINT_DIR/$1.done"
    log "✓ Checkpoint: $1"
}

run_step() {
    local name=$1
    shift
    if step_done "$name"; then
        log "Skipping $name (already completed)"
        return 0
    fi
    "$@"
    mark_done "$name"
}

###############################################################################
# Save / load config
###############################################################################
save_config() {
    ensure_dir "$WORK_DIR" >/dev/null
    cat > "$CONFIG_FILE" <<EOF
# Pipeline configuration — auto-generated
FORWARD_READ="$FORWARD_READ"
REVERSE_READ="$REVERSE_READ"
REFERENCE_DIR="$REFERENCE_DIR"
THREADS="$THREADS"
PHRED="$PHRED"
ADAPTER_FILE="$ADAPTER_FILE"
LEADING="$LEADING"
TRAILING="$TRAILING"
SLIDINGWINDOW="$SLIDINGWINDOW"
MINLEN="$MINLEN"
MAX_OVERLAP="$MAX_OVERLAP"
ASSEMBLER="$ASSEMBLER"
BUSCO_LINEAGE="$BUSCO_LINEAGE"
PAD_SCRIPT="$PAD_SCRIPT"
PAD_READ_LEN="$PAD_READ_LEN"
DIST_LOW="$DIST_LOW"
DIST_HIGH="$DIST_HIGH"
ALIGN_THRESHOLD="$ALIGN_THRESHOLD"
SKIP_PHASE2="$SKIP_PHASE2"
EOF
}

load_config() {
    local cfg="$WORK_DIR/.pipeline_config"
    [[ -f "$cfg" ]] || return 0
    # shellcheck disable=SC1090
    source "$cfg"
    log "Loaded config from previous run: $cfg"
}

###############################################################################
# Input validation
###############################################################################
validate_inputs() {
    [[ -n "$FORWARD_READ" ]] || { usage; die "Forward read (-1) is required."; }
    [[ -n "$REVERSE_READ" ]] || { usage; die "Reverse read (-2) is required."; }
    [[ -f "$FORWARD_READ" ]] || die "Forward read not found: $FORWARD_READ"
    [[ -f "$REVERSE_READ" ]] || die "Reverse read not found: $REVERSE_READ"

    [[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "Threads must be a positive integer (got: $THREADS)."
    [[ "$PHRED" == "33" || "$PHRED" == "64" ]] || die "Phred must be 33 or 64 (got: $PHRED)."
    [[ "$ASSEMBLER" == "skesa" || "$ASSEMBLER" == "spades" ]] || die "Assembler must be 'skesa' or 'spades'."
    [[ -z "$ADAPTER_FILE" || -f "$ADAPTER_FILE" ]] || die "Adapter file not found: $ADAPTER_FILE"

    if [[ "$SKIP_PHASE2" -eq 0 ]]; then
        [[ -n "$REFERENCE_DIR" ]] || die "Reference directory (-g) required for Phase 2 (or pass -S)."
        [[ -d "$REFERENCE_DIR" ]] || die "Reference directory not found: $REFERENCE_DIR"

        shopt -s nullglob
        CANDIDATE_REFS=(
            "$REFERENCE_DIR"/*.fa
            "$REFERENCE_DIR"/*.fna
            "$REFERENCE_DIR"/*.fasta
            "$REFERENCE_DIR"/*.fa.gz
            "$REFERENCE_DIR"/*.fna.gz
            "$REFERENCE_DIR"/*.fasta.gz
        )
        shopt -u nullglob

        [[ ${#CANDIDATE_REFS[@]} -gt 0 ]] || \
            die "No FASTA files found in: $REFERENCE_DIR"

        [[ -n "$PAD_SCRIPT" ]] || die "pad_reads_1.py (-C) required for Phase 2 (or pass -S)."
        [[ -f "$PAD_SCRIPT" ]] || die "pad_reads_1.py not found: $PAD_SCRIPT"
    fi
}

###############################################################################
# Tool resolution
###############################################################################
declare -A TOOL_PATH

resolve_tool() {
    local name=$1
    if [[ -n "${TOOL_PATH[$name]:-}" ]]; then
        echo "${TOOL_PATH[$name]}"
        return 0
    fi
    if command -v "$name" &>/dev/null; then
        TOOL_PATH[$name]=$(command -v "$name")
    else
        warn "$name not found on PATH."
        local p
        read -rp "Absolute path to $name: " p
        [[ -x "$p" ]] || die "Not executable: $p"
        TOOL_PATH[$name]=$p
    fi
    echo "${TOOL_PATH[$name]}"
}

###############################################################################
# Phase 1.1 — FastQC
###############################################################################
run_fastqc() {
    log "=== Phase 1.1: FastQC ==="
    local fastqc_bin; fastqc_bin=$(resolve_tool fastqc)
    local out; out=$(ensure_dir "$WORK_DIR/fastqc_out")

    export LC_ALL=C
    "$fastqc_bin" --outdir "$out" -f fastq -t "$THREADS" \
        "$FORWARD_READ" "$REVERSE_READ"

    log "FastQC complete -> $out"
}

###############################################################################
# Phase 1.2 — MultiQC
###############################################################################
run_multiqc() {
    log "=== Phase 1.2: MultiQC ==="
    local multiqc_bin; multiqc_bin=$(resolve_tool multiqc)
    local out; out=$(ensure_dir "$WORK_DIR/multiqc_out")

    "$multiqc_bin" "$WORK_DIR/fastqc_out" -o "$out"

    log "MultiQC complete -> $out"
}

###############################################################################
# Phase 1.3 — Trimmomatic (wrapper on PATH)
###############################################################################
run_trimmomatic() {
    log "=== Phase 1.3: Trimmomatic ==="
    command -v java &>/dev/null || die "Java is required for Trimmomatic."

    local out; out=$(ensure_dir "$WORK_DIR/trimmomatic_out")

    local -a adapter_args=()
    if [[ -n "$ADAPTER_FILE" ]]; then
        adapter_args=(ILLUMINACLIP:"$ADAPTER_FILE":2:30:10)
    fi

    if command -v trimmomatic &>/dev/null; then
        trimmomatic PE \
            -threads "$THREADS" \
            -phred"$PHRED" \
            "$FORWARD_READ" "$REVERSE_READ" \
            "$out/output_1P.fq" "$out/output_1U.fq" \
            "$out/output_2P.fq" "$out/output_2U.fq" \
            "${adapter_args[@]}" \
            LEADING:"$LEADING" \
            TRAILING:"$TRAILING" \
            SLIDINGWINDOW:"$SLIDINGWINDOW" \
            MINLEN:"$MINLEN"
    else
        local trim_jar
        read -rp "Absolute path to trimmomatic.jar: " trim_jar
        [[ -f "$trim_jar" ]] || die "trimmomatic.jar not found: $trim_jar"
        java -jar "$trim_jar" PE \
            -threads "$THREADS" \
            -phred"$PHRED" \
            "$FORWARD_READ" "$REVERSE_READ" \
            "$out/output_1P.fq" "$out/output_1U.fq" \
            "$out/output_2P.fq" "$out/output_2U.fq" \
            "${adapter_args[@]}" \
            LEADING:"$LEADING" \
            TRAILING:"$TRAILING" \
            SLIDINGWINDOW:"$SLIDINGWINDOW" \
            MINLEN:"$MINLEN"
    fi

    log "Trimmomatic complete -> $out"
}

###############################################################################
# Phase 1.4 — FLASH
###############################################################################
run_flash() {
    log "=== Phase 1.4: FLASH ==="
    local flash_bin; flash_bin=$(resolve_tool flash)
    local out; out=$(ensure_dir "$WORK_DIR/flash_out")

    "$flash_bin" \
        --max-overlap "$MAX_OVERLAP" \
        --threads "$THREADS" \
        --output-prefix out \
        --output-directory "$out" \
        "$WORK_DIR/trimmomatic_out/output_1P.fq" \
        "$WORK_DIR/trimmomatic_out/output_2P.fq"

    log "FLASH complete -> $out"
}

###############################################################################
# Phase 1.5 — Unicycler
###############################################################################
run_unicycler() {
    log "=== Phase 1.5: Unicycler ==="
    local unicycler_bin; unicycler_bin=$(resolve_tool unicycler)
    resolve_tool spades >/dev/null

    local out; out=$(ensure_dir "$WORK_DIR/unicycler_out")

    "$unicycler_bin" \
        -1 "$WORK_DIR/trimmomatic_out/output_1P.fq" \
        -2 "$WORK_DIR/trimmomatic_out/output_2P.fq" \
        -s "$WORK_DIR/flash_out/out.extendedFrags.fastq" \
        -t "$THREADS" \
        -o "$out/assembly"

    log "Unicycler complete -> $out/assembly"
}

###############################################################################
# QUAST (reusable)
###############################################################################
run_quast() {
    local label=$1
    local assembly=$2
    local outdir=$3

    log "=== QUAST ($label) ==="
    local quast_bin; quast_bin=$(resolve_tool quast)

    [[ -s "$assembly" ]] || die "Assembly FASTA not found or empty: $assembly"
    ensure_dir "$outdir" >/dev/null

    "$quast_bin" -o "$outdir/quast_output" -t "$THREADS" --min-contig 0 "$assembly"

    log "QUAST ($label) complete -> $outdir/quast_output"
}

###############################################################################
# Phase 2.1 — plentyofbugs
###############################################################################
run_plentyofbugs() {
    log "=== Phase 2.1: plentyofbugs ==="
    log "Candidate references in $REFERENCE_DIR: ${#CANDIDATE_REFS[@]}"

    local pob_bin; pob_bin=$(resolve_tool plentyofbugs)
    resolve_tool mash >/dev/null
    resolve_tool seqtk >/dev/null
    resolve_tool "$ASSEMBLER" >/dev/null

    local out="$WORK_DIR/plentyofbugs_out"
    if [[ -d "$out" ]]; then
        log "Removing stale plentyofbugs_out/ so plentyofbugs can run fresh."
        rm -rf "$out"
    fi

    "$pob_bin" \
        --assembler "$ASSEMBLER" \
        -f "$WORK_DIR/trimmomatic_out/output_1P.fq" \
        -r "$WORK_DIR/trimmomatic_out/output_2P.fq" \
        -g "$REFERENCE_DIR" \
        -o "$out"

    log "plentyofbugs complete -> $out"
}

###############################################################################
# Extract chosen reference path
###############################################################################
get_best_reference() {
    local best_ref_file="$WORK_DIR/plentyofbugs_out/best_reference"
    [[ -f "$best_ref_file" ]] || die "best_reference not found: $best_ref_file"
    local ref_path
    ref_path=$(awk 'NR==1 {print $1}' "$best_ref_file")
    [[ -f "$ref_path" ]] || die "Reference listed in best_reference not found: $ref_path"
    echo "$ref_path"
}

###############################################################################
# Alignment rate (correct bitmask handling)
###############################################################################
ALIGN_TOTAL="NA"
ALIGN_ALIGNED="NA"
ALIGN_RATE="NA"

compute_alignment_rate() {
    local sam=$1
    local total aligned

    if command -v samtools &>/dev/null; then
        total=$(samtools view -c "$sam")
        aligned=$(samtools view -F 4 -c "$sam")
    elif awk 'BEGIN{exit !(and(1,1)==1)}' 2>/dev/null; then
        total=$(grep -vc '^@' "$sam" || true)
        aligned=$(awk '!/^@/ && and($2,4)==0' "$sam" | wc -l)
    else
        die "Need samtools or GNU awk for alignment-rate calculation."
    fi

    ALIGN_TOTAL=${total:-0}
    ALIGN_ALIGNED=${aligned:-0}

    awk -v a="$ALIGN_ALIGNED" -v t="$ALIGN_TOTAL" \
        'BEGIN { if (t > 0) printf "%.2f", (a/t)*100; else print "0.00" }'
}

###############################################################################
# Phase 2.2 — Bowtie2 + threshold (asks user in BOTH cases)
###############################################################################
run_bowtie2() {
    log "=== Phase 2.2: Bowtie2 alignment ==="
    local bt2_bin; bt2_bin=$(resolve_tool bowtie2)
    local bt2_build; bt2_build=$(resolve_tool bowtie2-build)

    local ref_path; ref_path=$(get_best_reference)
    log "Best reference from plentyofbugs: $ref_path"

    local out; out=$(ensure_dir "$WORK_DIR/bowtie2_out")
    local idx="$out/reference_index"

    "$bt2_build" -f "$ref_path" "$idx"

    "$bt2_bin" \
        -p "$THREADS" \
        -x "$idx" \
        -1 "$WORK_DIR/trimmomatic_out/output_1P.fq" \
        -2 "$WORK_DIR/trimmomatic_out/output_2P.fq" \
        -S "$out/out.sam"

    local rate
    rate=$(compute_alignment_rate "$out/out.sam")
    ALIGN_RATE="$rate"

    cat > "$out/alignment_report.txt" <<EOF
Bowtie2 Alignment Report
========================
SAM file                : $out/out.sam
Reference               : $ref_path
Total alignment records : $ALIGN_TOTAL
Aligned records         : $ALIGN_ALIGNED
Alignment rate          : ${rate}%
Threshold               : ${ALIGN_THRESHOLD}%
EOF

    log "Alignment rate: ${rate}% (threshold: ${ALIGN_THRESHOLD}%)"
    cat "$out/alignment_report.txt"

    if awk -v r="$rate" -v t="$ALIGN_THRESHOLD" 'BEGIN { exit !(r >= t) }'; then
        echo -e "${green}Alignment rate is ${rate}%, which is >= ${ALIGN_THRESHOLD}%.${reset}"
        echo -e "${yellow}Do you want to continue the pipeline? (yes/no)${reset}"
        read -rp "Please provide your answer: " continue_pipeline
        if [[ "$continue_pipeline" == "yes" || "$continue_pipeline" == "y" ]]; then
            log "Continuing pipeline..."
        else
            log "Pipeline stopped by user."
            exit 0
        fi
    else
        echo -e "${red}Alignment rate is ${rate}%, which is < ${ALIGN_THRESHOLD}%.${reset}"
        echo -e "${yellow}The alignment rate is below the recommended threshold.${reset}"
        echo -e "${yellow}Do you still want to continue the pipeline? (yes/no)${reset}"
        read -rp "Please provide your answer: " continue_pipeline
        if [[ "$continue_pipeline" == "yes" || "$continue_pipeline" == "y" ]]; then
            warn "Continuing pipeline despite alignment rate < ${ALIGN_THRESHOLD}%..."
        else
            log "Pipeline stopped by user."
            exit 0
        fi
    fi
}

###############################################################################
# Phase 2.3 — AlignGraph
###############################################################################
run_aligngraph() {
    log "=== Phase 2.3: AlignGraph ==="
    local seqtk_bin; seqtk_bin=$(resolve_tool seqtk)
    local ag_bin; ag_bin=$(resolve_tool AlignGraph)

    local ref_path; ref_path=$(get_best_reference)

    local contigs="$WORK_DIR/unicycler_out/assembly/assembly.fasta"
    [[ -s "$contigs" ]] || die "Unicycler assembly not found or empty: $contigs"

    local out; out=$(ensure_dir "$WORK_DIR/reference_based_assembly")

    "$seqtk_bin" seq -A "$WORK_DIR/trimmomatic_out/output_1P.fq" > "$out/output_1P.fa"
    "$seqtk_bin" seq -A "$WORK_DIR/trimmomatic_out/output_2P.fq" > "$out/output_2P.fa"

    python "$PAD_SCRIPT" "$out/output_1P.fa" "$out/padded_out1.fa" "$PAD_READ_LEN"
    python "$PAD_SCRIPT" "$out/output_2P.fa" "$out/padded_out2.fa" "$PAD_READ_LEN"

    (
        cd "$out"
        "$ag_bin" \
            --read1 padded_out1.fa \
            --read2 padded_out2.fa \
            --contig "$contigs" \
            --genome "$ref_path" \
            --distanceLow "$DIST_LOW" \
            --distanceHigh "$DIST_HIGH" \
            --extendedContig sample_extendedcontig.fasta \
            --remainingContig sample_remainingcontig.fasta
    )

    log "AlignGraph complete -> $out"

    # --- Fallback logic ---
    local ext="$out/sample_extendedcontig.fasta"
    local rem="$out/sample_remainingcontig.fasta"

    if [[ -s "$ext" ]]; then
        ALIGNGRAPH_STATUS="OK"
        ALIGNGRAPH_USED_FILE="$ext"
        log "AlignGraph produced extended contigs: $ext"
    elif [[ -s "$rem" ]]; then
        ALIGNGRAPH_STATUS="FALLBACK"
        ALIGNGRAPH_USED_FILE="$rem"
        warn "extendedContig.fasta is empty — falling back to remainingContig.fasta."
        warn "Downstream QUAST/BUSCO will use: $rem"
    else
        ALIGNGRAPH_STATUS="FAILED"
        ALIGNGRAPH_USED_FILE=""
        die "AlignGraph produced neither extended nor remaining contigs."
    fi
}

###############################################################################
# Phase 2.4 — BUSCO (uses ALIGNGRAPH_USED_FILE)
###############################################################################
run_busco() {
    log "=== Phase 2.4: BUSCO ==="
    local busco_bin; busco_bin=$(resolve_tool busco)

    local out; out=$(ensure_dir "$WORK_DIR/busco_out")

    # Choose input based on AlignGraph status
    local input=""
    if [[ -n "$ALIGNGRAPH_USED_FILE" && -s "$ALIGNGRAPH_USED_FILE" ]]; then
        input="$ALIGNGRAPH_USED_FILE"
    else
        # Fallback chain
        if [[ -s "$WORK_DIR/reference_based_assembly/sample_extendedcontig.fasta" ]]; then
            input="$WORK_DIR/reference_based_assembly/sample_extendedcontig.fasta"
        elif [[ -s "$WORK_DIR/reference_based_assembly/sample_remainingcontig.fasta" ]]; then
            input="$WORK_DIR/reference_based_assembly/sample_remainingcontig.fasta"
        else
            input="$WORK_DIR/unicycler_out/assembly/assembly.fasta"
        fi
    fi
    [[ -s "$input" ]] || die "No non-empty assembly FASTA found for BUSCO."

    log "BUSCO input: $input"

    (
        cd "$out"
        "$busco_bin" \
            -i "$input" \
            -m genome \
            -l "$BUSCO_LINEAGE" \
            -c "$THREADS" \
            -o busco_run
    )

    log "BUSCO complete -> $out/busco_run"
}

###############################################################################
# QUAST metric extractor
###############################################################################
quast_metric() {
    local report=$1
    local metric=$2
    [[ -f "$report" ]] || { echo "NA"; return; }
    awk -F'\t' -v m="$metric" '$1==m {print $2; exit}' "$report" 2>/dev/null | head -1
}

###############################################################################
# Summary report
###############################################################################
write_summary_report() {
    log "=== Writing summary report ==="
    ensure_dir "$WORK_DIR" >/dev/null

    local fq1_size fq2_size
    fq1_size=$(du -h "$FORWARD_READ" 2>/dev/null | awk '{print $1}')
    fq2_size=$(du -h "$REVERSE_READ" 2>/dev/null | awk '{print $1}')

    local denovo_dir="$WORK_DIR/unicycler_out/assembly"
    local refguided_dir="$WORK_DIR/reference_based_assembly"

    local denovo_quast="$WORK_DIR/quast_out/quast_output/report.tsv"
    local refguided_quast="$WORK_DIR/quast2_out/quast_output/report.tsv"

    local best_ref="NA"
    if [[ -f "$WORK_DIR/plentyofbugs_out/best_reference" ]]; then
        best_ref=$(awk 'NR==1 {print $1}' "$WORK_DIR/plentyofbugs_out/best_reference")
    fi

    local busco_short="NA"
    if compgen -G "$WORK_DIR/busco_out/busco_run/short_summary*.txt" > /dev/null; then
        busco_short=$(grep -hE "C:.*S:.*D:.*F:.*M:" \
            "$WORK_DIR/busco_out/busco_run/short_summary"*.txt 2>/dev/null | head -1 | sed 's/^[[:space:]]*//')
    fi

    # Which file was used for ref-guided QC?
    local refguided_used="$ALIGNGRAPH_USED_FILE"
    [[ -n "$refguided_used" ]] || refguided_used="$refguided_dir/sample_remainingcontig.fasta"

    {
        echo "======================================================================"
        echo "                 GENOME ASSEMBLY PIPELINE — SUMMARY REPORT"
        echo "======================================================================"
        echo "Run finished   : $(date)"
        echo "Working dir    : $WORK_DIR"
        echo ""
        echo "----------------------------------------------------------------------"
        echo " INPUT"
        echo "----------------------------------------------------------------------"
        echo "Forward reads  : $FORWARD_READ   ($fq1_size)"
        echo "Reverse reads  : $REVERSE_READ   ($fq2_size)"
        echo ""
        echo "----------------------------------------------------------------------"
        echo " PHASE 1 — QC AND DE NOVO ASSEMBLY"
        echo "----------------------------------------------------------------------"
        echo "FastQC output        : $WORK_DIR/fastqc_out"
        echo "MultiQC report       : $WORK_DIR/multiqc_out/multiqc_report.html"
        echo ""
        echo "Trimmomatic params   : threads=$THREADS phred=$PHRED"
        echo "                       adapter=$ADAPTER_FILE"
        echo "                       LEADING=$LEADING TRAILING=$TRAILING"
        echo "                       SLIDINGWINDOW=$SLIDINGWINDOW MINLEN=$MINLEN"
        echo "Trimmed reads        : $WORK_DIR/trimmomatic_out"
        echo ""
        echo "FLASH merged reads   : $WORK_DIR/flash_out/out.extendedFrags.fastq"
        echo "                       max-overlap = $MAX_OVERLAP"
        echo ""
        echo "De novo assembly     : $denovo_dir/assembly.fasta"
        echo "QUAST (de novo)      : $WORK_DIR/quast_out/quast_output/report.tsv"
        echo "  # contigs          : $(quast_metric "$denovo_quast" '# contigs')"
        echo "  Largest contig     : $(quast_metric "$denovo_quast" 'Largest contig')"
        echo "  Total length       : $(quast_metric "$denovo_quast" 'Total length')"
        echo "  N50                : $(quast_metric "$denovo_quast" 'N50')"
        echo "  L50                : $(quast_metric "$denovo_quast" 'L50')"
        echo "  GC (%)             : $(quast_metric "$denovo_quast" 'GC (%)')"
        echo ""

        if [[ "$SKIP_PHASE2" -eq 1 ]]; then
            echo "----------------------------------------------------------------------"
            echo " PHASE 2 — SKIPPED (-S)"
            echo "----------------------------------------------------------------------"
        else
            echo "----------------------------------------------------------------------"
            echo " PHASE 2 — REFERENCE-BASED EVALUATION"
            echo "----------------------------------------------------------------------"
            echo "Candidate refs dir   : $REFERENCE_DIR  (${#CANDIDATE_REFS[@]} FASTA files)"
            echo "Best reference       : $best_ref"
            echo "                       (chosen by plentyofbugs via Mash)"
            echo "plentyofbugs output  : $WORK_DIR/plentyofbugs_out"
            echo "Assembler used       : $ASSEMBLER"
            echo ""
            echo "Bowtie2 alignment"
            echo "  SAM                : $WORK_DIR/bowtie2_out/out.sam"
            echo "  Report             : $WORK_DIR/bowtie2_out/alignment_report.txt"
            echo "  Total records      : $ALIGN_TOTAL"
            echo "  Aligned records    : $ALIGN_ALIGNED"
            echo "  Alignment rate     : ${ALIGN_RATE}%"
            echo "  Threshold          : ${ALIGN_THRESHOLD}%"
            echo ""
            echo "AlignGraph output    : $refguided_dir"
            echo "  Status             : $ALIGNGRAPH_STATUS"
            echo "  Extended contigs   : $refguided_dir/sample_extendedcontig.fasta"
            echo "  Remaining contigs  : $refguided_dir/sample_remainingcontig.fasta"
            echo "  Used for QC        : $refguided_used"
            echo "  distanceLow        : $DIST_LOW"
            echo "  distanceHigh       : $DIST_HIGH"
            echo ""
            echo "QUAST (ref-guided)   : $WORK_DIR/quast2_out/quast_output/report.tsv"
            echo "  # contigs          : $(quast_metric "$refguided_quast" '# contigs')"
            echo "  Largest contig     : $(quast_metric "$refguided_quast" 'Largest contig')"
            echo "  Total length       : $(quast_metric "$refguided_quast" 'Total length')"
            echo "  N50                : $(quast_metric "$refguided_quast" 'N50')"
            echo "  L50                : $(quast_metric "$refguided_quast" 'L50')"
            echo "  GC (%)             : $(quast_metric "$refguided_quast" 'GC (%)')"
            echo ""
            echo "BUSCO completeness   : $WORK_DIR/busco_out/busco_run"
            echo "  Lineage            : $BUSCO_LINEAGE"
            echo "  Summary            : $busco_short"
            echo ""

            if [[ "$ALIGNGRAPH_STATUS" == "FALLBACK" ]]; then
                echo "NOTE: AlignGraph produced no extended contigs (0-byte extendedContig)."
                echo "      The pipeline automatically used remainingContig for QUAST/BUSCO."
                echo "      This indicates the de novo assembly was already complete"
                echo "      relative to the reference; no scaffolding was necessary."
                echo ""
            fi

            echo "----------------------------------------------------------------------"
            echo " COMPARISON: de novo vs reference-guided"
            echo "----------------------------------------------------------------------"
            local dn_n50 rf_n50 dn_c rf_c dn_len rf_len
            dn_n50=$(quast_metric "$denovo_quast" 'N50')
            rf_n50=$(quast_metric "$refguided_quast" 'N50')
            dn_c=$(quast_metric "$denovo_quast" '# contigs')
            rf_c=$(quast_metric "$refguided_quast" '# contigs')
            dn_len=$(quast_metric "$denovo_quast" 'Total length')
            rf_len=$(quast_metric "$refguided_quast" 'Total length')
            printf "                         %-14s %s\n" "de novo" "reference-guided"
            printf "  # contigs           : %-14s %s\n" "$dn_c"   "$rf_c"
            printf "  Total length        : %-14s %s\n" "$dn_len" "$rf_len"
            printf "  N50                 : %-14s %s\n" "$dn_n50" "$rf_n50"
        fi

        echo ""
        echo "----------------------------------------------------------------------"
        echo " KEY OUTPUT FILES"
        echo "----------------------------------------------------------------------"
        echo "Final de novo asm    : $denovo_dir/assembly.fasta"
        if [[ "$SKIP_PHASE2" -eq 0 ]]; then
            echo "Final ref-guided asm : $refguided_used"
        fi
        echo "Full log             : $LOG_FILE"
        echo "This report          : $REPORT_FILE"
        echo ""
        echo "======================================================================"
    } > "$REPORT_FILE"

    log "Summary report written -> $REPORT_FILE"
    echo ""
    cat "$REPORT_FILE"
}

###############################################################################
# Main
###############################################################################
main() {
    parse_args "$@"

    local CURRENT_DIR
    CURRENT_DIR=$(pwd)

    if [[ -n "$RESUME_DIR" ]]; then
        WORK_DIR="$RESUME_DIR"
    else
        WORK_DIR="$CURRENT_DIR/assembly_pipeline_$(date +%Y%m%d_%H%M%S)"
    fi
    CHECKPOINT_DIR="$WORK_DIR/.checkpoints"
    LOG_FILE="$WORK_DIR/pipeline.log"
    REPORT_FILE="$WORK_DIR/summary_report.txt"
    CONFIG_FILE="$WORK_DIR/.pipeline_config"

    if [[ -n "$RESUME_DIR" ]]; then
        ensure_dir "$WORK_DIR" >/dev/null
        local cli_fwd="$FORWARD_READ"
        local cli_rev="$REVERSE_READ"
        local cli_ref="$REFERENCE_DIR"
        local cli_skip="$SKIP_PHASE2"
        load_config
        [[ -n "$cli_fwd" ]] && FORWARD_READ="$cli_fwd"
        [[ -n "$cli_rev" ]] && REVERSE_READ="$cli_rev"
        [[ -n "$cli_ref" ]] && REFERENCE_DIR="$cli_ref"
        [[ "$cli_skip" -eq 1 ]] && SKIP_PHASE2=1
    fi

    validate_inputs

    ensure_dir "$WORK_DIR" >/dev/null
    ensure_dir "$CHECKPOINT_DIR" >/dev/null
    : >> "$LOG_FILE"

    log "Pipeline started: $(date)"
    log "Working dir   : $WORK_DIR"
    if [[ -n "$RESUME_DIR" ]]; then
        log "RESUME MODE — skipping completed steps"
    fi

    save_config

    # ---- Phase 1 ----
    run_step "01_fastqc"       run_fastqc
    run_step "02_multiqc"      run_multiqc
    run_step "03_trimmomatic"  run_trimmomatic
    run_step "04_flash"        run_flash
    run_step "05_unicycler"    run_unicycler
    run_step "06_quast_denovo" run_quast "de-novo" \
        "$WORK_DIR/unicycler_out/assembly/assembly.fasta" \
        "$WORK_DIR/quast_out"

    log "Phase 1 complete."

    # ---- Phase 2 ----
    if [[ "$SKIP_PHASE2" -eq 1 ]]; then
        log "Phase 2 skipped (-S)."
        write_summary_report
        log "Pipeline finished: $(date)"
        exit 0
    fi

    run_step "07_plentyofbugs" run_plentyofbugs
    run_step "08_bowtie2"      run_bowtie2
    run_step "09_aligngraph"   run_aligngraph

    # After AlignGraph: if we skipped it (checkpoint present), recompute status
    if step_done "09_aligngraph" && [[ "$ALIGNGRAPH_STATUS" == "not run" ]]; then
        local ext="$WORK_DIR/reference_based_assembly/sample_extendedcontig.fasta"
        local rem="$WORK_DIR/reference_based_assembly/sample_remainingcontig.fasta"
        if [[ -s "$ext" ]]; then
            ALIGNGRAPH_STATUS="OK"
            ALIGNGRAPH_USED_FILE="$ext"
        elif [[ -s "$rem" ]]; then
            ALIGNGRAPH_STATUS="FALLBACK"
            ALIGNGRAPH_USED_FILE="$rem"
        fi
    fi

    # Choose ref-guided input for QUAST
    local refguided_input="$WORK_DIR/reference_based_assembly/sample_extendedcontig.fasta"
    if [[ ! -s "$refguided_input" ]]; then
        refguided_input="$WORK_DIR/reference_based_assembly/sample_remainingcontig.fasta"
        [[ -s "$refguided_input" ]] || refguided_input="$WORK_DIR/unicycler_out/assembly/assembly.fasta"
    fi

    run_step "10_quast_refguided" run_quast "reference-guided" \
        "$refguided_input" \
        "$WORK_DIR/quast2_out"

    run_step "11_busco" run_busco

    log "Phase 2 complete."

    write_summary_report

    log "Pipeline finished: $(date)"
    log "Results: $WORK_DIR"
}

main "$@"
