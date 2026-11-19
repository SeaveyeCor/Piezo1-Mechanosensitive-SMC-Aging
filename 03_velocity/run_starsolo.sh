#!/bin/bash
# =============================================================================
# RNA Velocity Pipeline — STARsolo Alignment (25GB RAM)
# =============================================================================
# Aligns 10x Chromium FASTQs with STARsolo to produce BAM files
# suitable for RNA velocity (velocyto/scVelo).
#
# BEFORE RUNNING:
#   1. Run: bash run_starsolo.sh install   (one-time setup)
#   2. Run: bash run_starsolo.sh index     (build genome index, ~45min)
#   3. Run: bash run_starsolo.sh align     (align samples, ~4-8hr each)
#   4. Run: bash run_starsolo.sh velocyto  (generate loom files)
#
# Memory strategy (25GB RAM):
#   - Full (non-sparse) suffix array for fastest alignment
#   - 8 threads
#   - 20GB BAM sort RAM, 24GB index generation
# =============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
WORK_DIR="data/velocity"
FASTQ_DIR="${WORK_DIR}"
GENOME_DIR="${WORK_DIR}/genome"
INDEX_DIR="${WORK_DIR}/star_index"
OUTPUT_DIR="${WORK_DIR}/starsolo_output"
VELOCYTO_DIR="${WORK_DIR}/velocyto_output"

# Resource limits (tuned for 32GB system)
THREADS=8
SORT_RAM=20000000000        # 20GB for BAM sorting
INDEX_RAM=24000000000       # 24GB for genome index generation
STAR_GENOME_LOAD="NoSharedMemory"  # safest for WSL

# Mouse genome (GRCm39 / mm39 — matches 10x reference)
GENOME_URL="https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.gtf.gz"

# 10x barcode whitelist (v3 chemistry — 3M barcodes)
WHITELIST_URL="https://github.com/cellranger/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"

# ── Helper functions ──────────────────────────────────────────────────────────
log() { echo "[$(date '+%H:%M:%S')] $*"; }
check_ram() {
    local avail_gb=$(free -g | awk '/Mem:/{print $7}')
    log "Available RAM: ${avail_gb}GB"
    if [ "$avail_gb" -lt 16 ]; then
        echo ""
        echo "WARNING: Only ${avail_gb}GB available. STAR index generation needs ~24GB."
        echo "Close other applications or increase WSL memory allocation."
        echo ""
        read -p "Continue anyway? (y/n) " -n 1 -r
        echo
        [[ $REPLY =~ ^[Yy]$ ]] || exit 1
    fi
}

# ── Discover samples from FASTQ filenames ─────────────────────────────────────
discover_samples() {
    # Pattern: {SAMPLE}_GEX_R_S{n}_L{lane}_{read}_001.fastq.gz
    # Extract unique sample IDs (e.g., 2181, 2404)
    ls "${FASTQ_DIR}"/*_GEX_R_*_R1_001.fastq.gz 2>/dev/null \
        | sed 's|.*/||; s|_GEX_R_.*||' \
        | sort -u
}

# ── STEP 1: Install dependencies ─────────────────────────────────────────────
do_install() {
    log "Installing STAR, samtools, velocyto dependencies..."

    sudo apt-get update -qq
    sudo apt-get install -y -qq \
        build-essential \
        samtools \
        python3-pip \
        python3-venv \
        wget \
        pigz

    # Install STAR 2.7.11b
    if ! command -v STAR &>/dev/null; then
        log "Installing STAR..."
        cd /tmp
        wget -q https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip
        unzip -q -o STAR_2.7.11b.zip
        sudo cp STAR_2.7.11b/Linux_x86_64_static/STAR /usr/local/bin/
        rm -rf STAR_2.7.11b STAR_2.7.11b.zip
        log "STAR installed: $(STAR --version 2>&1 | head -1)"
    else
        log "STAR already installed: $(STAR --version 2>&1 | head -1)"
    fi

    # Python venv for velocyto
    if [ ! -d "${WORK_DIR}/venv" ]; then
        log "Creating Python venv for velocyto..."
        python3 -m venv "${WORK_DIR}/venv"
        source "${WORK_DIR}/venv/bin/activate"
        pip install -q --upgrade pip
        pip install -q velocyto numpy scipy scvelo anndata
        deactivate
    else
        log "Python venv already exists"
    fi

    log "Installation complete."
}

# ── STEP 2: Download genome + build STAR index ───────────────────────────────
do_index() {
    check_ram
    mkdir -p "${GENOME_DIR}" "${INDEX_DIR}"

    # Download genome FASTA
    if [ ! -f "${GENOME_DIR}/genome.fa" ]; then
        log "Downloading mouse genome (GRCm39)..."
        wget -q --show-progress -O "${GENOME_DIR}/genome.fa.gz" "${GENOME_URL}"
        log "Decompressing genome..."
        pigz -d "${GENOME_DIR}/genome.fa.gz"
    else
        log "Genome FASTA already exists"
    fi

    # Download GTF
    if [ ! -f "${GENOME_DIR}/genes.gtf" ]; then
        log "Downloading gene annotations..."
        wget -q --show-progress -O "${GENOME_DIR}/genes.gtf.gz" "${GTF_URL}"
        pigz -d "${GENOME_DIR}/genes.gtf.gz"
    else
        log "GTF already exists"
    fi

    # Download 10x barcode whitelist
    if [ ! -f "${GENOME_DIR}/3M-february-2018.txt" ]; then
        log "Downloading 10x v3 barcode whitelist..."
        wget -q -O "${GENOME_DIR}/3M-february-2018.txt.gz" "${WHITELIST_URL}" \
            || {
                # Fallback: extract from STAR's built-in
                log "Downloading whitelist from alternate source..."
                wget -q -O "${GENOME_DIR}/3M-february-2018.txt.gz" \
                    "https://raw.githubusercontent.com/ENCODE-DCC/atac-seq-pipeline/master/genome/10x/3M-february-2018.txt.gz"
            }
        pigz -d "${GENOME_DIR}/3M-february-2018.txt.gz"
    else
        log "Barcode whitelist already exists"
    fi

    # Build STAR index
    if [ ! -f "${INDEX_DIR}/SA" ]; then
        log "Building STAR genome index (this takes ~45min)..."
        STAR \
            --runMode genomeGenerate \
            --runThreadN ${THREADS} \
            --genomeDir "${INDEX_DIR}" \
            --genomeFastaFiles "${GENOME_DIR}/genome.fa" \
            --sjdbGTFfile "${GENOME_DIR}/genes.gtf" \
            --sjdbOverhang 90 \
            --limitGenomeGenerateRAM ${INDEX_RAM} \
            --outTmpDir "${WORK_DIR}/_STARtmp_index"
        log "Genome index built successfully."
    else
        log "STAR index already exists"
    fi
}

# ── STEP 3: Align samples (one at a time for RAM safety) ─────────────────────
do_align() {
    check_ram
    mkdir -p "${OUTPUT_DIR}"

    local samples=($(discover_samples))
    if [ ${#samples[@]} -eq 0 ]; then
        log "ERROR: No FASTQ files found matching pattern *_GEX_R_*_R1_001.fastq.gz"
        exit 1
    fi

    log "Found ${#samples[@]} samples: ${samples[*]}"

    for sample in "${samples[@]}"; do
        local sample_out="${OUTPUT_DIR}/${sample}"

        if [ -f "${sample_out}/Aligned.sortedByCoord.out.bam" ]; then
            log "Sample ${sample} already aligned — skipping"
            continue
        fi

        log "═══════════════════════════════════════════════════════"
        log "Aligning sample: ${sample}"
        log "═══════════════════════════════════════════════════════"

        mkdir -p "${sample_out}"

        # Collect R1 (barcode+UMI) and R3 (cDNA) across all lanes
        # 10x dual-index: I1=i7 index, R1=barcode+UMI, R2=i5 index, R3=cDNA
        local r1_files=$(ls ${FASTQ_DIR}/${sample}_GEX_R_*_R1_001.fastq.gz 2>/dev/null | sort | tr '\n' ',')
        local r3_files=$(ls ${FASTQ_DIR}/${sample}_GEX_R_*_R3_001.fastq.gz 2>/dev/null | sort | tr '\n' ',')

        # Remove trailing comma
        r1_files="${r1_files%,}"
        r3_files="${r3_files%,}"

        if [ -z "$r1_files" ] || [ -z "$r3_files" ]; then
            log "WARNING: Missing R1 or R3 files for ${sample} — skipping"
            continue
        fi

        log "R1 (barcode+UMI): ${r1_files}"
        log "R3 (cDNA):        ${r3_files}"

        # Clean up any previous failed run
        rm -rf "${WORK_DIR}/_STARtmp_${sample}"

        STAR \
            --runMode alignReads \
            --runThreadN ${THREADS} \
            --genomeDir "${INDEX_DIR}" \
            --genomeLoad ${STAR_GENOME_LOAD} \
            --readFilesIn ${r3_files} ${r1_files} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --outBAMsortingBinsN 50 \
            --limitBAMsortRAM ${SORT_RAM} \
            --outTmpDir "${WORK_DIR}/_STARtmp_${sample}" \
            --outFileNamePrefix "${sample_out}/" \
            \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist "${GENOME_DIR}/3M-february-2018.txt" \
            --soloCBstart 1 --soloCBlen 16 \
            --soloUMIstart 17 --soloUMIlen 12 \
            --soloFeatures Gene Velocyto \
            --soloOutFileNames "solo/" "genes.tsv" "barcodes.tsv" "matrix.mtx" \
            --soloCellFilter EmptyDrops_CR \
            \
            --outSAMunmapped Within \
            --clipAdapterType CellRanger4 \
            --outFilterScoreMin 30 \
            --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloUMIfiltering MultiGeneUMI_CR \
            --soloUMIdedup 1MM_CR

        # Index the BAM
        log "Indexing BAM for ${sample}..."
        samtools index -@ ${THREADS} "${sample_out}/Aligned.sortedByCoord.out.bam"

        log "Sample ${sample} complete."
        log "BAM: ${sample_out}/Aligned.sortedByCoord.out.bam"
        log "Solo: ${sample_out}/solo/"

        # Report stats
        if [ -f "${sample_out}/Log.final.out" ]; then
            echo ""
            grep -E "Uniquely mapped|Number of input reads|Average mapped length" \
                "${sample_out}/Log.final.out" || true
            echo ""
        fi
    done

    log "All samples aligned."
}

# ── STEP 4: Run velocyto to generate loom files (alternative to STARsolo Velocyto) ──
do_velocyto() {
    source "${WORK_DIR}/venv/bin/activate"
    mkdir -p "${VELOCYTO_DIR}"

    local samples=($(discover_samples))

    for sample in "${samples[@]}"; do
        local bam="${OUTPUT_DIR}/${sample}/Aligned.sortedByCoord.out.bam"
        local loom="${VELOCYTO_DIR}/${sample}.loom"

        if [ ! -f "$bam" ]; then
            log "WARNING: No BAM for ${sample} — skipping velocyto"
            continue
        fi

        if [ -f "$loom" ]; then
            log "Loom file exists for ${sample} — skipping"
            continue
        fi

        log "Running velocyto on ${sample}..."
        velocyto run10x \
            -@ ${THREADS} \
            --samtools-memory 4000 \
            "${OUTPUT_DIR}/${sample}" \
            "${GENOME_DIR}/genes.gtf"

        # Move loom to output dir
        mv "${OUTPUT_DIR}/${sample}/velocyto/"*.loom "${loom}" 2>/dev/null || true
        log "Loom: ${loom}"
    done

    deactivate
    log "Velocyto complete. Loom files in: ${VELOCYTO_DIR}/"
}

# ── STEP 5: Summary / status check ───────────────────────────────────────────
do_status() {
    echo ""
    echo "=== VELOCITY PIPELINE STATUS ==="
    echo ""

    echo "FASTQ files:"
    ls -lh "${FASTQ_DIR}"/*_GEX_R_*.fastq.gz 2>/dev/null | awk '{print "  "$NF, $5}' || echo "  None found"
    echo ""

    echo "Samples detected:"
    for s in $(discover_samples); do
        echo "  ${s}"
        local lanes=$(ls ${FASTQ_DIR}/${s}_GEX_R_*_R1_001.fastq.gz 2>/dev/null | wc -l)
        echo "    Lanes: ${lanes}"
        if [ -f "${OUTPUT_DIR}/${s}/Aligned.sortedByCoord.out.bam" ]; then
            local bam_size=$(ls -lh "${OUTPUT_DIR}/${s}/Aligned.sortedByCoord.out.bam" | awk '{print $5}')
            echo "    BAM: ✓ (${bam_size})"
        else
            echo "    BAM: ✗"
        fi
        if [ -f "${VELOCYTO_DIR}/${s}.loom" ]; then
            echo "    Loom: ✓"
        elif [ -d "${OUTPUT_DIR}/${s}/solo/Velocyto" ]; then
            echo "    STARsolo Velocyto: ✓ (spliced/unspliced matrices in solo/Velocyto/)"
        else
            echo "    Velocity: ✗"
        fi
    done

    echo ""
    echo "RAM: $(free -h | awk '/Mem:/{print $7}') available / $(free -h | awk '/Mem:/{print $2}') total"
    echo "Disk: $(df -h . 2>/dev/null | awk 'NR==2{print $4}' || echo 'N/A') free"
    echo ""
}

# ── Main dispatcher ───────────────────────────────────────────────────────────
case "${1:-status}" in
    install)  do_install ;;
    index)    do_index ;;
    align)    do_align ;;
    velocyto) do_velocyto ;;
    status)   do_status ;;
    all)
        do_install
        do_index
        do_align
        log "Alignment complete. STARsolo Velocyto matrices in solo/Velocyto/."
        log "Run 'bash run_starsolo.sh velocyto' if you also want .loom files."
        ;;
    *)
        echo "Usage: bash run_starsolo.sh {install|index|align|velocyto|status|all}"
        exit 1
        ;;
esac
