#!/bin/bash
# validate_ispcr.sh -- Compare amplirust results against UCSC isPCR
#
# Validates that amplirust produces biologically correct results by comparing
# against the gold-standard UCSC isPCR tool using E. coli K-12 MG1655 with
# 16S universal primers (27F/1492R). E. coli K-12 has 7 copies of the 16S
# rRNA gene, providing a rich multi-hit validation case.
#
# Usage: bash scripts/validate_ispcr.sh
#
# Prerequisites:
#   - isPcr binary in PATH (see install instructions below)
#   - Rust toolchain (cargo) for building amplirust
#   - Internet connection for genome download (first run only)
#
# This is a LOCAL validation script, NOT a CI job. Downloading genomes in
# CI is fragile and slow.

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

CACHE_DIR="/tmp/amplirust_validation"
GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
GENOME="${CACHE_DIR}/ecoli_k12_mg1655.fasta"
ISPCR_QUERY="${CACHE_DIR}/primers_ispcr.txt"
ISPCR_OUTPUT="${CACHE_DIR}/ispcr_output.bed"
AMPLIRUST_OUTPUT="${CACHE_DIR}/amplirust_output.tsv"
AMPLIRUST_FASTA="${CACHE_DIR}/amplirust_output.fasta"

# Resolved 16S primers (no IUPAC ambiguity codes)
# 27F: M resolved to A; 1492R: Y resolved to T
FWD_PRIMER="AGAGTTTGATCATGGCTCAG"
REV_PRIMER="TACGGTTACCTTGTTACGACTT"
PRIMER_NAME="27F_1492R"

# Tolerances for comparison
POSITION_TOLERANCE=5   # +/- 5bp for start coordinate comparison
LENGTH_TOLERANCE=20    # +/- 20bp for product length comparison
EXPECTED_PRODUCTS=7    # E. coli K-12 has 7 copies of 16S rRNA
EXPECTED_LENGTH=1465   # Approximate 16S amplicon length (27F to 1492R)

# ============================================================================
# Functions
# ============================================================================

info() {
    echo "[INFO] $*"
}

warn() {
    echo "[WARN] $*" >&2
}

error() {
    echo "[ERROR] $*" >&2
    exit 1
}

check_ispcr() {
    if command -v isPcr &>/dev/null; then
        info "isPcr found: $(command -v isPcr)"
        return 0
    fi

    echo ""
    echo "============================================================"
    echo "  isPcr binary not found in PATH"
    echo "============================================================"
    echo ""
    echo "Install options:"
    echo ""
    echo "  Option A: Direct download (Linux x86_64)"
    echo "    wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/isPcr"
    echo "    chmod +x isPcr"
    echo "    sudo mv isPcr /usr/local/bin/"
    echo ""
    echo "  Option B: Bioconda"
    echo "    conda install bioconda::ispcr"
    echo ""
    echo "  Option C: macOS"
    echo "    wget https://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/blat/isPcr"
    echo "    chmod +x isPcr"
    echo "    sudo mv isPcr /usr/local/bin/"
    echo ""
    return 1
}

download_genome() {
    mkdir -p "${CACHE_DIR}"

    if [ -f "${GENOME}" ]; then
        info "Genome already cached: ${GENOME}"
        return 0
    fi

    info "Downloading E. coli K-12 MG1655 genome from NCBI..."
    local gz_file="${GENOME}.gz"

    if ! curl -fSL -o "${gz_file}" "${GENOME_URL}"; then
        error "Failed to download genome from ${GENOME_URL}"
    fi

    info "Decompressing genome..."
    gunzip -f "${gz_file}"
    info "Genome cached: ${GENOME} ($(wc -c < "${GENOME}") bytes)"
}

run_ispcr() {
    info "Running isPCR..."

    # Create isPCR query file (tab-separated: name, forward, reverse)
    printf '%s\t%s\t%s\n' "${PRIMER_NAME}" "${FWD_PRIMER}" "${REV_PRIMER}" > "${ISPCR_QUERY}"

    # Run isPCR with BED output (0-based half-open coordinates, matching amplirust)
    isPcr "${GENOME}" "${ISPCR_QUERY}" "${ISPCR_OUTPUT}" -out=bed -maxSize=5000

    local count
    count=$(wc -l < "${ISPCR_OUTPUT}" | tr -d ' ')
    info "isPCR found ${count} products"
}

run_amplirust() {
    info "Running amplirust..."

    local script_dir
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    local project_dir
    project_dir="$(dirname "${script_dir}")"

    local primers_csv="${project_dir}/tests/data/primers_16s_resolved.csv"

    if [ ! -f "${primers_csv}" ]; then
        error "Resolved primers file not found: ${primers_csv}"
    fi

    # Run amplirust with TSV output for coordinate comparison
    cargo run --release --manifest-path "${project_dir}/Cargo.toml" -- \
        -i "${GENOME}" \
        -p "${primers_csv}" \
        -o "${AMPLIRUST_FASTA}" \
        --tsv "${AMPLIRUST_OUTPUT}" \
        --max-len 5000 \
        -k 0 \
        --min-identity 0.9 \
        -q

    local count
    count=$(tail -n +2 "${AMPLIRUST_OUTPUT}" | wc -l | tr -d ' ')
    info "amplirust found ${count} products"
}

compare_results() {
    info "Comparing results..."
    echo ""
    echo "============================================================"
    echo "  Validation Results: amplirust vs UCSC isPCR"
    echo "============================================================"
    echo ""

    local pass=true

    # --- Product count comparison ---
    local ispcr_count amplirust_count
    ispcr_count=$(wc -l < "${ISPCR_OUTPUT}" | tr -d ' ')
    amplirust_count=$(tail -n +2 "${AMPLIRUST_OUTPUT}" | wc -l | tr -d ' ')

    echo "Product count:"
    echo "  isPCR:     ${ispcr_count}"
    echo "  amplirust: ${amplirust_count}"
    echo "  expected:  ${EXPECTED_PRODUCTS}"

    if [ "${ispcr_count}" -ne "${amplirust_count}" ]; then
        warn "Product counts differ!"
        pass=false
    else
        echo "  -> MATCH"
    fi
    echo ""

    # --- Parse isPCR BED output (0-based half-open) ---
    # BED format: chrom, chromStart, chromEnd, name, score, strand
    local ispcr_starts=()
    local ispcr_lengths=()
    while IFS=$'\t' read -r _chrom start end _name _score _strand; do
        ispcr_starts+=("${start}")
        local len=$((end - start))
        ispcr_lengths+=("${len}")
    done < "${ISPCR_OUTPUT}"

    # --- Parse amplirust TSV output ---
    # TSV has a header row; columns include start, end, length
    local amplirust_starts=()
    local amplirust_lengths=()
    local header_skipped=false

    while IFS=$'\t' read -r line_parts; do
        if [ "${header_skipped}" = false ]; then
            header_skipped=true
            continue
        fi
        # Parse the TSV line to extract start, end, and length fields
        local start end length
        start=$(echo "${line_parts}" | cut -f4)
        end=$(echo "${line_parts}" | cut -f5)
        length=$(echo "${line_parts}" | cut -f6)
        amplirust_starts+=("${start}")
        amplirust_lengths+=("${length}")
    done < "${AMPLIRUST_OUTPUT}"

    # --- Position comparison ---
    echo "Position comparison (tolerance: +/- ${POSITION_TOLERANCE}bp):"

    # Sort both sets of starts for comparison
    local ispcr_sorted amplirust_sorted
    ispcr_sorted=$(printf '%s\n' "${ispcr_starts[@]}" | sort -n)
    amplirust_sorted=$(printf '%s\n' "${amplirust_starts[@]}" | sort -n)

    local i=0
    local position_matches=0
    while IFS= read -r ispcr_start; do
        local amplirust_start
        amplirust_start=$(echo "${amplirust_sorted}" | sed -n "$((i + 1))p")

        if [ -z "${amplirust_start}" ]; then
            echo "  isPCR pos ${ispcr_start}: no matching amplirust product"
            pass=false
        else
            local diff=$((ispcr_start - amplirust_start))
            if [ "${diff}" -lt 0 ]; then
                diff=$((-diff))
            fi

            if [ "${diff}" -le "${POSITION_TOLERANCE}" ]; then
                echo "  isPCR=${ispcr_start}  amplirust=${amplirust_start}  diff=${diff}bp  -> MATCH"
                position_matches=$((position_matches + 1))
            else
                echo "  isPCR=${ispcr_start}  amplirust=${amplirust_start}  diff=${diff}bp  -> MISMATCH"
                pass=false
            fi
        fi

        i=$((i + 1))
    done <<< "${ispcr_sorted}"
    echo ""

    # --- Length comparison ---
    echo "Product length comparison (tolerance: +/- ${LENGTH_TOLERANCE}bp, expected ~${EXPECTED_LENGTH}bp):"

    local ispcr_len_sorted amplirust_len_sorted
    ispcr_len_sorted=$(printf '%s\n' "${ispcr_lengths[@]}" | sort -n)
    amplirust_len_sorted=$(printf '%s\n' "${amplirust_lengths[@]}" | sort -n)

    i=0
    local length_matches=0
    while IFS= read -r ispcr_len; do
        local amplirust_len
        amplirust_len=$(echo "${amplirust_len_sorted}" | sed -n "$((i + 1))p")

        if [ -z "${amplirust_len}" ]; then
            echo "  isPCR len ${ispcr_len}: no matching amplirust product"
            pass=false
        else
            local diff=$((ispcr_len - amplirust_len))
            if [ "${diff}" -lt 0 ]; then
                diff=$((-diff))
            fi

            if [ "${diff}" -le "${LENGTH_TOLERANCE}" ]; then
                echo "  isPCR=${ispcr_len}bp  amplirust=${amplirust_len}bp  diff=${diff}bp  -> MATCH"
                length_matches=$((length_matches + 1))
            else
                echo "  isPCR=${ispcr_len}bp  amplirust=${amplirust_len}bp  diff=${diff}bp  -> MISMATCH"
                pass=false
            fi
        fi

        i=$((i + 1))
    done <<< "${ispcr_len_sorted}"
    echo ""

    # --- Summary ---
    echo "============================================================"
    if [ "${pass}" = true ]; then
        echo "  RESULT: PASS"
        echo ""
        echo "  Both tools found ${ispcr_count} products with matching"
        echo "  positions (+/- ${POSITION_TOLERANCE}bp) and lengths (+/- ${LENGTH_TOLERANCE}bp)."
    else
        echo "  RESULT: FAIL"
        echo ""
        echo "  Differences found between isPCR and amplirust output."
        echo "  Review details above."
    fi
    echo "============================================================"

    if [ "${pass}" = false ]; then
        return 1
    fi
}

# ============================================================================
# Main
# ============================================================================

main() {
    echo ""
    echo "Amplirust vs UCSC isPCR Validation"
    echo "==================================="
    echo "Organism: E. coli K-12 MG1655"
    echo "Primers:  27F (${FWD_PRIMER}) / 1492R (${REV_PRIMER})"
    echo "Expected: ${EXPECTED_PRODUCTS} copies of 16S rRNA (~${EXPECTED_LENGTH}bp each)"
    echo ""

    if ! check_ispcr; then
        error "isPcr binary required for validation. See install instructions above."
    fi

    download_genome
    run_ispcr
    run_amplirust
    compare_results
}

main "$@"
