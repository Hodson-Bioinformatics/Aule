#!/bin/bash

# This lcr-script is adopted from Bruno's convert_maf_coords.sh script available at Morin Lab Scripts
# github repo. Some modifications are adopted to allow its implementation as a part of vcf2maf lcr-module.

#   Usage:
#       convert_maf_coords.sh input.maf hg38ToHg19.over.chain.gz output.maf [crossmap|liftover]
#
#   The default mode uses CrossMap.py to convert between genomic coordinate systems, because it
#   tends to drop fewer variants than liftOver. 
#
#   The script assumes that there are no pipe (|) characters in your MAF file. If this is not
#   the case, you can specify your own MAF column separator using the MAFCOLSEP variable.
#
#   By default, the script expects CrossMap.py or liftOver to be in your PATH. If that's not the
#   case, you can specify the program path using the CROSSMAP or LIFTOVER environment variables.
#
# ----------------------------------------------------------------------------------------------- #


set -euf -o pipefail

# Parameters
SCRIPT="$(basename $0)"
INPUT_MAF="$1"
INPUT_BED="$(mktemp /tmp/${SCRIPT}.input.XXXXXX)"
CHAIN_FILE="$2"
OUTPUT_BED="$(mktemp /tmp/${SCRIPT}.output.XXXXXX)"
OUTPUT_MAF="$3"
OUTPUT_UNMAPPED="${OUTPUT_MAF%.*}.unmapped.bed"
MODE="${4:-crossmap}"
MAFCOLSEP="${MAFCOLSEP:-|}"
DEBUG="${DEBUG:-0}"

# For colour
RED='\033[0;31m'
NC='\033[0m'

DATE=$(date --rfc-3339 seconds | cut -c1-19)
echo "[${DATE}] INFO: Checking whether there are \"${MAFCOLSEP}\" characters in the input" \
	"MAF file."
if fgrep -q \'${MAFCOLSEP}\' "${INPUT_MAF}"; then
	echo -e "${RED}[${DATE}] ERROR: ${INPUT_MAF} contains \"${MAFCOLSEP}\" characters, which" \
		"are used to separate MAF columns in the intermediate BED files. Please set the" \
		"MAFCOLSEP environment variable to change this.${NC}"; exit 1;
fi

DATE=$(date --rfc-3339 seconds | cut -c1-19)
echo "[${DATE}] INFO: Creating temporary BED file (${INPUT_BED})" \
	"from input MAF file."
awk 'BEGIN {FS=OFS="\t"} \
		{chrom = $5; start = $6; end = $7} \
		$10 == "SNP" {start = start - 1; end = end} \
		$10 == "DEL" {start = start - 1; end = end} \
		$10 == "INS" {start = start; end = end - 1} \
		$0 !~ /^(#|Hugo_Symbol)/ {row=gensub(/\t/, "'"${MAFCOLSEP}"'", "g", $0); print chrom, start, end, row}' \
		"${INPUT_MAF}" | \
	sed 's/ /_____/g' > "${INPUT_BED}"

if [[ "${MODE}" == "crossmap" ]]; then
	CROSSMAP="${CROSSMAP:-CrossMap.py}"
	command -v "${CROSSMAP}" >/dev/null 2>&1 || { \
		echo -e "${RED}[${DATE}] ERROR: ${CROSSMAP} is not in the PATH or the file doesn't" \
			"exist.${NC}"; exit 1; }
	DATE=$(date --rfc-3339 seconds | cut -c1-19)
	echo "[${DATE}] INFO: Converting genomic coordinates using ${CROSSMAP}."
	"${CROSSMAP}" bed "${CHAIN_FILE}" "${INPUT_BED}" "${OUTPUT_BED}"
	mv -f "${OUTPUT_BED}.unmap" "${OUTPUT_UNMAPPED}"
else
	DATE=$(date --rfc-3339 seconds | cut -c1-19)
	LIFTOVER="${LIFTOVER:-liftOver}"
	command -v "${LIFTOVER}" >/dev/null 2>&1 || { \
		echo -e "${RED}[${DATE}] ERROR: ${LIFTOVER} is not in the PATH or the file doesn't" \
			"exist.${NC}"; exit 1; }
	echo "[${DATE}] INFO: Converting genomic coordinates using ${LIFTOVER}." \
	"${LIFTOVER}" "${INPUT_BED}" "${CHAIN_FILE}" \
		"${OUTPUT_BED}" "${OUTPUT_UNMAPPED}"
fi

NUM_UNMAPPED=$(wc -l < "${OUTPUT_UNMAPPED}")
if [[ ${NUM_UNMAPPED} > 0 ]]; then
	DATE=$(date --rfc-3339 seconds | cut -c1-19)
	echo -e "${RED}[${DATE}] WARNING: ${NUM_UNMAPPED} variants were not properly mapped to the" \
		"target genomic coordinate system.${NC}"
fi

DATE=$(date --rfc-3339 seconds | cut -c1-19)
echo "[${DATE}] INFO: Converting output BED file (${OUTPUT_BED}) to MAF format."
head "${INPUT_MAF}" | grep "^Hugo_Symbol" "${INPUT_MAF}" > "${OUTPUT_MAF}"
awk 'BEGIN {FS=OFS="\t"} \
		{row = gensub(/['"${MAFCOLSEP}"']/, "\t", "g", $4)} \
		{print $1, $2, $3, row}' \
		"${OUTPUT_BED}" \
	# The output needs to be always chr-prefixed to work with vcf2maf 1.3 and reannotate_maf
	| awk -v chain="${CHAIN_FILE}" 'BEGIN {FS=OFS="\t"} \
		{if ($1 ~ /chr/) chrom =$1; else chrom ="chr"$1; start = $2; end = $3} \ 
		$13 == "SNP" {start = start + 1; end = end} \
		$13 == "DEL" {start = start + 1; end = end} \
		$13 == "INS" {start = start; end = end + 1} \
		{$8 = chrom; $9 = start; $10 = end; print $0}' \
	| cut -f4- \
	| awk -v chain="${CHAIN_FILE}" 'BEGIN {FS=OFS="\t"}; { if (chain ~ /hg38ToHg19/) $4 ="GRCh37"; else $4="GRCh38"; print ($0)}' \
	| sed 's/_____/ /g' \
  | perl -F'\t' -ane 'print if $F[4] =~ /^(chr)*[\dX]{1,2}$/' >> "${OUTPUT_MAF}"

if [[ ${DEBUG} == 0 ]]; then
	DATE=$(date --rfc-3339 seconds | cut -c1-19)
	echo "[${DATE}] INFO: Cleaning up temporary BED files."
	rm -f "${INPUT_BED}" "${OUTPUT_BED}"
fi

DATE=$(date --rfc-3339 seconds | cut -c1-19)
echo -e "[${DATE}] INFO: Done! You can find the following output files:
    MAF File:          ${OUTPUT_MAF}
    Unmapped Variants: ${OUTPUT_UNMAPPED}"
