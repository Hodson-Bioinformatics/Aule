#!/usr/bin/env bash

# This script will use UCSC liftOver to convert seg or bed files between hg19/hg38 genome builds.
# The scripts accepts input files with any number of columns.
# The argument 'mode' will determine whether it is `SEG` or 'BED' files to convert.
# The argument 'chain' is expecting the path to liftOver chain file.
# The argument 'header' accepts 'YES' or 'NO' to indicate whether header is present in the input file.
# The argument <minmatch> expects a float in the range 0..1 to pass to the liftover minMatch parameter.

# Usage:
# bash liftover.sh <mode> <input file> <output file> <chain> <header> <minmatch>
# Eexample:
# bash liftover.sh BED input.bed output.bed path/to/chain YES 0.9775
# 

set -e

# Read variables to store the arguments from command line
MODE="$1"
INPUT_FILE="$2"
OUTPUT_FILE="$3"
CHAIN="$4"
HEADER="$5"
MINMATCH="$6"


# First, check that proper mode is specified, rearrange columns for seg file, and collapse extra columns together
if [[ "$MODE" == *"SEG"* ]]; then
    echo "Running in the $MODE mode ..."
    cat $INPUT_FILE \
    | perl -pale 'BEGIN { $"="\t"; } $_ = "@F[1,2,3,0,4..$#F]"' \
    | perl -ne '@a=split("\t");$ncol=3;while($ncol>0){$first = shift @a;print "$first\t";$ncol--;}print join "|", @a;' \
    > $OUTPUT_FILE.collapsed
elif [[ "$MODE" == *"BED"* ]]; then
    echo "Running in the $MODE mode ..."
    cat $INPUT_FILE \
    | perl -ne '@a=split("\t");$ncol=3;while($ncol>0){$first = shift @a;print "$first\t";$ncol--;}print join "|", @a;' \
    > $OUTPUT_FILE.collapsed
else
    echo "You specified mode $MODE, which is not supported. Please provide SEG or BED file."
    exit 1 # terminate and indicate error
fi

# Second, if the input file has header, save it temporarily in a separate file
# Prepend with chr if chromosomes are not prepended, othervise liftOver silently outputs empty file
if [[ "$HEADER" == *"YES"* ]]; then
    echo "Header is specified as $HEADER, handling it separately ..."
    head -1 $OUTPUT_FILE.collapsed > $OUTPUT_FILE.header
    tail -n+2 $OUTPUT_FILE.collapsed \
    | perl -lane 'if ( /^chr/ ) { print } else { s/^/chr/;print }' \
    > $OUTPUT_FILE.bed
elif [[ "$HEADER" == *"NO"* ]]; then
    echo "Header is specified as $HEADER, the first entry of $MODE file will not be treated as header. Processing ..."
    cat $OUTPUT_FILE.collapsed \
    | perl -lane 'if ( /^chr/ ) { print } else { s/^/chr/;print }' \
    > $OUTPUT_FILE.bed
else
    echo "You specified header $HEADER, which is not recognized. Please specify YES or NO."
    rm $OUTPUT_FILE.collapsed
    exit 1 # terminate and indicate error
fi

# Now, run the liftOver
echo "Running UCSC liftOver with minimum match of converted regions set to $MINMATCH ..."
UNMAPPED="${OUTPUT_FILE%.*}.unmapped.bed"
liftOver -minMatch=$MINMATCH $OUTPUT_FILE.bed $CHAIN $OUTPUT_FILE.lifted-temp.bed $UNMAPPED

# Next, if the input file had header, merge it back to the lifted file
if [[ "$HEADER" == *"YES"* ]]; then
    cat $OUTPUT_FILE.header > $OUTPUT_FILE.merged
    cat $OUTPUT_FILE.lifted-temp.bed \
    | sort -k1,1 -k2,2n -V \
    >> $OUTPUT_FILE.merged
    rm $OUTPUT_FILE.header # this is only specific if input has header, so clean up this temp file here
else
    cat $OUTPUT_FILE.lifted-temp.bed \
    | sort -k1,1 -k2,2n -V \
    > $OUTPUT_FILE.merged
fi

# Now, split back all concatenated columns into the separate ones and rearrange back if it is SEG file
if [[ "$MODE" == *"SEG"* ]]; then
    cat $OUTPUT_FILE.merged \
    | perl -ne 's/\|/\t/g;print;' \
    | perl -pale 'BEGIN { $"="\t"; } $_ = "@F[3,0..2,4..$#F]"' \
    > $OUTPUT_FILE
else
    cat $OUTPUT_FILE.merged \
    | perl -ne 's/\|/\t/g;print;' \
    > $OUTPUT_FILE
fi

# Cleanup
echo "Completed successfully, cleaning up temp files ..."
rm $OUTPUT_FILE.collapsed
rm $OUTPUT_FILE.bed
rm $OUTPUT_FILE.lifted-temp.bed
rm $OUTPUT_FILE.merged

echo "Done!"
