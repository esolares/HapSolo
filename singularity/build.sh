#!/bin/bash
# Build the HapSolo Singularity image.
#
# Run from the singularity/ directory:
#     cd singularity/
#     ./build.sh
#
# Required: minimap2 and miniprot binaries placed in binaries/<tool>/
# Optional: blat (uncomment in hapsolo.def first if you want to include it)

set -e

cd "$(dirname "$0")"

# Verify required binaries are present
missing=0
for bin in minimap2 miniprot; do
    if [ ! -f "binaries/$bin/$bin" ]; then
        echo "Error: binaries/$bin/$bin not found."
        echo "       See binaries/$bin/README.md for download instructions."
        missing=1
    fi
done
if [ "$missing" -ne 0 ]; then
    echo
    echo "Place the required binaries and run this script again."
    exit 1
fi

# Build
OUTPUT="${1:-hapsolo.sif}"
echo "Building $OUTPUT ..."
sudo singularity build "$OUTPUT" hapsolo.def

echo
echo "Build complete: $OUTPUT"
echo
echo "Test it:"
echo "    singularity run $OUTPUT --help"
