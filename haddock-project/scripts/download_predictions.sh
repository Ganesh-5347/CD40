#!/bin/bash
# Download structure prediction CIF files from GitHub Release and
# place them directly into haddock-project/inputs/.
#
# Run from the repo root:
#   bash haddock-project/scripts/download_predictions.sh
#
# This populates:
#   haddock-project/inputs/
#   ├── af3/    (190 CIFs — 38 candidates × 5 seeds)
#   ├── boltz/  (190 CIFs — 38 candidates × 5 samples)
#   └── af2/    (950 CIFs — 38 candidates × 5 architectures × 5 seeds)

URL="https://github.com/Ganesh-5347/CD40/releases/download/v1.0/structure_predictions.tar.gz"
OUT="structure_predictions.tar.gz"

echo "Downloading structure predictions (~82 MB)..."
curl -L -o "$OUT" "$URL"

echo "Extracting..."
tar xzf "$OUT"

# Move into pipeline input directory
mkdir -p haddock-project/inputs
mv structure_predictions/af3 haddock-project/inputs/
mv structure_predictions/af2 haddock-project/inputs/
mv structure_predictions/boltz haddock-project/inputs/
rm -rf structure_predictions "$OUT"

echo "Done. $(find haddock-project/inputs/ -name '*.cif' | wc -l) CIF files in haddock-project/inputs/."

# Clone HADDOCK3 (required for refinement step)
if [ ! -d "haddock3" ]; then
    echo "Cloning HADDOCK3..."
    git clone https://github.com/haddocking/haddock3.git
    echo "HADDOCK3 cloned. Install with: pip install -e haddock3/"
else
    echo "haddock3/ already exists, skipping."
fi
