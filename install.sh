#!/bin/bash
set -euo pipefail

# Copy environment template to root
cp example/environment.sh environment.sh

# Source it to get LCG (temporarily allow unset variables in LCG scripts)
set +u
source environment.sh
set -u

# Create venv with access to LCG system packages
python3 -m venv --system-site-packages .venv
source .venv/bin/activate

# Install package in editable mode
python3 -m pip install --upgrade --upgrade-strategy eager --editable .

# Recompile correctionlib against our LCG
python3 -m pip install --force-reinstall --no-cache-dir --no-binary correctionlib --no-deps correctionlib

# Patch the placeholder path in environment.sh
VENV_PATH="$(realpath .venv)"
sed -i "s|# source /INSERT/ABSOLUTE/PATH/TO/YOUR/ENVIRONMENT/bin/activate|source ${VENV_PATH}/bin/activate|" environment.sh

echo "Done. Run 'source environment.sh' to activate."