#!/bin/bash
set -euo pipefail

# Check we're in the repo root by looking for setup.cfg with the expected content
if [ ! -f setup.cfg ] || ! grep -q "name = pepper" setup.cfg 2>/dev/null; then
    echo "error: install.sh must be run from the pepper repo root" >&2
    exit 1
fi


# Copy environment template to root only if it doesn't exist, so that users can customize it without affecting the repo
# Idempotent re-runs become safe and venv can be updated by re-running the script.
if [ ! -f environment.sh ]; then
    cp example/environment.sh environment.sh
fi

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