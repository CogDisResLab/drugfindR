#!/bin/bash
set -euo pipefail

# Check if there commitizen binary on PATH
if ! [ -x "$(command -v cz)" ]; then
    echo "Commitizen is not installed. Please install it with 'npm install -g commitizen'."
    exit 1
fi

# Get the most recent tag
tag=$(git describe --abbrev=0 --tags)

# Update the changelog. Assume this is not the first tag
cz changelog
