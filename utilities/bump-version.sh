#!/bin/bash
set -euo pipefail

# Get the commit message of the latest commit
commit_message=$(git log -1 --pretty=%B)

# Check if the commit message contains "Bump version"
if [[ $commit_message == *"Bump version"* ]]; then
    echo "Already bumped version. Skipping..."
else
    # Run bumpversion
    bumpversion patch
fi
