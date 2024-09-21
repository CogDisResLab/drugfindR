#!/usr/bin/env bash

set -euo pipefail

# This script creates a tag for the current commit
# It reads the current version string from the DESCRIPTION file

DESCRIPTION_FILE="DESCRIPTION"

# Extract the version from the DESCRIPTION file
VERSION=$(grep '^Version:' "$DESCRIPTION_FILE" | awk '{print $2}')
TAG="v$VERSION"

# Check if the tag already exists
if git rev-parse "$TAG" >/dev/null 2>&1; then
  echo "Tag $TAG already exists. Exiting."
  exit 0
fi

# Create a signed tag if it doesn't already exist
echo "Creating tag for version $VERSION"
git tag -s -a "$TAG" -m "build: New Version $VERSION"

echo "Tag $TAG created"
