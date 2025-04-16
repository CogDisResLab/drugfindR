#!/usr/bin/env bash

set -euo pipefail

# --- Configuration ---
DESCRIPTION_FILE="DESCRIPTION"

# --- Extract version from DESCRIPTION file ---
# Assumes the format: Version: X.Y.Z
VERSION=$(grep '^Version:' "$DESCRIPTION_FILE" | awk '{print $2}')
TAG="v$VERSION"

# --- Safety check: ensure we're on a valid commit (not during rebase/merge) ---
if ! git rev-parse HEAD >/dev/null 2>&1; then
  echo "⚠️ Not on a valid commit. Skipping tag creation."
  exit 1
fi

# --- Check if the tag already exists ---
if git rev-parse "$TAG" >/dev/null 2>&1; then
  echo "🏷️ Tag $TAG already exists. Skipping."
  exit 0
fi

# --- Create a signed, annotated tag ---
echo "🏷️ Creating signed tag: $TAG"
git tag -s -a "$TAG" -m "build: New Version $VERSION"

echo "✅ Tag $TAG created successfully."
