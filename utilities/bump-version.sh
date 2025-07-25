#!/bin/bash
set -euo pipefail

# --- Configuration ---

DESCRIPTION_FILE="DESCRIPTION"
CODEMETA_FILE="codemeta.json"

# --- Get base version (e.g., 0.99) from DESCRIPTION file ---
# Assumes DESCRIPTION contains: Version: X.Y.Z
BASE_VERSION=$(grep '^Version:' "$DESCRIPTION_FILE" | awk '{print $2}' | cut -d. -f1,2)

# --- Use total commit count as patch number ---
PATCH_NUMBER=$(git rev-list --all --count HEAD)

# --- Compose new version: MAJOR.MINOR.PATCH ---
NEW_VERSION="$BASE_VERSION.$PATCH_NUMBER"

# --- Get current version from DESCRIPTION ---
CURRENT_VERSION=$(grep '^Version:' "$DESCRIPTION_FILE" | awk '{print $2}' || true)

# --- If already up to date, skip ---
if [[ "$CURRENT_VERSION" == "$NEW_VERSION" ]]; then
    echo "Version is already up-to-date: $NEW_VERSION"
    exit 0
fi

echo "🔧 Updating version from $CURRENT_VERSION to $NEW_VERSION"

# --- Set environment variable to skip CI jobs if needed ---
export SKIP="bump-version,codemeta-json-updated"

# --- Update DESCRIPTION file ---
# macOS `sed` needs backup extension, GNU sed does not
if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' "s/^Version: $CURRENT_VERSION/Version: $NEW_VERSION/" "$DESCRIPTION_FILE"
else
    sed -i "s/^Version: $CURRENT_VERSION/Version: $NEW_VERSION/" "$DESCRIPTION_FILE"
fi

# --- Update codemeta.json using jq ---
# Overwrites the version key
tmpfile=$(mktemp)
jq --arg new_version "$NEW_VERSION" '.version = $new_version' "$CODEMETA_FILE" >"$tmpfile"
mv "$tmpfile" "$CODEMETA_FILE"

# --- Stage updated files for git commit ---
git add "$DESCRIPTION_FILE" "$CODEMETA_FILE"

echo "✅ Version updated and staged: $NEW_VERSION"
