#!/usr/bin/env bash
set -euo pipefail

# Ensure required tools are available
command -v gh >/dev/null 2>&1 || {
  echo "::error::GitHub CLI (gh) is not installed"
  exit 1
}

# Determine the tag associated with the current commit
TAG_VERSION=$(git tag --points-at HEAD)

if [[ -z "$TAG_VERSION" ]]; then
  echo "::warning::No tag found on HEAD. Falling back to latest reachable tag."
  TAG_VERSION=$(git describe --tags --abbrev=0)
fi

if [[ -z "$TAG_VERSION" ]]; then
  echo "::error::No tag could be determined. Cannot create release."
  exit 1
fi

echo "Creating or overwriting release for tag: $TAG_VERSION"

# Attempt to create the release; overwrite if it exists
gh release delete "$TAG_VERSION" --yes || echo "::notice::No existing release to delete"
gh release create "$TAG_VERSION" \
  --title "Release $TAG_VERSION" \
  --notes "Release $TAG_VERSION" \
  --target "$(git rev-parse HEAD)"

# Output tag to GitHub Actions
echo "tag_version=$TAG_VERSION" >>"$GITHUB_OUTPUT"
