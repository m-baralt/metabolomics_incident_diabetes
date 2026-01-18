#!/bin/bash

set -e

OWNER="m-baralt"
REPO="metabolomics_incident_diabetes"
TAG="v1.0.0"
DEST="data"

mkdir -p "$DEST"
mkdir -p "figures"
mkdir -p "processed_files"
mkdir -p "results"
mkdir -p "results/discovery"

echo "Fetching release asset list for $OWNER/$REPO $TAG..."

URLS=$(wget -qO- "https://api.github.com/repos/$OWNER/$REPO/releases/tags/$TAG" \
       | grep "browser_download_url" \
       | cut -d '"' -f 4)

if [ -z "$URLS" ]; then
    echo "Error: No assets found for release $TAG"
    exit 1
fi

echo "Downloading files into $DEST/"

for url in $URLS; do
    echo "  $(basename "$url")"
    wget -q -P "$DEST" "$url"
done

echo "All files downloaded successfully."

