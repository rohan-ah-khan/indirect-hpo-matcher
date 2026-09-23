#!/usr/bin/env bash
set -e

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
destination="${1:-$script_dir/references}"
mkdir -p "$destination"

wget -O "$destination/genes_to_phenotype.txt" \
    https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt
wget -O "$destination/hp.json" \
    https://purl.obolibrary.org/obo/hp.json
