#!/usr/bin/bash
set -eou pipefail

in_gfa=$1
out_prefix=$2
threads=8

out_dir=$(dirname "$out_prefix")

mkdir -p $out_dir

vg_anchor get_anchors --dictionary test/large_test/sentinel_to_anchors.pkl --graph test/large_test/chr20 --alignment path/to/alignment.gaf --output path/to/output

vg_anchor verify-output --anchors path/to/output/anchors.json --fastq reads/used/for/alignment.fastq