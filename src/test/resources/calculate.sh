#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

ls | sort

/oicr/local/analysis/sw//samtools/samtools-0.1.19/bin/samtools view -H *.bam | grep '^@RG' | sort

/oicr/local/analysis/sw//samtools/samtools-0.1.19/bin/samtools flagstat *.bam | sort
