#!/bin/bash
set -o nounset
set -o pipefail

cd $1

ls | sort

find -name *.sam -exec grep '^@RG' {} \; | sort
find -name *.bam -exec /oicr/local/analysis/sw/samtools/samtools-1.2/bin/samtools view -H {} \; | grep '^@RG' | sort
find -name *.cram -exec /oicr/local/analysis/sw/samtools/samtools-1.2/bin/samtools view -H {} \; | grep '^@RG' | sort

find -name *.sam -exec /oicr/local/analysis/sw/samtools/samtools-1.2/bin/samtools flagstat {} \; | sort
find -name *.bam -exec /oicr/local/analysis/sw/samtools/samtools-1.2/bin/samtools flagstat {} \; | sort
find -name *.cram -exec /oicr/local/analysis/sw/samtools/samtools-1.2/bin/samtools flagstat {} \; | sort
