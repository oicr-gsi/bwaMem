#!/bin/bash
set -o nounset
set -o pipefail

cd $1

module load samtools/1.9 2>/dev/null

ls | sort

find -name *.bam -exec samtools view -H {} \; | grep '^@RG' | sort

find -name *.bam -exec samtools flagstat {} \; | sort

