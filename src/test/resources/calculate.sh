#!/bin/bash

set -o pipefail

for i in `find . -type f`; do

        if [[ "$i" == *bam ]]
            then
            LINES=$(samtools view ${i} | wc -l) || exit $?;
            if [ $LINES = '0' ];  then echo there are no lines in the BAM file; exit 1; fi
        fi;

        md5sum $i;

done
