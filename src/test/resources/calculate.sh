#/bin/bash
DIR="$(cd "$( dirname "$0")" && pwd)";
echo $DIR;
cd $1
for i in `find . -type f`; do
	md5sum $i;
	if [[ "$i" == *bam ]]
	    then
	    #Workflow_Bundle_BWA_1.2_SeqWare_0.13.6.13
	    SAMTOOLS=$DIR/../Workflow_Bundle_${workflow-directory-name}_${version}_SeqWare_${seqware-version}/Workflow_Bundle_${workflow-directory-name}/${version}/bin/samtools-${samtools-version}/samtools
	    LINES=`${SAMTOOLS} view ${i} | wc -l`;
	    if [ $LINES = '0' ];  then echo there are no lines in the BAM file; exit 1; fi
	fi;
done
