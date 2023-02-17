# bwaMem

BwaMem Workflow version 2.0

## Overview

## Dependencies

* [bwa 0.7.12](https://github.com/lh3/bwa/archive/0.7.12.tar.gz)
* [samtools 1.9](https://github.com/samtools/samtools/archive/0.1.19.tar.gz)
* [cutadapt 1.8.3](https://cutadapt.readthedocs.io/en/v1.8.3/)
* [slicer 0.3.0](https://github.com/OpenGene/slicer/archive/v0.3.0.tar.gz)
* [barcodex-rs 0.1.2](https://github.com/oicr-gsi/barcodex-rs/archive/v0.1.2.tar.gz)
* [rust 1.2](https://www.rust-lang.org/tools/install)


## Usage

### Cromwell
```
java -jar cromwell.jar run bwaMem.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`fastqR1`|File|Fastq file for read 1
`reference`|String|The genome reference build. For example: hg19, hg38, mm10
`runBwaMem.readGroups`|String|Array of readgroup lines


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqR2`|File?|None|Fastq file for read 2
`outputFileNamePrefix`|String|"output"|Prefix for output file
`numChunk`|Int|1|Number of chunks to split fastq file [1, no splitting]
`doUMIextract`|Boolean|false|If true, UMI will be extracted before alignment [false]
`doTrim`|Boolean|false|If true, adapters will be trimmed before alignment [false]


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`countChunkSize.jobMemory`|Int|16|Memory allocated for this job
`countChunkSize.timeout`|Int|48|Hours before task timeout
`slicerR1.modules`|String|"slicer/0.3.0"|Required environment modules
`slicerR1.jobMemory`|Int|16|Memory allocated for this job
`slicerR1.timeout`|Int|48|Hours before task timeout
`slicerR2.modules`|String|"slicer/0.3.0"|Required environment modules
`slicerR2.jobMemory`|Int|16|Memory allocated for this job
`slicerR2.timeout`|Int|48|Hours before task timeout
`extractUMIs.umiList`|String|"umiList"|Reference file with valid UMIs
`extractUMIs.outputPrefix`|String|"extractUMIs_output"|Specifies the start of the output files
`extractUMIs.pattern1`|String|"(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"|UMI RegEx pattern 1
`extractUMIs.pattern2`|String|"(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"|UMI RegEx pattern 2
`extractUMIs.modules`|String|"barcodex-rs/0.1.2 rust/1.45.1"|Required environment modules
`extractUMIs.jobMemory`|Int|24|Memory allocated for this job
`extractUMIs.timeout`|Int|12|Time in hours before task timeout
`adapterTrimming.modules`|String|"cutadapt/1.8.3"|Required environment modules
`adapterTrimming.doUMItrim`|Boolean|false|If true, do umi trimming
`adapterTrimming.umiLength`|Int|5|The number of bases to trim. If the given length is positive, the bases are removed from the beginning of each read. If it is negative, the bases are removed from the end
`adapterTrimming.trimMinLength`|Int|1|Minimum length of reads to keep
`adapterTrimming.trimMinQuality`|Int|0|Minimum quality of read ends to keep
`adapterTrimming.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|Adapter sequence to trim from read 1
`adapterTrimming.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|Adapter sequence to trim from read 2
`adapterTrimming.addParam`|String?|None|Additional cutadapt parameters
`adapterTrimming.jobMemory`|Int|16|Memory allocated for this job
`adapterTrimming.timeout`|Int|48|Hours before task timeout
`runBwaMem.addParam`|String?|None|Additional BWA parameters
`runBwaMem.threads`|Int|8|Requested CPU threads
`runBwaMem.jobMemory`|Int|32|Memory allocated for this job
`runBwaMem.timeout`|Int|96|Hours before task timeout
`bamMerge.jobMemory`|Int|32|Memory allocated indexing job
`bamMerge.modules`|String|"samtools/1.9"|Required environment modules
`bamMerge.timeout`|Int|72|Hours before task timeout
`indexBam.jobMemory`|Int|12|Memory allocated indexing job
`indexBam.modules`|String|"samtools/1.9"|Modules for running indexing job
`indexBam.timeout`|Int|48|Hours before task timeout
`adapterTrimmingLog.jobMemory`|Int|12|Memory allocated indexing job
`adapterTrimmingLog.timeout`|Int|48|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`bwaMemBam`|File|output merged bam aligned to genome
`bwaMemIndex`|File|output index file for bam aligned to genome
`log`|File?|a summary log file for adapter trimming
`cutAdaptAllLogs`|File?|a file containing all logs for adapter trimming for each fastq chunk


## Commands
 This section lists command(s) run by bwaMem workflow
 
 * Running WORKFLOW
 
 === Description here ===.
 
 <<<
         set -euo pipefail
         totalLines=$(zcat ~{fastqR1} | wc -l)
         python -c "from math import ceil; print int(ceil(($totalLines/4.0)/~{numChunk})*4)"
     >>>
 <<<
         set -euo pipefail
         slicer -i ~{fastqR} -l ~{chunkSize} --gzip 
     >>>
 <<<
             set -euo pipefail
 
             barcodex-rs --umilist ~{umiList} --prefix ~{outputPrefix} --separator "__" inline \
             --pattern1 '~{pattern1}' --r1-in ~{fastq1} \
             ~{if (defined(fastq2)) then "--pattern2 '~{pattern2}' --r2-in ~{fastq2} " else ""}
 
             cat ~{outputPrefix}_UMI_counts.json > umiCounts.txt
 
             tr [,] ',\n' < umiCounts.txt | sed 's/[{}]//' > tmp.txt
             echo "{$(sort -i tmp.txt)}" > new.txt
             tr '\n' ',' < new.txt | sed 's/,$//' > ~{outputPrefix}_UMI_counts.json
         >>>
 <<<
         set -euo pipefail
 
         cutadapt -q ~{trimMinQuality} \
                 -m ~{trimMinLength} \
                 -a ~{adapter1} \
                 -o ~{resultFastqR1} \
                 ~{if (defined(fastqR2)) then "-A ~{adapter2} -p ~{resultFastqR2} " else ""} \
                 ~{if (doUMItrim) then "-u ~{umiLength} -U ~{umiLength} " else ""} \
                 ~{addParam} \
                 ~{fastqR1} \
                 ~{fastqR2} > ~{resultLog}
 
     >>>
 <<<
         set -euo pipefail
         mkdir -p ~{tmpDir}
         bwa mem -M \
             -t ~{threads} ~{addParam}  \
             -R  ~{readGroups} \
             ~{bwaRef} \
             ~{read1s} \
             ~{read2s} \
         | \
         samtools sort -O bam -T ~{tmpDir} -o ~{resultBam} - 
     >>>
 <<<
         set -euo pipefail
         samtools merge \
         -c \
         ~{resultMergedBam} \
         ~{sep=" " bams} 
     >>>
 <<<
         set -euo pipefail
         samtools index ~{inputBam} ~{resultBai}
     >>>
 <<<
         set -euo pipefail
         awk 'BEGINFILE {print "###################################\n"}{print}' ~{sep=" " inputLogs} > ~{allLog}
 
         totalBP=$(cat ~{allLog} | grep "Total basepairs processed:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
 
         bpQualitytrimmed=$(cat ~{allLog} | grep "Quality-trimmed:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g; s/ (.*)//' | awk '{x+=$1}END{print x}')
         percentQualitytrimmed=$(awk -v A="${bpQualitytrimmed}" -v B="${totalBP}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')
 
         bpTotalWritten=$(cat ~{allLog} | grep "Total written (filtered):" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g; s/ (.*)//' | awk '{x+=$1}END{print x}')
         percentBPWritten=$(awk -v A="${bpTotalWritten}" -v B="${totalBP}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')
 
         echo -e "This is a cutadapt summary from ~{numChunk} fastq chunks\n" > ~{log}
 
         if ! ~{singleEnded} ; then
           totalRead=$(cat ~{allLog} | grep "Total read pairs processed:" | cut -d":" -f2 | sed 's/ //g; s/,//g' | awk '{x+=$1}END{print x}')
           adapterR1=$(cat ~{allLog} | grep " Read 1 with adapter:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
           percentAdapterR1=$(awk -v A="${adapterR1}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')
           adapterR2=$(cat ~{allLog} | grep " Read 2 with adapter:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
           percentAdapterR2=$(awk -v A="${adapterR2}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')
 
           shortPairs=$(cat ~{allLog} | grep "Pairs that were too short:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
           percentShortPairs=$(awk -v A="${shortPairs}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')
 
           pairsWritten=$(cat ~{allLog} | grep "Pairs written (passing filters): " | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
           percentpairsWritten=$(awk -v A="${pairsWritten}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')
 
           bpR1=$(cat ~{allLog} | grep -A 2 "Total basepairs processed:" | grep "Read 1:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
           bpR2=$(cat ~{allLog} | grep -A 2 "Total basepairs processed:" | grep "Read 2:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
 
           bpQualitytrimmedR1=$(cat ~{allLog} | grep -A 2 "Quality-trimmed:" | grep "Read 1:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
           bpQualitytrimmedR2=$(cat ~{allLog} | grep -A 2 "Quality-trimmed:" | grep "Read 2:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
 
           bpWrittenR1=$(cat ~{allLog} | grep -A 2 "Total written (filtered):" | grep "Read 1:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
           bpWrittenR2=$(cat ~{allLog} | grep -A 2 "Total written (filtered):" | grep "Read 2:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
 
           echo -e "Total read pairs processed:\t${totalRead}" >> ~{log}
           echo -e "  Read 1 with adapter:\t${adapterR1} (${percentAdapterR1}%)" >> ~{log}
           echo -e "  Read 2 with adapter:\t${adapterR2} (${percentAdapterR2}%)" >> ~{log}
           echo -e "Pairs that were too short:\t${shortPairs} (${percentShortPairs}%)" >> ~{log}
           echo -e "Pairs written (passing filters):\t${pairsWritten} (${percentpairsWritten}%)\n\n" >> ~{log}
           echo -e "Total basepairs processed:\t${totalBP} bp" >> ~{log}
           echo -e "  Read 1:\t${bpR1} bp" >> ~{log}
           echo -e "  Read 2:\t${bpR2} bp" >> ~{log}
           echo -e "Quality-trimmed:\t${bpQualitytrimmed} bp (${percentQualitytrimmed}%)" >> ~{log}
           echo -e "  Read 1:\t${bpQualitytrimmedR1} bp" >> ~{log}
           echo -e "  Read 2:\t${bpQualitytrimmedR2} bp" >> ~{log}
           echo -e "Total written (filtered):\t${bpTotalWritten} bp (${percentBPWritten}%)" >> ~{log}
           echo -e "  Read 1:\t${bpWrittenR1} bp" >> ~{log}
           echo -e "  Read 2:\t${bpWrittenR2} bp" >> ~{log}
 
         else 
           totalRead=$(cat ~{allLog} | grep "Total reads processed:" | cut -d":" -f2 | sed 's/ //g; s/,//g' | awk '{x+=$1}END{print x}')
           adapterR=$(cat ~{allLog} | grep "Reads with adapters:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
           percentAdapterR=$(awk -v A="${adapterR}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')
 
           shortReads=$(cat ~{allLog} | grep "Reads that were too short:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
           percentShortReads=$(awk -v A="${shortReads}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')
 
           ReadsWritten=$(cat ~{allLog} | grep "Reads written (passing filters): " | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
           percentreadsWritten=$(awk -v A="${ReadsWritten}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')                 
 
           echo -e "Total reads processed:\t${totalRead}" >> ~{log}
           echo -e "Reads with adapters:\t${adapterR} (${percentAdapterR}%)" >> ~{log}
           echo -e "Reads that were too short:\t${shortReads} (${percentShortReads}%)" >> ~{log}
           echo -e "Reads written (passing filters):\t${ReadsWritten} (${percentreadsWritten}%)\n\n" >> ~{log}
           echo -e "Total basepairs processed:\t${totalBP} bp" >> ~{log}
           echo -e "Quality-trimmed:\t${bpQualitytrimmed} bp (${percentQualitytrimmed}%)" >> ~{log}
           echo -e "Total written (filtered):\t${bpTotalWritten} bp (${percentBPWritten}%)" >> ~{log}
         fi
     >>>
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
