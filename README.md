# bwaMem

BwaMem Workflow version 2.0

## Overview

## Dependencies

* [bwa 0.7.12](https://github.com/lh3/bwa/archive/0.7.12.tar.gz)
* [samtools 1.9](https://github.com/samtools/samtools/archive/0.1.19.tar.gz)
* [cutadapt 1.8.3](https://cutadapt.readthedocs.io/en/v1.8.3/)
* [slicer 0.3.0](https://github.com/OpenGene/slicer/archive/v0.3.0.tar.gz)


## Usage

### Cromwell
```
java -jar cromwell.jar run bwaMem.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`fastqR1`|File|fastq file for read 1
`readGroups`|String|Complete read group header line
`runBwaMem.modules`|String|Required environment modules
`runBwaMem.bwaRef`|String|The reference genome to align the sample with by BWA


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqR2`|File?|None|fastq file for read 2
`outputFileNamePrefix`|String|"output"|Prefix for output file
`numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`doTrim`|Boolean|false|if true, adapters will be trimmed before alignment
`trimMinLength`|Int|1|minimum length of reads to keep [1]
`trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]


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
`adapterTrimming.modules`|String|"cutadapt/1.8.3"|Required environment modules
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
 
 This section lists command(s) run by fastqc workflow
 
 * Running BWA
 
 bwaMem workflow runs the following command (excerpt from .wdl file).
 
 * countChunkSize
 
 Depending on the file size job will be scattered in a number of chunks.
 
 ```
 totalLines=$(zcat FASTQ_FILE | wc -l)
 python -c "from math import ceil; print int(ceil(($totalLines/4.0)/~{numChunk})*4)"
 ```
 
 * slicer
 
 CHUNK_SIZE calcualted in the previous step.
 
 ```
  slicer -i FASTQ_FILE -l CHUNK_SIZE --gzip
 
 ```
 
 * adapterTrimming
 
 Optional adapter trimming
 
 ```
 cutadapt -q TRIM_MIN_QUALITY
             -m TRIM_MIN_LENGTH
             -a ADAPTER_1
             -o FASTQ_FILER1_TRIMMED
             -A ADAPTER_2 -p FATSQ_FILER2_TRIMMED
             ADDITIONAL_PARAM
             FASTQ_R1
             FASTQ_R2 > RESULT_LOG
 ```
 
 * run BwaMem
 
 Main command, need to pass READ_GROUP_STRING (Read group information) and BWA_REF (directory with BWA reference files)
 
 ```
 bwa mem -M
         -t THEADS ADDITIONAL_PARAM
         -R  READ_GROUP_STRING
             BWA_REF
             FASTQ_R1
             FASTQ_R2
         |
         samtools sort -O bam -T TMP_DIR -o RESULT_BAM -
 ```
 
 * bamMerge
 
 Merging chunks into final bam file.
 
 ```
 samtools merge
         -c
         RESULT_MERGE_BAM
         BAM_FILE1 BAM_FILE2 ...
 
 ```
 
 * indexBam
 
 Indexing with samtools, name of the index file is specified with INDEX_BAI
 
 ```
  samtools index INPUT_BAM INDEX_BAI
 
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
