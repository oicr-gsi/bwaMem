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
`fastqR2`|File|fastq file for read 2
`readGroups`|String|Complete read group header line
`runBwaMem.modules`|String|Required environment modules
`runBwaMem.bwaRef`|String|The reference genome to align the sample with by BWA


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
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
`slicerR1.modules`|String|"slicer/0.3.0"|
`slicerR1.jobMemory`|Int|16|
`slicerR1.timeout`|Int|48|
`slicerR2.modules`|String|"slicer/0.3.0"|
`slicerR2.jobMemory`|Int|16|
`slicerR2.timeout`|Int|48|
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


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with wdl_doc_gen (https://github.com/oicr-gsi/wdl_doc_gen/)_
