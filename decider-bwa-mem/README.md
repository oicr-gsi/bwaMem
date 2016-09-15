##BWA-MEM Decider

Version 1.1, SeqWare version 1.1.0

###Overview

This decider launches the [BWA-MEM Workflow](../workflow-bwa-mem) on paired-end FastQ files. Pairs are identified by IUS accession, and the decider relies on hints in the filename to determine whether a file is read 1 or read 2. File details including IUS, run name, lane, barcode, library, platform, sample name, and group ID are pulled from the metadata to include in the generated INI files.

Please see [basic deciders](https://seqware.github.io/docs/6-pipeline/basic_deciders) for additional information.

###Validation and Filtering

The following are causes for the decider to skip passing a file on to a workflow run.
* Library template type: The decider only accepts files with library template types WG, EX, and TS
* Missing tissue type value: Every fastq file should have a tissue type in LIMS \('n' for unknown is acceptable\)
* Tissue type is Xenograph, but the file doesn't appear to be a Xenome output \(the decider checks for "xenome" in the filename\)
* Read 1 or read 2 missing: The workflow is designed for aligning paired end reads, and requires two files

###SAM Read Group Header

Using file metadata, the decider will set read group header values in this format:
* ID: runName-barcode_lane
* LB: librarySample
* PL: sequencerRunPlatformName
* PU: \(same as ID\)
* SM: donor_tissueOrigin_tissueType\[_groupId\]

Example result:

```
'@RG    ID:130110_SN804_0104_AD1NPCACXX-NoIndex_6    LB:ASHPC_0005_Pa_R_PE_700_WG    PL:ILLUMINA    PU:130110_SN804_0104_AD1NPCACXX-NoIndex_6    SM:PCSI_0350_Pa_n_2'
```

###Compile

```
mvn clean install
```

###Usage
```
java -jar Decider.jar --study-name \<study-name\> --wf-accession \<bwa-mem-workflow-accession\>
```

###Options

**Required**
See [basic deciders](https://seqware.github.io/docs/6-pipeline/basic_deciders) for general decider options. No additional options are strictly required.

**Optional**

Parameter | Type | Description \[default\]
----------|------|-------------
verbose | none | Log all SeqWare (debug) information
outputPath | path | Path to store the file output files \[./\]
output-filename | string | Specific filename to use
manual-output | none | Set output path manually
output-format | string | Alignment output format. May be SAM, BAM, or CRAM \[BAM\]
adapter-trimming | none | Enable to trim adapters
trim-memory | int | RAM in MB to allocate for CutAdapt job \[16384\]
trim-min-length | int | Minimum length of reads to keep \[0\]
trim-min-quality | int | Minimum quality of read ends to keep \[0\]
read1-adapter-trim | string | Adapter sequence to trim from read 1 \[AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG\]
read2-adapter-trim | string | Adapter sequence to trim from read 2 \[AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\]
read1-trim-params | string | Additional cutadapt parameters for read 1
read2-trim-params | string | Additional cutadapt parameters for read 2
bwa-memory | int | RAM in MB to allocate for align/sort/convert job \[16384\]
bwa-threads | int | Threads to use for BWA-MEM \[8\]
bwa-pacbio | none | Enable BWA PacBio mode
bwa-ont2d | none | Enable BWA ONT mode
bwa-no-mark-secondary | none | Disable marking of supplementary alignments as secondary. This will break compatibility with Picard
bwa-params | string | Additional BWA-MEM parameters
samtools-memory | int | RAM in MB to allocate for Samtools index job if output format is BAM or CRAM

##Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .
