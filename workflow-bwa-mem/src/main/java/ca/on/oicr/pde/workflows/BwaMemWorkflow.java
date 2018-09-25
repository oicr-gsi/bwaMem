package ca.on.oicr.pde.workflows;

import java.util.Map;

import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;

/**
 * This workflow takes paired-end FastQ input files and performs alignment, sorting, conversion,
 * and indexing. It makes use of CutAdapt, BWA-MEM, and Samtools, to build the output SAM, BAM with index,
 * or CRAM with index.
 *
 * @author dcooke
 *
 */
public class BwaMemWorkflow extends OicrWorkflow {

    private static final String FASTQ_MIMETYPE = "chemical/seq-na-fastq";
    private static final String FASTQ_GZIP_MIMETYPE = "chemical/seq-na-fastq-gzip";

    private static enum AlignmentFormat {
        SAM(".sam", "application/sam", null, null),
        BAM(".bam", "application/bam", ".bai", "application/bam-index"),
        CRAM(".cram", "application/cram", ".crai", "application/cram-index");

        private final String fileExtension;
        private final String mimeType;
        private final String indexExtension;
        private final String indexType;

        AlignmentFormat(String fileExtension, String mimeType, String indexExtension, String indexType) {
            this.fileExtension = fileExtension;
            this.mimeType = mimeType;
            this.indexExtension = indexExtension;
            this.indexType = indexType;
        }

        public String getFileExtension() {
            return fileExtension;
        }

        public String getMimeType() {
            return mimeType;
        }

        public String getIndexExtension() {
            return indexExtension;
        }

        public String getIndexType() {
            return indexType;
        }
    }

    private String dataDir = null;
    private String tempDir = null;
    private String reference_path = null;
    private SqwFile read1;
    private SqwFile read2;
    private String bwa, samtools, jre;
    private String outputFilePath = null;
    private String outputIndexPath = null;
    private String cutadaptCmd = null;
    private SqwFile cutadaptLogFile;
    private SqwFile outputFile;
    private SqwFile outputIndex;
    private AlignmentFormat outputFormat;

    private String queue;

    private void init() {
        final String binDir = getWorkflowBaseDir() + "/bin/";
        bwa = binDir + getProperty("bwa");
        samtools = binDir + getProperty("samtools");
        jre = binDir + getProperty("bundled_jre") + "/bin/java";
        outputFormat = AlignmentFormat.valueOf(getProperty("output_format"));
        // initialize cutadaptCmd which is a command that all cutadapt jobs will use
        if (getPropertyOrNull("module") != null) {
            cutadaptCmd = getProperty("module") + "; ";
        }
        if (cutadaptCmd != null) {
            cutadaptCmd += getProperty("python") + " " + getProperty("cutadapt") + " -q " + getProperty("trim_min_quality") + " -m " + getProperty("trim_min_length");
        } else {
            cutadaptCmd = getProperty("python") + " " + getProperty("cutadapt") + " -q " + getProperty("trim_min_quality") + " -m " + getProperty("trim_min_length");
        }
    }

    @Override
    public void setupDirectory() {
        // this is first method called, so init from here
        init();

        // create the final output directory
        dataDir = getProperty("data_dir");
        if (!dataDir.endsWith("/")) {
            dataDir += "/";
        }
        this.addDirectory(dataDir);

        tempDir = getProperty("tmp_dir");
        if (!tempDir.endsWith("/")) {
            tempDir += "/";
        }
        this.addDirectory(tempDir);
    }

    @Override
    public Map<String, SqwFile> setupFiles() {
        String input1_path = null;
        String input2_path = null;

        input1_path = getProperty("input_file_1");
        input2_path = getPropertyOrNull("input_file_2");
        reference_path = getProperty("input_reference");

        String outputFileName = getOutputFileName();

        outputFilePath = dataDir + outputFileName;

        // register input files
        read1 = this.createFile("file_in_0");
        read1.setSourcePath(input1_path);
        read1.setType(input1_path.endsWith(".gz") ? FASTQ_GZIP_MIMETYPE : FASTQ_MIMETYPE);
        read1.setIsInput(true);

        if (input2_path != null) {
            read2 = this.createFile("file_in_1");
            read2.setSourcePath(input2_path);
            read2.setType(input1_path.endsWith(".gz") ? FASTQ_GZIP_MIMETYPE : FASTQ_MIMETYPE);
            read2.setIsInput(true);
        }

        // register output files
        outputFile = createOutputFile(outputFilePath, outputFormat.getMimeType(), Boolean.valueOf(getProperty("manual_output")));
        if (outputFormat.getIndexType() != null) {
            outputIndexPath = outputFilePath + outputFormat.getIndexExtension();
            outputIndex = createOutputFile(outputIndexPath, outputFormat.getIndexType(), Boolean.valueOf(getProperty("manual_output")));
        }

        return this.getFiles();
    }

    private String getOutputFileName() {
        String outputFileName = getPropertyOrNull("output_file_name");
        if (outputFileName == null) {
            String groupId = getPropertyOrNull("group_id");
            StringBuilder sb = new StringBuilder()
                    .append("SWID_")
                    .append(getProperty("ius_accession"))
                    .append("_")
                    .append(getProperty("rg_library"))
                    .append("_");

            if (groupId != null) {
                sb.append(groupId).append("_");
            }

            sb.append(getProperty("sequencer_run_name"))
                    .append("_")
                    .append(getProperty("barcode"))
                    .append("_L00")
                    .append(getProperty("lane"))
                    .append("_001.annotated")
                    .append(outputFormat.getFileExtension());
            outputFileName = sb.toString();
        }

        return outputFileName;

    }

    @Override
    public void buildWorkflow() {
        queue = getPropertyOrNull("queue");
        String r1 = null;
        String r2 = null;
        String basename1 = null;
        String basename2 = null;

        r1 = read1.getProvisionedPath();
        basename1 = r1.substring(r1.lastIndexOf("/") + 1, r1.lastIndexOf(".fastq.gz"));

        if (read2 != null) {
            r2 = read2.getProvisionedPath();
            basename2 = r2.substring(r2.lastIndexOf("/") + 1, r2.lastIndexOf(".fastq.gz"));
        }

        // cutadapt - save log files
        Job trimLogJob01 = getCutadaptLogFile(r1, "r1", basename1);

        if (queue != null) {
            trimLogJob01.setQueue(queue);
        }

        if (read2 != null) {
            Job trimLogJob02 = null;
            trimLogJob02 = getCutadaptLogFile(r2, "r2", basename2);
            if (queue != null) {
                trimLogJob02.setQueue(queue);
            }
        }

        // cutadapt (optional) trim reads and pass trimmed output files to align job
        Job trimJob01 = null;
        if (Boolean.valueOf(getProperty("do_trim"))) {
            String trimmed1 = this.dataDir + basename1 + ".trim.fastq.gz";
            if (read2 == null) {
                trimJob01 = getCutadaptJob(r1, null, trimmed1, null);
                r1 = trimmed1;
            } else {
                String trimmed2 = this.dataDir + basename2 + ".trim.fastq.gz";
                trimJob01 = getCutadaptJob(r1, r2, trimmed1, trimmed2);
                r1 = trimmed1;
                r2 = trimmed2;
            }
            if (queue != null) {
                trimJob01.setQueue(queue);
            }
        }

        // Align, sort, and convert
        Job alignJob = getAlignJob(r1, r2);
        if (trimJob01 != null) {
            alignJob.addParent(trimJob01);
        }
        if (queue != null) {
            alignJob.setQueue(queue);
        }

        // Index
        if (outputFormat != AlignmentFormat.SAM) {
            Job indexJob = getIndexJob();
            indexJob.addParent(alignJob);
            if (queue != null) {
                indexJob.setQueue(queue);
            }
        }
    }

    /**
     * Creates a job to trim adapters from a read (either read1 or read2) using cutadapt, discard trimmed output but save log files
     *
     * @param readPath  input file
     * @param whichRead r1 or r2
     * @param basename  basename in the output path
     *
     * @return
     */
    private Job getCutadaptLogFile(String readPath, String whichRead, String basename) {
        Job job = this.getWorkflow().createBashJob("cutadaptLog");
        job.setMaxMemory(getProperty("cutadapt_log_mem_mb"));

        String logFilePath = this.dataDir + basename + ".log";
        cutadaptLogFile = createOutputFile(logFilePath, "text/plain", Boolean.valueOf(getProperty("manual_output")));
        cutadaptLogFile.getAnnotations().put("tool", getProperty("cutadapt_log_file_tool_annotation"));
        cutadaptLogFile.getAnnotations().put("type", getProperty("cutadapt_log_file_type_annotation"));

        Command cmd = job.getCommand();
        cmd.addArgument(cutadaptCmd);
        if (whichRead == "r1") {
            cmd.addArgument(" -a " + getProperty("r1_adapter_trim") + " " + getProperty("cutadapt_r1_log_other_params"));
        } else {
            cmd.addArgument(" -a " + getProperty("r2_adapter_trim") + " " + getProperty("cutadapt_r2_log_other_params"));
        }
        cmd.addArgument(" -o " + "/dev/null" + " " + readPath);
        cmd.addArgument("> " + logFilePath);
        job.addFile(cutadaptLogFile);

        return job;
    }

    /**
     * Creates a job to trim adapters from the reads using cutadapt, save trimmed output files
     * If read2 is not null, pass read1 and read2 together in one command so cutadapt will check if they match as expected
     *
     * @param read1Path input file
     * @param read2Path input file
     * @param trimmed1  output file path for read1
     * @param trimmed2  output file path for read2
     *
     * @return
     */
    private Job getCutadaptJob(String read1Path, String read2Path, String trimmed1, String trimmed2) {
        Job job = this.getWorkflow().createBashJob("cutadapt");
        job.setMaxMemory(getProperty("trim_mem_mb"));
        Command cmd = job.getCommand();
        cmd.addArgument(cutadaptCmd + " " + getProperty("cutadapt_additional_params") + " -a " + getProperty("r1_adapter_trim"));
        if (read2Path == null) {
            cmd.addArgument(" -o " + trimmed1 + " " + read1Path);
        } else {
            // -A parameter is only available in cutadapt v1.8+ 
            cmd.addArgument(" -A " + getProperty("r2_adapter_trim"));
            cmd.addArgument(" -o " + trimmed1 + " -p " + trimmed2 + " " + read1Path + " " + read2Path);
        }

        return job;
    }

    /**
     * Creates a job to run alignment, sorting, and conversion to the desired output format.
     * Runs bwa mem | samtools sort | samtools view (if necessary)
     *
     * @param read1
     * @param read2
     *
     * @return the job
     */
    private Job getAlignJob(String read1, String read2) {
        Job job = this.getWorkflow().createBashJob("bwa_mem");
        job.addFile(outputFile);
        job.setMaxMemory(getProperty("bwa_mem_mb"));

        Command cmd = job.getCommand();
        cmd.addArgument("set -e; set -o pipefail;");
        cmd.addArgument(bwa + " mem");
        if (Boolean.valueOf(getProperty("bwa_mark_secondary_alignments"))) {
            cmd.addArgument("-M");
        }
        cmd.addArgument(getBwaSpecialMode());
        String threads = getPropertyOrNull("bwa_threads");
        if (threads != null) {
            cmd.addArgument("-t " + threads);
        }
        String otherParams = getPropertyOrNull("bwa_other_params");
        if (otherParams != null) {
            cmd.addArgument(otherParams);
        }
        cmd.addArgument("-R " + getReadGroupHeader());
        cmd.addArgument(reference_path);
        cmd.addArgument(read1);
        if (read2 != null) {
            cmd.addArgument(read2);
        }

        // pipe to samtools sort
        cmd.addArgument("|");
        cmd.addArgument(samtools + " sort");

        switch (outputFormat) {
            case SAM:
                cmd.addArgument("-O sam -T " + this.tempDir);
                break;
            case BAM:
                cmd.addArgument("-O bam -T " + this.tempDir);
                break;
            case CRAM:
                cmd.addArgument("-O bam -l 0 -T " + this.tempDir + " -"); // uncompressed bam (samtools sort cannot output to cram)

                // pipe to samtools view to convert to CRAM
                cmd.addArgument("|");
                cmd.addArgument(samtools + " view");
                cmd.addArgument("-T " + reference_path);
                cmd.addArgument("-C");
                break;
        }
        // Either ending - samtools sort or samtools view - requires output path and input source (stdin)
        cmd.addArgument("-o " + outputFilePath);
        cmd.addArgument("-");

        return job;
    }

    /**
     * Creates a job to index the output BAM or CRAM
     *
     * @return the job
     */
    private Job getIndexJob() {
        Job job = this.getWorkflow().createBashJob("samtools_index");

        Command cmd = job.getCommand();
        cmd.addArgument(samtools + " index");
        cmd.addArgument(outputFilePath);

        job.addFile(outputIndex);
        job.setMaxMemory(getProperty("samtools_index_mem_mb"));
        return job;
    }

    /**
     * Checks settings for BWA modes
     *
     * @return the BWA parameter to append (e.g. "-x pacbio" if pacbio mode is set to true) or an empty String if none
     */
    private String getBwaSpecialMode() {
        String modeParam = null;
        if (Boolean.valueOf(getPropertyOrNull("bwa_pacbio_mode"))) {
            modeParam = "-x pacbio";
        }
        if (Boolean.valueOf(getPropertyOrNull("bwa_ont_mode"))) {
            if (modeParam != null) {
                throw new IllegalArgumentException("bwa_pacbio_mode and bwa_ont_mode cannot both be used.");
            }
            modeParam = "-x ont2d";
        }
        return modeParam == null ? "" : modeParam;
    }

    /**
     * Creates the SAM read group header, using values from the INI.
     *
     * @return a String formatted similar to "'@RG ID:? LB:? PL:? PU:? SM:?'" with elements separated by tabs
     */
    private String getReadGroupHeader() {
        String unit = getProperty("rg_platform_unit");

        StringBuilder sb = new StringBuilder();
        sb.append("'@RG\\tID:");
        sb.append(unit);
        sb.append("\\tLB:");
        sb.append(getProperty("rg_library"));
        sb.append("\\tPL:");
        sb.append(getProperty("rg_platform"));
        sb.append("\\tPU:");
        sb.append(unit);
        sb.append("\\tSM:");
        sb.append(getProperty("rg_sample_name"));
        sb.append("'");

        return sb.toString();
    }

    /**
     * Convenience method. Checks for an optional property in the INI and returns either the non-empty
     * String value, or null.
     *
     * @return the value if it exists and is non-empty; null otherwise
     */
    private String getPropertyOrNull(String key) {
        if (!hasPropertyAndNotNull(key)) {
            return null;
        }
        String val = getProperty(key);
        return val.isEmpty() ? null : val;
    }

}
