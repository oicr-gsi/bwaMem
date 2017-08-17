package ca.on.oicr.pde.deciders;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.module.ReturnValue.ExitStatus;
import net.sourceforge.seqware.common.util.Log;

/**
 * This decider finds FastQ files, groups them by IUS Accession, and schedules BWA workflow runs for each pair.
 * The decider performs some validation, and reads file metadata to provide details to the workflow.
 *
 * @author dcooke,mlaszloffy
 *
 */
public class BwaDecider extends OicrDecider {

    private static final String ARG_REFERENCE = "reference";
    private static final String ARG_OUT_FILE = "output-filename";
    private static final String ARG_MANUAL_OUT = "manual-output";
    private static final String ARG_USE_CUTADAPT = "adapter-trimming";
    private static final String ARG_CUTADAPT_MEMORY = "trim-memory";
    private static final String ARG_CUTADAPT_MIN_LENGTH = "trim-min-length";
    private static final String ARG_CUTADAPT_MIN_QUALITY = "trim-min-quality";
    private static final String ARG_CUTADAPT_R1 = "read1-adapter-trim";
    private static final String ARG_CUTADAPT_R2 = "read2-adapter-trim";
    private static final String ARG_CUTADAPT_R1_PARAMS = "read1-trim-params";
    private static final String ARG_CUTADAPT_R2_PARAMS = "read2-trim-params";
    private static final String ARG_BWA_MEMORY = "bwa-memory";

    // Default INI settings here are overridden by command-line parameters in init() and then maintained for all files/groups
    // Details available in metadata are used in preference to both
    // Input
    private String inputReference = null;
    private String inputFile1, inputFile2;

    // Output
    private String outputFileName = "";
    private String iusAccession = "";
    private String groupId = "";
    private String runName = "";
    private String lane = "";
    private String barcode = "NoIndex";
    private boolean manualOutput = false;

    // CutAdapt
    private boolean doForceTrim = false;
    private boolean doTrim = false;
    private int trimMemoryMb = 16384;
    private int trimMinLength = 0;
    private int trimMinQuality = 0;
    private String read1AdapterTrim = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
    private String read2AdapterTrim = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
    private String cutadapt1OtherParams = "";
    private String cutadapt2OtherParams = "";

    // BWA
    private int bwaMemoryMb = 16384;
    
    // Read group header (BWA), pulled from file metadata
    private String rgLibrary = "";
    private String rgPlatform = "";
    private String rgPlatformUnit = "";
    private String rgSample = "";

    private final Map<String, String> modelToPlatform = new HashMap<>();

    public BwaDecider() {
        super();
        // Input args
        parser.accepts(ARG_REFERENCE, "indexed reference genome").withRequiredArg();

        // Output args
        parser.accepts(ARG_OUT_FILE, "Optional: Default filename is created from the IUS accession, library, run name, barcode, and lane.").withRequiredArg();
        parser.accepts(ARG_MANUAL_OUT, "Optional: Set output path manually.");

        // CutAdapt args
        parser.accepts(ARG_USE_CUTADAPT);
        parser.accepts(ARG_CUTADAPT_MEMORY, "Optional: Memory (MB) to allocate for cutadapt (default: 16384)");
        parser.accepts(ARG_CUTADAPT_MIN_LENGTH, "Optional: Minimum length of reads to keep (default: 0)").withRequiredArg();
        parser.accepts(ARG_CUTADAPT_MIN_QUALITY, "Optional: Minimum quality of read ends to keep (default: 0)").withRequiredArg();
        parser.accepts(ARG_CUTADAPT_R1, "Adapter sequence to trim from read 1.").withRequiredArg();
        parser.accepts(ARG_CUTADAPT_R2, "Adapter sequence to trim from read 2.").withRequiredArg();
        parser.accepts(ARG_CUTADAPT_R1_PARAMS, "Optional: Additional cutadapt parameters for read 1.").withRequiredArg();
        parser.accepts(ARG_CUTADAPT_R2_PARAMS, "Optional: Additional cutadapt parameters for read 2.").withRequiredArg();

        // BWA args
        parser.accepts(ARG_BWA_MEMORY, "Optional: Memory (MB) to allocate for bwa mem (default: 16384)");

        //populate model to platform map with known sequencer run models and their corresponding platform (@RG PL)
        //GATK expects one of: ILLUMINA,SLX,SOLEXA,SOLID,454,LS454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN
        modelToPlatform.put("HiSeq", "illumina");
        modelToPlatform.put("ILLUMINA", "illumina");
        modelToPlatform.put("Illumina HiSeq 2500", "illumina");
        modelToPlatform.put("NextSeq 550", "illumina");
        modelToPlatform.put("Illumina MiSeq", "illumina");
        modelToPlatform.put("PacBio RS", "pacbio");
        modelToPlatform.put("454", "454");
        modelToPlatform.put("Genome Analyzer", "illumina");
        modelToPlatform.put("SOLiD", "solid");
        modelToPlatform.put("Illumina_HiSeq_2500", "illumina");
        modelToPlatform.put("NextSeq_550", "illumina");
        modelToPlatform.put("Illumina_MiSeq", "illumina");
    }

    @Override
    public ReturnValue init() {
        // Parse, test, and store command-line parameters. Runs through once per decider invocation (all files)
        Log.debug("INIT");
        this.setMetaType(Arrays.asList("chemical/seq-na-fastq", "chemical/seq-na-fastq-gzip"));
        this.setHeadersToGroupBy(Arrays.asList(Header.IUS_SWA));
        this.setNumberOfFilesPerGroup(2);

        if (this.options.has("group-by")) {
            Log.error("Argument --group-by is not supported");
            System.exit(1);
        }

        // Input
        if (this.options.has(ARG_REFERENCE)) {
            this.inputReference = options.valueOf(ARG_REFERENCE).toString();
        }

        // Output
        if (this.options.has(ARG_OUT_FILE)) {
            outputFileName = options.valueOf(ARG_OUT_FILE).toString();
        }
        if (this.options.has(ARG_MANUAL_OUT)) {
            this.manualOutput = true;
        }

        // CutAdapt
        if (this.options.has(ARG_USE_CUTADAPT)) {
            this.doForceTrim = true;
        }
        if (this.options.has(ARG_CUTADAPT_MEMORY)) {
            this.trimMemoryMb = Integer.valueOf(options.valueOf(ARG_CUTADAPT_MEMORY).toString());
        }
        if (this.options.has(ARG_CUTADAPT_MIN_LENGTH)) {
            this.trimMinLength = Integer.valueOf(options.valueOf(ARG_CUTADAPT_MIN_LENGTH).toString());
        }
        if (this.options.has(ARG_CUTADAPT_MIN_QUALITY)) {
            this.trimMinQuality = Integer.valueOf(options.valueOf(ARG_CUTADAPT_MIN_QUALITY).toString());
        }
        if (this.options.has(ARG_CUTADAPT_R1)) {
            this.read1AdapterTrim = options.valueOf(ARG_CUTADAPT_R1).toString();
        }
        if (this.options.has(ARG_CUTADAPT_R2)) {
            this.read2AdapterTrim = options.valueOf(ARG_CUTADAPT_R2).toString();
        }
        if (this.options.has(ARG_CUTADAPT_R1_PARAMS)) {
            this.cutadapt1OtherParams = options.valueOf(ARG_CUTADAPT_R1_PARAMS).toString();
        }
        if (this.options.has(ARG_CUTADAPT_R2_PARAMS)) {
            this.cutadapt2OtherParams = options.valueOf(ARG_CUTADAPT_R2_PARAMS).toString();
        }

        // BWA
        if (this.options.has(ARG_BWA_MEMORY)) {
            this.bwaMemoryMb = Integer.valueOf(options.valueOf(ARG_BWA_MEMORY).toString());
        }

        return super.init();
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        // Validate individual files and read metadata
        Log.debug("CHECK FILE DETAILS:" + fm);
        if (!super.checkFileDetails(returnValue, fm)) {
            return false;
        }

        FileAttributes attribs = new FileAttributes(returnValue, returnValue.getFiles().get(0));

        //  Skip if library_source_template_type isn't WG, EX, or TS
        String filePath = fm.getFilePath();
        String templateType = attribs.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE);
        if (!"WG".equals(templateType) && !"EX".equals(templateType) && !"TS".equals(templateType)) {
            Log.debug("Skipping " + filePath + " due to incompatible library template type " + templateType);
            return false;
        }
        this.doTrim = this.doForceTrim;
        String doTrimValue = String.valueOf(this.doTrim);
        if (doTrimValue.equals("false")) {
            if (templateType.equals("EX") || templateType.equals("TS")) {
                this.doTrim = true;
            }
        }
        // If xenograft, check if it's a Xenome output
        String tissueType = attribs.getLimsValue(Lims.TISSUE_TYPE);
        if (tissueType == null) {
            Log.debug("Skipping " + filePath + " due to missing Tissue Type");
            return false;
        } else if (tissueType.equals("X") && !filePath.contains("xenome")) {
            Log.debug("Skipping " + filePath + " because it is not a Xenome output");
            return false;
        }

        // Get additional metadata
        this.iusAccession = returnValue.getAttribute(Header.IUS_SWA.getTitle());
        this.groupId = attribs.getLimsValue(Lims.GROUP_ID);
        if (groupId == null) {
            groupId = "";
        }
        this.runName = returnValue.getAttribute(Header.SEQUENCER_RUN_NAME.getTitle());
        this.lane = returnValue.getAttribute(Header.LANE_NUM.getTitle());

        this.barcode = attribs.getBarcode();
        this.rgLibrary = getRGLB(attribs);
        this.rgPlatform = returnValue.getAttribute("Sequencer Run Platform Name");
        this.rgSample = getRGSM(attribs);
        this.rgPlatformUnit = this.runName
                + "-"
                + this.barcode
                + "_"
                + this.lane;

        return true;
    }

    /**
     * Constructs a String for use in the SAM read group header SM field
     *
     * @param fa metadata for the sample file
     *
     * @return a String in the format: {donor}_{tissue origin}_{tissue type}[_group id]
     */
    private String getRGSM(FileAttributes fa) {
        String groupId = fa.getLimsValue(Lims.GROUP_ID);

        StringBuilder sb = new StringBuilder()
                .append(fa.getDonor())
                .append("_")
                .append(fa.getLimsValue(Lims.TISSUE_ORIGIN))
                .append("_")
                .append(fa.getLimsValue(Lims.TISSUE_TYPE));

        if (groupId != null) {
            sb.append("_").append(groupId);
        }

        return sb.toString();
    }
    
    /**
     * Constructs a String for use in the SAM read group header LB field
     *
     * @param fa metadata for the sample file
     *
     * @return a String in the format: {library}[_group id]
     */
    private String getRGLB(FileAttributes fa) {
        String groupId = fa.getLimsValue(Lims.GROUP_ID);

        StringBuilder sb = new StringBuilder()
                .append(fa.getLibrarySample());

        if (groupId != null) {
            sb.append("_").append(groupId);
        }

        return sb.toString();
    }

    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        // Make sure there are two appropriate files for paired end
        this.inputFile1 = null;
        this.inputFile2 = null;

        if (commaSeparatedFilePaths.contains(",")) {
            String[] fqFilesArray = commaSeparatedFilePaths.split(",");

            for (String file : fqFilesArray) {
                int mate = idMate(file);
                switch (mate) {
                    case 1:
                        if (this.inputFile1 != null) {
                            Log.error("More than one file found for read 1: " + inputFile1 + ", " + file);
                            return new ReturnValue(ExitStatus.INVALIDFILE);
                        }
                        this.inputFile1 = file;
                        break;
                    case 2:
                        if (this.inputFile2 != null) {
                            Log.error("More than one file found for read 2: " + inputFile2 + ", " + file);
                            return new ReturnValue(ExitStatus.INVALIDFILE);
                        }
                        this.inputFile2 = file;
                        break;
                    default:
                        Log.error("Cannot identify " + file + " end (read 1 or 2)");
                        return new ReturnValue(ExitStatus.INVALIDFILE);
                }
            }
        } else {
            Log.error("Missing file: Two files required for workflow");
            return new ReturnValue(ExitStatus.INVALIDFILE);
        }

        if (modelToPlatform.containsKey(this.rgPlatform)) {
            this.rgPlatform = modelToPlatform.get(this.rgPlatform);
        } else {
            Log.warn("Sequencer run model = [" + this.rgPlatform + "] platform is missing");
            return new ReturnValue(ExitStatus.FAILURE);
        }

        return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        Log.debug("INI FILE:" + commaSeparatedFilePaths);

        Map<String, String> iniFileMap = super.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);

        // Inputs
        iniFileMap.put("input_file_1", this.inputFile1);
        iniFileMap.put("input_file_2", this.inputFile2);
        if (this.inputReference != null) {
            iniFileMap.put("input_reference", this.inputReference);
        }

        // Output
        iniFileMap.put("outputFileName", this.outputFileName);
        iniFileMap.put("ius_accession", this.iusAccession);
        iniFileMap.put("group_id", this.groupId);
        iniFileMap.put("sequencer_run_name", this.runName);
        iniFileMap.put("barcode", this.barcode);
        iniFileMap.put("lane", this.lane);

        iniFileMap.put("manual_output", String.valueOf(this.manualOutput));

        // CutAdapt
        iniFileMap.put("do_trim", String.valueOf(this.doTrim));
        iniFileMap.put("trim_mem_mb", String.valueOf(this.trimMemoryMb));
        iniFileMap.put("trim_min_quality", String.valueOf(this.trimMinQuality));
        iniFileMap.put("trim_min_length", String.valueOf(this.trimMinLength));
        iniFileMap.put("r1_adapter_trim", this.read1AdapterTrim);
        iniFileMap.put("r2_adapter_trim", this.read2AdapterTrim);
        iniFileMap.put("cutadapt_r1_other_params", this.cutadapt1OtherParams);
        iniFileMap.put("cutadapt_r2_other_params", this.cutadapt2OtherParams);

        // BWA
        iniFileMap.put("bwa_aln_mem_mb", String.valueOf(this.bwaMemoryMb));

        // Read Group Header (added by BWA)
        iniFileMap.put("rg_library", this.rgLibrary);
        iniFileMap.put("rg_platform", this.rgPlatform);
        iniFileMap.put("rg_platform_unit", this.rgPlatformUnit);
        iniFileMap.put("rg_sample_name", this.rgSample);

        return iniFileMap;
    }

    public static void main(String args[]) {
        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(BwaDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }
}
