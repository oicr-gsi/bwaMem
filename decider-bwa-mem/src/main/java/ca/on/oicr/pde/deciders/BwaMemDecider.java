package ca.on.oicr.pde.deciders;

import ca.on.oicr.pde.deciders.FileAttributes;
import ca.on.oicr.pde.deciders.Lims;
import ca.on.oicr.pde.deciders.OicrDecider;

import java.util.*;

import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.module.ReturnValue.ExitStatus;
import net.sourceforge.seqware.common.util.Log;

/**
 * This decider finds FastQ files, groups them by IUS Accession, and schedules BWA-MEM workflow runs for each pair. 
 * The decider performs some validation, and reads file metadata to provide details to the workflow.
 * 
 * @author dcooke
 *
 */
public class BwaMemDecider extends OicrDecider {
	
	private static enum AlignmentFormat {SAM, BAM, CRAM};
	
	// Default INI settings here are overridden by command-line parameters in init() and then maintained for all files/groups
	// Details available in metadata are used in preference to both
	
	// Input
	private String inputReference = null;
	private String inputFile1, inputFile2;
	
	// Output
	private String outputFileName = "";
	private AlignmentFormat outputFormat = AlignmentFormat.BAM;
	private String iusAccession = "";
	private String runName = "";
	private String lane = "";
	private String barcode = "NoIndex";
	private boolean manualOutput = false;
	
	// CutAdapt
	private boolean doTrim = false;
	private int trimMemoryMb = 16384;
	private String read1AdapterTrim = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG";
	private String read2AdapterTrim = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
	private String cutadapt1OtherParams = "";
	private String cutadapt2OtherParams = "";
	
	// BWA
	private int bwaMemoryMb = 16384;
	private int bwaThreads = 8;
	private boolean bwaPacbioMode = false;
	private boolean bwaOntMode = false;
	private String bwaAdditionalParams = "";
	
	// samtools index
	private int samtoolsMemoryMb = 3072;
	
	// Read group header (BWA), pulled from file metadata
	private String rgLibrary = "";
	private String rgPlatform = "";
	private String rgPlatformUnit = "";
	private String rgSample = "";
	
	
	public BwaMemDecider() {
		super();
		// Input args
		parser.accepts("reference", "indexed reference genome").withRequiredArg();
        
		// Output args
		parser.accepts("verbose", "Optional: Log all SeqWare information.");
		parser.accepts("output-filename", "Optional: Default filename is created from the IUS accession, library, run name, barcode, and lane.").withRequiredArg();
		parser.accepts("manual-output", "Optional: Set output path manually.");
		parser.accepts("output-format", "Optional: Alignment output format. Options are SAM, BAM (default), and CRAM.").withRequiredArg();
		
		// CutAdapt args
		parser.accepts("adapter-trimming");
		parser.accepts("trim-memory", "Optional: Memory (MB) to allocate for cutadapt (default: 16384)");
		parser.accepts("read1-adapter-trim", "Adapter sequence to trim from read 1.").withRequiredArg();
		parser.accepts("read2-adapter-trim", "Adapter sequence to trim from read 2.").withRequiredArg();
		parser.accepts("read1-trim-params", "Optional: Additional cutadapt parameters for read 1.").withRequiredArg();
		parser.accepts("read2-trim-params", "Optional: Additional cutadapt parameters for read 2.").withRequiredArg();
		
		// BWA args
		parser.accepts("bwa-memory", "Optional: Memory (MB) to allocate for bwa mem (default: 16384)");
		parser.accepts("bwa-threads", "Optional: Threads to use for bwa mem (default: 8).").withRequiredArg();
		parser.accepts("bwa-pacbio", "Optional: Execute BWA in PacBio mode.");
		parser.accepts("bwa-ont2d", "Optional: Execute BWA in ONT mode.");
		parser.accepts("bwa-params", "Optional: Additional bwa mem parameters.").withRequiredArg();
		
		// Samtools index args
		parser.accepts("samtools-memory", "Optional: Memory (MB) to allocate for Samtools index (default: 16384)");
	}
	
	@Override
	public ReturnValue init() {
		// Parse, test, and store command-line parameters. Runs through once per decider invocation (all files)
		if (this.options.has("verbose")) {
			Log.setVerbose(true);
		}
		Log.debug("INIT");
		this.setMetaType(Arrays.asList("chemical/seq-na-fastq", "chemical/seq-na-fastq-gzip"));
		
		if (this.options.has("group-by")) {
            Log.error("Grouping by anything other than IUS_SWA does not make much sense for this workflow, ignoring...");
        }
		
		// Input
		if (this.options.has("reference")){
			this.inputReference = options.valueOf("reference").toString();
		}
		
		// Output
		if (this.options.has("output-filename")) {
			outputFileName = options.valueOf("output-filename").toString();
		}
		if (this.options.has("manual-output")) {
			this.manualOutput = true;
		}
		if (this.options.has("output-format")) {
			outputFormat = AlignmentFormat.valueOf(options.valueOf("output-format").toString());
		}
		
		// CutAdapt
		if (this.options.has("adapter-trimming")) {
			this.doTrim = true;
		}
		if (this.options.has("trim-memory")) {
			this.trimMemoryMb = Integer.valueOf(options.valueOf("trim-memory").toString());
		}
		if (this.options.has("read1-adapter-trim")) {
			this.read1AdapterTrim = options.valueOf("read1-adapter-trim").toString();
		}
		if (this.options.has("read2-adapter-trim")) {
			this.read2AdapterTrim = options.valueOf("read2-adapter-trim").toString();
		}
		if (this.options.has("read1-trim-params")) {
			this.cutadapt1OtherParams = options.valueOf("read1-trim-params").toString();
		}
		if (this.options.has("read2-trim-params")) {
			this.cutadapt2OtherParams = options.valueOf("read2-trim-params").toString();
		}
		
		// BWA
		if (this.options.has("bwa-memory")) {
			this.bwaMemoryMb = Integer.valueOf(options.valueOf("bwa-memory").toString());
		}
		if (this.options.has("bwa-threads")) {
			this.bwaThreads = Integer.valueOf(options.valueOf("bwa-threads").toString());
		}
		if (this.options.has("bwa-pacbio")) {
			this.bwaPacbioMode = Boolean.valueOf(options.valueOf("bwa-pacbio").toString());
		}
		if (this.options.has("bwa-ont2d")) {
			this.bwaOntMode = Boolean.valueOf(options.valueOf("bwa-ont2d").toString());
		}
		if (this.options.has("bwa-params")) {
			this.bwaAdditionalParams = options.valueOf("bwa-params").toString();
		}
		
		// Samtools index
		if (this.options.has("samtools-memory")) {
			this.samtoolsMemoryMb = Integer.valueOf(options.valueOf("samtools-memory").toString());
		}
		
		return super.init();
	}
	
	@Override
	public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
		Log.debug("SEPARATE_FILES: "+vals.size()+"files, groupBy="+groupBy);
		this.setNumberOfFilesPerGroup(2);
		
		// Group by IUS Accession
		Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();
		for (ReturnValue rv : vals) {
			String ius = rv.getAttribute(Header.IUS_SWA.getTitle());
			List<ReturnValue> list = map.get(ius);
			if (list == null) list = new ArrayList<>();
			list.add(rv);
			map.put(ius, list);
		}
		return map;
	}

	@Override
	protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
		// Validate individual files and read metadata
		Log.debug("CHECK FILE DETAILS:" + fm);
		if (!super.checkFileDetails(returnValue, fm)) return false;
		
		FileAttributes attribs = new FileAttributes(returnValue, returnValue.getFiles().get(0));
		
		// If xenograft, check if it's a Xenome output
		String filePath = fm.getFilePath();
		if (attribs.getLimsValue(Lims.TISSUE_TYPE).equals("X") && !filePath.contains("xenome")) {
			Log.debug("Skipping "+filePath+" because it is not a Xenome output");
			return false;
		}
		
		//Get additional metadata
		this.iusAccession = returnValue.getAttribute(Header.IUS_SWA.getTitle());
		this.runName = returnValue.getAttribute(Header.SEQUENCER_RUN_NAME.getTitle());
		this.lane = returnValue.getAttribute(Header.LANE_NUM.getTitle());
		
		this.barcode = attribs.getBarcode();
		this.rgLibrary = attribs.getLibrarySample();
		this.rgPlatform = returnValue.getAttribute("Sequencer Run Platform Name");
		this.rgSample = getRGSM(attribs);
		this.rgPlatformUnit = this.runName 
				+ "-" 
				+ this.barcode
				+ "_" 
				+ this.lane;
		
		return true;
	}
	
	private String getRGSM(FileAttributes fa) {
		StringBuilder sb = new StringBuilder()
				.append(fa.getDonor())
				.append("_")
				.append(fa.getLimsValue(Lims.TISSUE_ORIGIN))
				.append("_")
				.append(fa.getLimsValue(Lims.TISSUE_TYPE));
		
		String groupId = fa.getLimsValue(Lims.GROUP_ID);
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
			String [] fqFilesArray = commaSeparatedFilePaths.split(",");
			
			for (String file : fqFilesArray) {
				int mate = idMate(file);
				switch (mate) {
				case 1:
					if (this.inputFile1 != null) {
						Log.error("More than one file found for read 1: "+inputFile1+", "+file);
						return new ReturnValue(ExitStatus.INVALIDFILE);
					}
					this.inputFile1 = file;
					break;
				case 2:
					if (this.inputFile2 != null) {
						Log.error("More than one file found for read 2: "+inputFile2+", "+file);
						return new ReturnValue(ExitStatus.INVALIDFILE);
					}
					this.inputFile2 = file;
					break;
				default:
					Log.error("Cannot identify "+file+" end (read 1 or 2)");
					return new ReturnValue(ExitStatus.INVALIDFILE);
				}
			}
		} else {
			Log.error("Missing file: Two files required for workflow");
			return new ReturnValue(ExitStatus.INVALIDFILE);
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
		if (this.inputReference != null) iniFileMap.put("input_reference", this.inputReference);
		
		// Output
		iniFileMap.put("output_format", this.outputFormat.toString());
		
		iniFileMap.put("output_file_name", this.outputFileName);
		iniFileMap.put("ius_accession", this.iusAccession);
		iniFileMap.put("sequencer_run_name", this.runName);
		iniFileMap.put("barcode", this.barcode);
		iniFileMap.put("lane", this.lane);
		
		iniFileMap.put("manual_output", String.valueOf(this.manualOutput));
		
		// CutAdapt
		iniFileMap.put("do_trim", String.valueOf(this.doTrim));
		iniFileMap.put("trim_mem_mb", String.valueOf(this.trimMemoryMb));
		iniFileMap.put("r1_adapter_trim", this.read1AdapterTrim);
		iniFileMap.put("r2_adapter_trim", this.read2AdapterTrim);
		iniFileMap.put("cutadapt_r1_other_params", this.cutadapt1OtherParams);
		iniFileMap.put("cutadapt_r2_other_params", this.cutadapt2OtherParams);
		
		// BWA
		iniFileMap.put("bwa_mem_mb", String.valueOf(this.bwaMemoryMb));
		iniFileMap.put("bwa_threads", String.valueOf(this.bwaThreads));
		iniFileMap.put("bwa_other_params", this.bwaAdditionalParams);
		iniFileMap.put("bwa_pacbio_mode", String.valueOf(this.bwaPacbioMode));
		iniFileMap.put("bwa_ont_mode", String.valueOf(this.bwaOntMode));
		
		// Read Group Header (added by BWA)
		iniFileMap.put("rg_library", this.rgLibrary);
		iniFileMap.put("rg_platform", this.rgPlatform);
		iniFileMap.put("rg_platform_unit", this.rgPlatformUnit);
		iniFileMap.put("rg_sample_name", this.rgSample);
		
		// Samtools index
		iniFileMap.put("samtools_index_mem_mb", String.valueOf(this.samtoolsMemoryMb));
		
		return iniFileMap;
	}

	public static void main(String args[]){

		List<String> params = new ArrayList<String>();
		params.add("--plugin");
		params.add(BwaMemDecider.class.getCanonicalName());
		params.add("--");
		params.addAll(Arrays.asList(args));
		System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
		net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

	}
}
