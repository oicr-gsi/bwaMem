package ca.on.oicr.pde.deciders;

import ca.on.oicr.pde.deciders.FileAttributes;
import ca.on.oicr.pde.deciders.Lims;
import ca.on.oicr.pde.deciders.OicrDecider;

import java.util.*;

import org.apache.commons.lang3.StringUtils;

import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.module.ReturnValue.ExitStatus;
import net.sourceforge.seqware.common.util.Log;

public class BwaMemDecider extends OicrDecider {
	
	// Default INI settings here are overridden by command-line parameters in init() and then maintained for all files/groups
	// Metadata details (read group header) are used in preference to these
	
	// Input
	private String inputReference = "${workflow_bundle_dir}/Workflow_Bundle_BwaMem/1.0/data/reference/bwa-0.7.9/hg19_random.fa";
	private String inputFile1, inputFile2;
	
	// Output
	private String outputPrefix = "./";
	private String outputDir = "seqware-results";
	private String outputFileName = "";
	private String iusAccession = "";
	private String library = "";
	private String runName = "";
	private String lane = "";
	private String barcode = "NoIndex"; // TODO: this never changes. Should it?
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
	private String bwaAdditionalParams = "";
	
	// Picard
	private int picardMemoryMb = 3072;
	
	// Read group header (BWA), pulled from file metadata
	private String RGLB = "";
	private String RGPL = "";
	private String RGPU = "";
	private String RGSM = "";
	
	
	public BwaMemDecider() {
		super();
		// Input
		parser.accepts("reference", "indexed reference genome").withRequiredArg();
        
		// Output
		parser.accepts("verbose", "Optional: Log all SeqWare information.");
		parser.accepts("output-prefix", "Optional: Path where the files should be copied to after analysis.").withRequiredArg();
		parser.accepts("output-dir", "Optional: Folder to put the output into relative to the output-prefix.").withRequiredArg();
		parser.accepts("output-filename", "Optional: Template type for grouping samples.").withRequiredArg();
		
		parser.accepts("ius-accession", "Optional: Used in filename IF output-filename is not specified").withRequiredArg();
		parser.accepts("library", "Optional: Used in filename IF output-filename is not specified").withRequiredArg();
		parser.accepts("sequencer-run-name", "Optional: Used in filename IF output-filename is not specified").withRequiredArg();
		parser.accepts("barcode", "Optional: Used in filename IF output-filename is not specified").withRequiredArg();
		parser.accepts("lane", "Optional: Used in filename IF output-filename is not specified").withRequiredArg();
		
		parser.accepts("manual-output", "Optional: Set output path manually.");
		
		// CutAdapt
		parser.accepts("adapter-trimming");
		parser.accepts("trim-memory", "Optional: Memory (MB) to allocate for cutadapt (default: 16384)");
		parser.accepts("read1-adapter-trim", "Adapter sequence to trim from read 1.").withRequiredArg();
		parser.accepts("read2-adapter-trim", "Adapter sequence to trim from read 2.").withRequiredArg();
		parser.accepts("read1-trim-params", "Optional: Additional cutadapt parameters for read 1.").withRequiredArg();
		parser.accepts("read2-trim-params", "Optional: Additional cutadapt parameters for read 2.").withRequiredArg();
		
		// BWA
		parser.accepts("bwa-memory", "Optional: Memory (MB) to allocate for bwa mem (default: 16384)");
		parser.accepts("bwa-threads", "Optional: Threads to use for bwa mem (default: 8).").withRequiredArg();
		parser.accepts("bwa-params", "Optional: Additional bwa mem parameters.").withRequiredArg();
		
		// Picard
		parser.accepts("picard-memory", "Optional: Memory (MB) to allocate for Picard (default: 3072)");
	}
	
	@Override
	public ReturnValue init() {
		// Parse, test, and store command-line parameters. Runs through once per decider invocation (all files)
		// Return super.init() or new ReturnValue with non-success exit code;
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
		if (this.options.has("output-prefix")) {
			outputPrefix = options.valueOf("output-prefix").toString();
			if (!outputPrefix.endsWith("/")) {
				outputPrefix += "/";
			}
		}
		if (this.options.has("output-dir")) {
			outputDir = options.valueOf("output-dir").toString();
		}
		if (this.options.has("output-filename")) {
			outputFileName = options.valueOf("output-filename").toString();
		}
		else {
			if (this.options.has("ius-accession")) {
				iusAccession = options.valueOf("ius_accession").toString();
			}
			if (this.options.has("library")) {
				library = options.valueOf("library").toString();
			}
			if (this.options.has("sequencer-run-name")) {
				runName = options.valueOf("sequencer-run-name").toString();
			}
			if (this.options.has("barcode")) {
				barcode = options.valueOf("barcode").toString();
			}
			if (this.options.has("lane")) {
				lane = options.valueOf("lane").toString();
			}
		}
		
		if (this.options.has("set-manual-path")) {
			this.manualOutput = true;
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
		if (this.options.has("bwa-params")) {
			this.bwaAdditionalParams = options.valueOf("bwa-params").toString();
		}
		
		// Picard
		if (this.options.has("picard-memory")) {
			this.picardMemoryMb = Integer.valueOf(options.valueOf("picard-memory").toString());
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
		Log.debug("CHECK FILE DETAILS:" + fm);
		if (!super.checkFileDetails(returnValue, fm)) return false;
		
		 //Get additional metadata
		this.iusAccession = returnValue.getAttribute(Header.IUS_SWA.getTitle());
		this.runName = returnValue.getAttribute(Header.SEQUENCER_RUN_NAME.getTitle());
		this.lane = returnValue.getAttribute(Header.LANE_NUM.getTitle());

		FileAttributes rv = new FileAttributes(returnValue, returnValue.getFiles().get(0));
		String groupId = StringUtils.defaultIfBlank(rv.getLimsValue(Lims.GROUP_ID), "");
		this.RGLB = rv.getLibrarySample() + groupId;
		this.RGPL = "illumina";
		this.RGSM = rv.getDonor() + groupId;
        
		this.RGPU = rv.getSequencerRun() 
				+ "-" 
				+ rv.getBarcode() 
				+ "_" 
				+ this.lane;
		
		return true;
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
						Log.error("More than one file found for read 1");
						return new ReturnValue(ExitStatus.INVALIDFILE);
					}
					this.inputFile1 = file;
					break;
				case 2:
					if (this.inputFile2 != null) {
						Log.error("More than one file found for read 2");
						return new ReturnValue(ExitStatus.INVALIDFILE);
					}
					this.inputFile2 = file;
					break;
				default:
					Log.error("Cannot identify "+file+" end (read 1 or 2).");
					return new ReturnValue(ExitStatus.INVALIDFILE);
				}
			}
		} else {
			Log.error("Missing file: Two files required");
			return new ReturnValue(ExitStatus.INVALIDFILE);
		}
		
		return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
	}
	
	@Override
	protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
		Log.debug("INI FILE:" + commaSeparatedFilePaths);
		
		Map<String, String> iniFileMap = new TreeMap<String, String>();
		
		// Inputs
		iniFileMap.put("input_file_1", this.inputFile1);
		iniFileMap.put("input_file_2", this.inputFile2);
		iniFileMap.put("input_reference", this.inputReference);
		iniFileMap.put("cutadapt_r1_other_params", this.cutadapt1OtherParams);
		iniFileMap.put("cutadapt_r2_other_params", this.cutadapt2OtherParams);
		
		// Output
		iniFileMap.put("output_prefix",this.outputPrefix);
		iniFileMap.put("output_dir", this.outputDir);
		
		iniFileMap.put("outputFileName", this.outputFileName);
		iniFileMap.put("ius_accession", this.iusAccession);
		iniFileMap.put("rg_library", this.library);
		iniFileMap.put("sequencer_run_name", this.runName);
		iniFileMap.put("barcode", this.barcode);
		iniFileMap.put("lane", this.lane);
		
		iniFileMap.put("manual_output", String.valueOf(this.manualOutput));
		
		// CutAdapt
		iniFileMap.put("do_trim", String.valueOf(this.doTrim));
		iniFileMap.put("trim_mem_mb", String.valueOf(this.trimMemoryMb));
		iniFileMap.put("r1_adapter_trim", this.read1AdapterTrim);
		iniFileMap.put("r2_adapter_trim", this.read2AdapterTrim);
		
		// BWA
		iniFileMap.put("bwa_mem_mb", String.valueOf(this.bwaMemoryMb));
		iniFileMap.put("bwa_threads", String.valueOf(this.bwaThreads));
		iniFileMap.put("bwa_other_params", this.bwaAdditionalParams);
		
		// Read Group Header (added by BWA)
		iniFileMap.put("rg_library", this.RGLB);
		iniFileMap.put("rg_platform", this.RGPL);
		iniFileMap.put("rg_platform_unit", this.RGPU);
		iniFileMap.put("rg_sample_name", this.RGSM);
		
		// Picard
		iniFileMap.put("picard_mem_mb", String.valueOf(this.picardMemoryMb));
		
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
