package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;

import java.util.Map;

import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

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
		};
		
		public String getFileExtension() {return fileExtension;}
		public String getMimeType() {return mimeType;}
		public String getIndexExtension() {return indexExtension;}
		public String getIndexType() {return indexType;}
	}
	
	private String dataDir = null;
	private String tempDir = null;
	private String reference_path = null;
	private SqwFile read1;
	private SqwFile read2;
	private String bwa, samtools, jre;
	private String outputFilePath = null;
	private String outputIndexPath = null;
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
	}
	

	@Override
	public void setupDirectory() {
		// this is first method called, so init from here
		init();
		
		// create the final output directory
		dataDir = getProperty("data_dir");
		if (!dataDir.endsWith("/")) dataDir += "/";
		this.addDirectory(dataDir);
		
		tempDir = getProperty("tmp_dir");
		if (!tempDir.endsWith("/")) tempDir += "/";
		this.addDirectory(tempDir);
	} 

	@Override
	public Map<String, SqwFile> setupFiles() {
		String input1_path = null;
		String input2_path = null;
		
		input1_path = getProperty("input_file_1");
		input2_path = getProperty("input_file_2");
		reference_path = getProperty("input_reference");
		
		String outputFileName = getOutputFileName();
		
		outputFilePath = dataDir + outputFileName;
		
		// register input files
		read1 = this.createFile("file_in_0");
		read1.setSourcePath(input1_path);
		read1.setType(input1_path.endsWith(".gz") ? FASTQ_GZIP_MIMETYPE : FASTQ_MIMETYPE);
		read1.setIsInput(true);
		
		read2 = this.createFile("file_in_1");
		read2.setSourcePath(input2_path);
		read2.setType(input1_path.endsWith(".gz") ? FASTQ_GZIP_MIMETYPE : FASTQ_MIMETYPE);
		read2.setIsInput(true);
		
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
			
			if (groupId != null) sb.append(groupId).append("_");
			
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
		
		String r1 = read1.getProvisionedPath();
		String r2 = read2.getProvisionedPath();
		String basename1 = r1.substring(r1.lastIndexOf("/")+1,r1.lastIndexOf(".fastq.gz"));
		String basename2 = r2.substring(r2.lastIndexOf("/")+1,r2.lastIndexOf(".fastq.gz"));
		
		
		// cutadapt (optional)
		Job trimJob01=null;
		if (Boolean.valueOf(getProperty("do_trim"))) {
			String trim1=this.dataDir+basename1+".trim.fastq.gz";
			String trim2=this.dataDir+basename2+".trim.fastq.gz";
			
			trimJob01 = getCutAdaptJob(r1, r2, trim1, trim2);
			if (queue != null) trimJob01.setQueue(queue);
			
			r1=trim1;
			r2=trim2;
		}
		
		// Align, sort, and convert
		Job alignJob = getAlignJob(r1, r2);
		if (trimJob01 != null) {
			alignJob.addParent(trimJob01);
		}
		if (queue != null) alignJob.setQueue(queue);
		
		// Index
		if (outputFormat != AlignmentFormat.SAM) {
			Job indexJob = getIndexJob();
			indexJob.addParent(alignJob);
			if (queue != null) indexJob.setQueue(queue);
		}
	}
	
	/**
	 * Creates a job to trim adapters from the reads using CutAdapt
	 * 
	 * @param read1Path input file
	 * @param read2Path input file
	 * @param read1TrimmedPath output filename for trimmed read
	 * @param read2TrimmedPath output filename for trimmed read
	 * @return
	 */
	private Job getCutAdaptJob(String read1Path, String read2Path, String read1TrimmedPath, String read2TrimmedPath) {
		Job job = this.getWorkflow().createBashJob("cutAdapt");
		job.setMaxMemory(getProperty("trim_mem_mb"));
		
		Command cmd = job.getCommand();
		cmd.addArgument(jre);
		cmd.addArgument("-Xmx500M");
		cmd.addArgument("-cp "+getWorkflowBaseDir() + "/classes:"+getWorkflowBaseDir() + "/lib/"+getProperty("bundled_seqware"));
		cmd.addArgument("net.sourceforge.seqware.pipeline.runner.PluginRunner -p net.sourceforge.seqware.pipeline.plugins.ModuleRunner -- ");
		cmd.addArgument("--module ca.on.oicr.pde.utilities.workflows.modules.CutAdaptModule --no-metadata -- ");
		cmd.addArgument("--fastq-read-1 "+read1Path+" --fastq-read-2 "+read2Path);
		cmd.addArgument("--output-read-1 "+read1TrimmedPath+" --output-read-2 "+read2TrimmedPath);
		cmd.addArgument("--cutadapt \""+getProperty("python")+ " " +getProperty("cutadapt") +"\"");
		if (!getProperty("trim_min_quality").isEmpty()) cmd.addArgument("--quality "+getProperty("trim_min_quality"));
		if (!getProperty("trim_min_length").isEmpty()) cmd.addArgument("--minimum-length "+getProperty("trim_min_length"));
		cmd.addArgument("--adapters-1 "+getProperty("r1_adapter_trim")+" --adapters-2 "+getProperty("r2_adapter_trim"));
		if (!getProperty("cutadapt_r1_other_params").isEmpty()) cmd.addArgument("--other-parameters-1 "+getProperty("cutadapt_r1_other_params"));
		if (!getProperty("cutadapt_r2_other_params").isEmpty()) cmd.addArgument("--other-parameters-2 "+getProperty("cutadapt_r2_other_params"));
		
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
		cmd.addArgument(read2);
		
		// pipe to samtools sort
		cmd.addArgument("|");
		cmd.addArgument(samtools + " sort");
		
		switch (outputFormat) {
		case SAM:
			cmd.addArgument("-O sam -T "+this.tempDir);
			break;
		case BAM:
			cmd.addArgument("-O bam -T "+this.tempDir);
			break;
		case CRAM:
			cmd.addArgument("-O bam -l 0 -T "+this.tempDir+" -"); // uncompressed bam (samtools sort cannot output to cram)
			
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
		if (!hasPropertyAndNotNull(key)) return null;
		String val = getProperty(key);
		return val.isEmpty() ? null : val;
	}
	
}
