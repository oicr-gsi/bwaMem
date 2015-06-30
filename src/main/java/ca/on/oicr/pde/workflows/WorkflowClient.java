package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;

import java.util.Map;

import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class WorkflowClient extends OicrWorkflow {

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
		public String getindexType() {return indexType;}
	}
	
	private String dataDir = null;
	private String reference_path = null;
	private SqwFile read1;
	private SqwFile read2;
	private String bwa, samtools, jre, picardSort;
	private String outputFilePath = null;
	private String outputIndexPath = null;
	private SqwFile outputFile;
	private SqwFile outputIndex;
	private AlignmentFormat outputFormat;
	
	private String queue; // TODO: add ini setting
	
	private void init() {
		final String binDir = getWorkflowBaseDir() + "/bin/";
		bwa = binDir + getProperty("bwa");
		samtools = binDir + getProperty("samtools");
		jre = binDir + getProperty("bundled_jre") + "/bin/java";
		picardSort = binDir + getProperty("picardsort");
		outputFormat = AlignmentFormat.valueOf(getProperty("output_format"));
	}
	

	@Override
	public void setupDirectory() {
		// this is first method called, so init from here
		init();
		
		// creates the final output
		dataDir = getPropertyOrNull("data_dir"); //TODO: default?
		if (!dataDir.endsWith("/")) dataDir += "/";
		this.addDirectory(dataDir);
	} 

	@Override
	public Map<String, SqwFile> setupFiles() {
		String input1_path = null;
		String input2_path = null;
		
		input1_path = getProperty("input_file_1");
		input2_path = getProperty("input_file_2");
		reference_path = getProperty("input_reference");

		queue = getOptionalProperty("queue", "");

		String outputFileName = getPropertyOrNull("output_file_name");
		if (outputFileName == null) {
			outputFileName = "SWID_" + getProperty("ius_accession") + "_" 
					+ getProperty("rg_library") + "_" + getProperty("sequencer_run_name") + "_" + getProperty("barcode") 
					+ "_L00" + getProperty("lane") + "_001.annotated" + outputFormat.getFileExtension();
		}
		
		outputFilePath = dataDir + outputFileName;
		
		// register input files
		read1 = this.createFile("file_in_0");
		read1.setSourcePath(input1_path);
		read1.setType("chemical/seq-na-fastq-gzip");
		read1.setIsInput(true);
		
		read2 = this.createFile("file_in_1");
		read2.setSourcePath(input2_path);
		read2.setType("chemical/seq-na-fastq-gzip");
		read2.setIsInput(true);
		
		// register output files
		outputFile = createOutputFile(outputFilePath, outputFormat.getMimeType(), Boolean.valueOf(getProperty("manual_output")));
		if (outputFormat != AlignmentFormat.SAM) {
//			outputIndexPath = outputFilePath.substring(0,outputFilePath.lastIndexOf("bam"))+"bai";
			outputIndexPath = outputFilePath + outputFormat.getIndexExtension();
			outputIndex = createOutputFile(outputIndexPath, outputFormat.getindexType(), Boolean.valueOf(getProperty("manual_output")));
		}
		
		return this.getFiles();
	}
	
	
	@Override
	public void buildWorkflow() {
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
			trimJob01.setQueue(queue);
			
			r1=trim1;
			r2=trim2;
		}
		
		// Align
		Job alignJob = getAlignJob(r1, r2);
		if (trimJob01 != null) {
			alignJob.addParent(trimJob01);
		}
		alignJob.setQueue(queue);
		
		// Index
		if (outputFormat != AlignmentFormat.SAM) {
			Job indexJob = getIndexJob();
			indexJob.addParent(alignJob);
			indexJob.setQueue(queue);
		}
	}
	
	private Job getCutAdaptJob(String read1Path, String read2Path, String read1TrimmedPath, String read2TrimmedPath) {
		Job cutAdaptJob = this.getWorkflow().createBashJob("cutAdapt");
		
		Command command = cutAdaptJob.getCommand();
		command.addArgument(jre);
		command.addArgument("-Xmx500M");
		command.addArgument("-cp "+getWorkflowBaseDir() + "/classes:"+getWorkflowBaseDir() + "/lib/"+getProperty("bundled_seqware"));
		command.addArgument("net.sourceforge.seqware.pipeline.runner.PluginRunner -p net.sourceforge.seqware.pipeline.plugins.ModuleRunner -- ");
		command.addArgument("--module ca.on.oicr.pde.utilities.workflows.modules.CutAdaptModule --no-metadata -- ");
		command.addArgument("--fastq-read-1 "+read1Path+" --fastq-read-2 "+read2Path);
		command.addArgument("--output-read-1 "+read1TrimmedPath+" --output-read-2 "+read2TrimmedPath);
		command.addArgument("--cutadapt \""+getProperty("python")+ " " +getProperty("cutadapt") +"\"");
		if (!getProperty("trim_min_quality").isEmpty()) command.addArgument("--quality "+getProperty("trim_min_quality"));
		if (!getProperty("trim_min_length").isEmpty()) command.addArgument("--minimum-length "+getProperty("trim_min_length"));
		command.addArgument("--adapters-1 "+getProperty("r1_adapter_trim")+" --adapters-2 "+getProperty("r2_adapter_trim"));
		if (!getProperty("cutadapt_r1_other_params").isEmpty()) command.addArgument("--other-parameters-1 "+getProperty("cutadapt_r1_other_params"));
		if (!getProperty("cutadapt_r2_other_params").isEmpty()) command.addArgument("--other-parameters-2 "+getProperty("cutadapt_r2_other_params"));
		cutAdaptJob.setMaxMemory(getProperty("trim_mem_mb"));
		if (!this.queue.isEmpty()) {
			cutAdaptJob.setQueue(this.queue);
		}
		return cutAdaptJob;	
	}
	
//	/**
//	 * Creates a job to run bwa mem | samtools view | SortSam (Picard), creating a .bam and .bai
//	 * 
//	 * @param read1
//	 * @param read2
//	 * 
//	 * @return the job
//	 */
//	private Job getAlignJob(String read1, String read2) {
//		Job job = this.getWorkflow().createBashJob("bwa_mem");
//		
//		Command command = job.getCommand();
//		command.addArgument("set -e; set -o pipefail;");
//		command.addArgument(bwa + " mem");
//		command.addArgument(getBwaSpecialMode());
//		String threads = getPropertyOrNull("bwa_threads");
//		if (threads != null) {
//			command.addArgument("-t " + threads);
//		}
//		String otherParams = getPropertyOrNull("bwa_other_params");
//		if (otherParams != null) {
//			command.addArgument(otherParams);
//		}
//		command.addArgument("-R " + getReadGroupHeader());
//		command.addArgument(reference_path);
//		command.addArgument(read1);
//		command.addArgument(read2);
//		command.addArgument("|");
//		
//		// pipe to samtools view (create .bam; don't write .sam to disk)
//		command.addArgument(samtools + " view");
//		command.addArgument("-bS");
//		command.addArgument("- |");
//		
//		// pipe to picard to sort and index
//		command.addArgument(jre);
//		command.addArgument("-Xmx" + getProperty("picard_mem_mb")+"M");
//		command.addArgument("-jar " +  picardSort);
//		command.addArgument("INPUT=/dev/stdin");
//		command.addArgument("OUTPUT=" + outputFilePath);
//		command.addArgument("SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true");
//		command.addArgument("TMP_DIR=" + getProperty("tmp_dir"));
//		
//		job.addFile(outputFile);
//		job.addFile(outputIndex);
//		job.setMaxMemory(getProperty("bwa_mem_mb"));
//		
//		return job;
//	}
	
	/**
	 * Creates a job to run bwa mem | samtools sort | samtools view (if necessary)
	 * 
	 * @param read1
	 * @param read2
	 * 
	 * @return the job
	 */
	private Job getAlignJob(String read1, String read2) {
		Job job = this.getWorkflow().createBashJob("bwa_mem");
		
		Command command = job.getCommand();
		command.addArgument("set -e; set -o pipefail;");
		command.addArgument(bwa + " mem");
		command.addArgument(getBwaSpecialMode());
		String threads = getPropertyOrNull("bwa_threads");
		if (threads != null) {
			command.addArgument("-t " + threads);
		}
		String otherParams = getPropertyOrNull("bwa_other_params");
		if (otherParams != null) {
			command.addArgument(otherParams);
		}
		command.addArgument("-R " + getReadGroupHeader());
		command.addArgument(reference_path);
		command.addArgument(read1);
		command.addArgument(read2);
		
		// pipe to samtools sort
		command.addArgument("|");
		command.addArgument(samtools + " sort");
		String tempDir = "tmp/"; // TODO: temp dir
		switch (outputFormat) {
		case SAM:
			command.addArgument("-O sam -T "+tempDir+" -");
			break;
		case BAM:
			command.addArgument("-O bam -T "+tempDir); // TODO: compression level
			break;
		case CRAM:
			command.addArgument("-O bam -l 0 -T "+tempDir+" -"); // uncompressed bam (samtools sort cannot output to cram)
			command.addArgument("|"); // convert to CRAM
			command.addArgument(samtools + "view");
			command.addArgument("-T " + reference_path);
			command.addArgument("-C");
			break;
		}
		command.addArgument("-o " + outputFilePath);
		command.addArgument("-");
		
		job.addFile(outputFile);
		job.setMaxMemory(getProperty("bwa_mem_mb"));
		
		return job;
	}
	
	private Job getIndexJob() {
		Job job = this.getWorkflow().createBashJob("samtools_index");
		
		Command command = job.getCommand();
		command.addArgument(samtools + " index");
		command.addArgument(outputFilePath);
		
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
