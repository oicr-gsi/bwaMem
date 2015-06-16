package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class WorkflowClient extends OicrWorkflow {

	private String dataDir = null;
	private String reference_path = null;
	private SqwFile read1;
	private SqwFile read2;
	private String bwa, samtools, jre, picardSort;
	private String outputBamPath = null;
	private String outputIndexPath = null;
	private SqwFile outputBam;
	private SqwFile outputBai;
	
	String queue;
	
	private void init() {
		final String binDir = getWorkflowBaseDir() + "/bin/";
//		cutadapt = binDir + getProperty("cutadapt");
		bwa = binDir + getProperty("bwa");
		samtools = binDir + getProperty("samtools");
		jre = binDir + getProperty("bundled_jre") + "/bin/java";
		picardSort = binDir + getProperty("picardsort");
	}
	

	@Override
	public void setupDirectory() {
		// this is first method called, so init from here
		init();
		
		// creates the final output
		dataDir = getPropertyOrNull("data_dir");
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

		String outputBamName = getPropertyOrNull("outputBamPath");
		if (outputBamName == null) {
			outputBamName = "SWID_" + getProperty("ius_accession") + "_" 
					+ getProperty("rg_library") + "_" + getProperty("sequencer_run_name") + "_" + getProperty("barcode") 
					+ "_L00" + getProperty("lane") + "_001.annotated.bam";
		}
		
		outputBamPath = dataDir + outputBamName;
		outputIndexPath = outputBamPath.substring(0,outputBamPath.lastIndexOf("bam"))+"bai";
		
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
		outputBam = createOutputFile(outputBamPath, "application/bam", Boolean.valueOf(getProperty("manual_output")));
		outputBai = createOutputFile(outputIndexPath, "application/bam-index", Boolean.valueOf(getProperty("manual_output")));
		
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
//		Job trimJob02=null;
		if (Boolean.valueOf(getProperty("do_trim"))) {
			String trim1=this.dataDir+basename1+".trim.fastq.gz";
			String trim2=this.dataDir+basename2+".trim.fastq.gz";
			
			trimJob01 = getCutAdaptJob(r1, r2, trim1, trim2);
//			trimJob01 = getCutAdaptJob("cutAdapt1", r1, trim1, getProperty("r1_adapter_trim"), getPropertyOrNull("trim_r1_other_params"));
			trimJob01.setQueue(queue);
//			trimJob02 = getCutAdaptJob("cutAdapt2", r2, trim2, getProperty("r2_adapter_trim"), getPropertyOrNull("trim_r2_other_params"));
//			trimJob02.setQueue(queue);
			
			r1=trim1;
			r2=trim2;
		}
		
		
		// Create bam (bwa mem | samtools view -bS)
		Job bamJob = getBamJob(r1, r2);
		if (trimJob01 != null) {
			bamJob.addParent(trimJob01);
//			bamJob.addParent(trimJob02);
		}
		bamJob.setQueue(queue);
		
		
		// Create bai (samtools index)
//		Job indexJob = getIndexJob();
//		indexJob.addParent(bamJob);
//		indexJob.setQueue(queue);
	}
	
//	/**
//	 * Creates a job to run CutAdapt on a read
//	 * 
//	 * @param jobName
//	 * @param readPath input file
//	 * @param readTrimmedPath output file
//	 * @param adapterSequence
//	 * @param additionalParams a String to add onto the command line formatted like "-arg val -arg val..." May be null
//	 * 
//	 * @return the job
//	 */
//	private Job getCutAdaptJob(String jobName, String readPath, String readTrimmedPath, String adapterSequence, String additionalParams) {
//		Job job = this.getWorkflow().createBashJob(jobName);
//		
//		Command command = job.getCommand();
//		command.addArgument(cutadapt);
//		command.addArgument("-a " + adapterSequence);
//		
//		String minLength = getPropertyOrNull("trim_min_length");
//		if (minLength != null) {
//			command.addArgument("-m " + minLength); // minimum read length (after trimming) to keep
//		}
//		
//		String minQuality = getPropertyOrNull("trim_min_quality") ;
//		if (minQuality != null) {
//			command.addArgument("-q " + minQuality); // trim ends below this quality before adapter removal
//		}
//		
//		if (additionalParams != null) {
//			command.addArgument(additionalParams);
//		}
//		
//		command.addArgument("-o " + readTrimmedPath); // output
//		command.addArgument(readPath); // input
//		
//		job.setMaxMemory(getProperty("trim_mem_mb"));
//		
//		return job;
//	}
	
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
	
	/**
	 * Creates a job to run bwa mem | samtools view | SortSam (Picard), creating a .bam and .bai
	 * 
	 * @param read1
	 * @param read2
	 * 
	 * @return the job
	 */
	private Job getBamJob(String read1, String read2) {
		Job job = this.getWorkflow().createBashJob("bwa_mem");
		
		Command command = job.getCommand();
		command.addArgument("set -e; set -o pipefail;");
		command.addArgument(bwa + " mem");
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
		command.addArgument("|");
		
		// pipe to samtools view (create .bam; don't write .sam to disk)
		command.addArgument(samtools + " view");
		command.addArgument("-bS");
		command.addArgument("- |");
		
		// pipe to picard to sort and index
		command.addArgument(jre);
		command.addArgument("-Xmx" + getProperty("picard_memory")+"M");
		command.addArgument("-jar " +  picardSort);
		command.addArgument("INPUT=/dev/stdin");
		command.addArgument("OUTPUT=" + outputBamPath);
		command.addArgument("SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true");
		command.addArgument("TMP_DIR=" + getProperty("tmp_dir"));
		
		job.addFile(outputBam);
		job.addFile(outputBai);
		job.setMaxMemory(getProperty("bwa_mem_mb"));
		
		return job;
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
	 * Creates a job to create the BAM index (.bai) via samtools
	 * 
	 * @return the job
	 */
//	private Job getIndexJob() {
//		Job job = this.getWorkflow().createBashJob("samtools_index");
//		
//		Command command = job.getCommand();
//		command.addArgument(samtools + " index");
//		command.addArgument(outputBamPath);
//		
//		job.addFile(outputBai);
//		job.setMaxMemory(getProperty("index_mem_mb"));
//		
//		return job;
//	}
	
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
