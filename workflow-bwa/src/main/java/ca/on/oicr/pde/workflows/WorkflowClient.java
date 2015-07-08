package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class WorkflowClient extends OicrWorkflow {

    String input1_path = null;
    String input2_path = null;
    String reference_path = null;
    String dataDir = "data/";
    String outputFileName = null;
    String outputIndexName = null;
    boolean manualOutput;
    String adapter_Trimming_activated = null;
    String read1_adapterTrim = null;
    String read2_adapterTrim = null;
    String trimmedFile_1;
    String trimmedFile_2;
    //BWA parameters
    String bwa;
    String RGID;
    String RGLB;
    String RGPL;
    String RGPU;
    String RGSM;
    String additionalPicardParams;
    int readTrimming; //aln
    int numOfThreads; //aln 
    int pairingAccuracy; //aln
    int maxInsertSize; //sampe
    String readGroup;//sampe
    String bwa_aln_params;
    String bwa_sampe_params;
    SqwFile read1;
    SqwFile read2;
    SqwFile outputFile;

    String queue; 
    
    @Override
    public Map<String, SqwFile> setupFiles() {

        try {

            input1_path = getProperty("input_file_1");
            input2_path = getProperty("input_file_2");
            reference_path = getProperty("input_reference");
            //outputDir = this.getMetadata_output_dir();
            //outputPrefix = this.getMetadata_output_file_prefix();

            bwa = getProperty("bwa");

            manualOutput = Boolean.valueOf(getProperty("manual_output"));

            RGID = getProperty("rg_platform_unit");
            RGLB = getProperty("rg_library");
            RGPL = getProperty("rg_platform");
            RGPU = getProperty("rg_platform_unit");
            RGSM = getProperty("rg_sample_name");
            additionalPicardParams = getOptionalProperty("additionalPicardParams", "");

            queue = getOptionalProperty("queue", "");

            if (hasPropertyAndNotNull("outputFileName") && !getProperty("outputFileName").isEmpty()) {
                outputFileName = getProperty("outputFileName");
            } else {
		outputFileName = "SWID_" + getProperty("ius_accession") + "_" 
             	+ getProperty("rg_library") + "_" + getProperty("sequencer_run_name") + "_" + getProperty("barcode") 
             	+ "_L00" + getProperty("lane") + "_001.annotated.bam";
            }

	    outputIndexName = outputFileName.substring(0,outputFileName.lastIndexOf("bam"))+"bai";


        } catch (Exception e) {
            Logger.getLogger(WorkflowClient.class.getName()).log(Level.SEVERE, null, e);
            System.exit(1);
        }

        // registers the first input file
        read1 = this.createFile("file_in_0");
        read1.setSourcePath(input1_path);
        read1.setType("chemical/seq-na-fastq-gzip");
        read1.setIsInput(true);

        // registers the second input file
        read2 = this.createFile("file_in_1");
        read2.setSourcePath(input2_path);
        read2.setType("chemical/seq-na-fastq-gzip");
        read2.setIsInput(true);

        outputFile = createOutputFile(this.dataDir + outputFileName, "application/bam", manualOutput);
	

        return this.getFiles();
    }

    @Override
    public void setupDirectory() {
        // creates the final output 
        this.addDirectory(dataDir);
    }

    private Job doTrim(String read1Path, String read2Path, String read1TrimmedPath, String read2TrimmedPath) {
        Job cutAdaptJob = this.getWorkflow().createBashJob("cutAdapt");
        Command command = cutAdaptJob.getCommand();

        command.addArgument(getProperty("bundled_jre"));
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




    @Override
    public void buildWorkflow() {
        String r1 = read1.getProvisionedPath();
        String r2 = read2.getProvisionedPath();
        String basename1 = r1.substring(r1.lastIndexOf("/")+1,r1.lastIndexOf(".fastq.gz"));
        String basename2 = r2.substring(r2.lastIndexOf("/")+1,r2.lastIndexOf(".fastq.gz"));

	Job trimJob=null;
	if (Boolean.valueOf(getProperty("do_trim"))) {
	    String trim1=basename1+".trim.fastq.gz";
            String trim2=basename2+".trim.fastq.gz";
            trimJob = doTrim(r1,r2,trim1,trim2);
            r1=trim1;
            r2=trim2;
        }

        Job job01 = this.getWorkflow().createBashJob("bwa_align1");
        Job job02 = this.getWorkflow().createBashJob("bwa_align2");


        // Job job01 = this.getWorkflow().createBashJob("bwa_align1");
        job01.getCommand().addArgument(bwa + " aln ");
        job01.getCommand().addArgument((this.parameters("aln") == null ? " " : this.parameters("aln"))
                + reference_path + (" ")
		+ r1 + " "
                + " > " + this.dataDir + "aligned_1.sai 2> aligned_1.err");
        job01.setMaxMemory(getProperty("bwa_aln_mem_mb"));
        job01.setQueue(queue);
        if (trimJob != null) {
            job01.addParent(trimJob);
        }
        
        job02.getCommand().addArgument(bwa + " aln ");
        job02.getCommand().addArgument((this.parameters("aln") == null ? " " : this.parameters("aln"))
                + reference_path + (" ")
                + r2 + " "
                + (" > " + this.dataDir + "aligned_2.sai  2> aligned_2.err"));
        job02.setMaxMemory(getProperty("bwa_aln_mem_mb"));
        job02.setQueue(queue);
        if (trimJob != null) {
            job02.addParent(trimJob);
        }



        Job job03 = this.getWorkflow().createBashJob("bwa_sampe");
        job03.getCommand().addArgument(bwa + " sampe ");

        job03.getCommand().addArgument(this.parameters("sampe").isEmpty() ? " " : this.parameters("sampe"));
	job03.getCommand().addArgument(reference_path);
	job03.getCommand().addArgument(this.dataDir + "aligned_1.sai");
	job03.getCommand().addArgument(this.dataDir + "aligned_2.sai");
	job03.getCommand().addArgument(r1 + " " + r2);
	job03.getCommand().addArgument(" > " + this.dataDir + outputFileName + ".norg 2> " +this.dataDir + outputFileName + ".norg.log");

        job03.addParent(job01);
        job03.addParent(job02);
        job03.setQueue(queue);
        job03.setMaxMemory(getProperty("bwa_sampe_mem_mb"));
       
	Job job04 = this.getWorkflow().createBashJob("addReadGroups");
	job04.getCommand().addArgument("java -Xmx2g -jar "
                + getProperty("picard_addreadgroups") +" "
                + " RGID=" + RGID
                + " RGLB=" + RGLB
                + " RGPL=" + RGPL
                + " RGPU=" + RGPU
                + " RGSM=" + RGSM
		+ " VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate CREATE_INDEX=true"
                + " " + additionalPicardParams
                + " I=" + this.dataDir + outputFileName + ".norg"
                + " O=" + this.dataDir + outputFileName + " >> "+this.dataDir+outputFileName + ".out 2>> "+this.dataDir+outputFileName +".err");
	job04.addParent(job03);
	job04.setQueue(queue);
	job04.setMaxMemory("8000");
	job04.addFile(outputFile);
	job04.addFile(createOutputFile(this.dataDir + outputIndexName, "application/bam-index", manualOutput));

    }

    public String parameters(final String setup) {

        String paramCommand = null;
        StringBuilder a = new StringBuilder();

        try {
            if (setup.equals("aln")) {

                if (hasPropertyAndNotNull("readTrimming") && !getProperty("readTrimming").isEmpty()) {
                    readTrimming = Integer.parseInt(getProperty("readTrimming"));
                    a.append(" -q ");
                    a.append(readTrimming);
                    a.append(" ");
                }

                if (hasPropertyAndNotNull("numOfThreads") && !getProperty("numOfThreads").isEmpty()) {
                    numOfThreads = Integer.parseInt(getProperty("numOfThreads"));
                    a.append(" -t ");
                    a.append(numOfThreads);
                    a.append(" ");
                }

                if (hasPropertyAndNotNull("pairingAccuracy")&& !getProperty("pairingAccuracy").isEmpty()) {
                    pairingAccuracy = Integer.parseInt(getProperty("pairingAccuracy"));
                    a.append(" -R ");
                    a.append(pairingAccuracy);
                    a.append(" ");
                }
                if (hasPropertyAndNotNull("bwa_aln_params") && !getProperty("bwa_aln_params").isEmpty()) {
                    bwa_aln_params = getProperty("bwa_aln_params");
                    a.append(" ");
                    a.append(bwa_aln_params);
                    a.append(" ");
                }
                paramCommand = a.toString();
                return paramCommand;
            }

            if (setup.equals("sampe")) {

                if (hasPropertyAndNotNull("maxInsertSize") && !getProperty("maxInsertSize").isEmpty()) {
                    maxInsertSize = Integer.parseInt(getProperty("maxInsertSize"));
                    a.append(" -a ");
                    a.append(maxInsertSize);
                    a.append(" ");
                }

                if (hasPropertyAndNotNull("readGroup") && !getProperty("readGroup").isEmpty()) {
                    a.append(" -r ");
                    a.append(readGroup);
                    a.append(" ");
                }

                if (hasPropertyAndNotNull("bwa_sampe_params") && !getProperty("bwa_sampe_params").isEmpty()) {
                    bwa_sampe_params = getProperty("bwa_sampe_params");
                    a.append(" ");
                    a.append(bwa_sampe_params);
                    a.append(" ");
                }
                paramCommand = a.toString();
                return paramCommand;
            }

//            if (setup.equals("view")) {
//
//                if (!getProperty("samtools_view_params").isEmpty()) {
//                    samtools_view_params = getProperty("samtools_view_params");
//                    a.append(" ");
//                    a.append(samtools_view_params);
//                    a.append(" ");
//                }
//                paramCommand = a.toString();
//                return paramCommand;
//            }
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
        return paramCommand;
    }
}
