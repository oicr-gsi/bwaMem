package com.github.seqware;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class WorkflowClient extends AbstractWorkflowDataModel {

        String input1_path = null;
        String input2_path = null;
        String reference_path = null;
        String outputPrefix = null;
        String outputDir = null;
        String finalOutputDir = null;
        //String filename; 
    @Override
    public Map<String, SqwFile> setupFiles() {
        
        try {
    	
        input1_path = getProperty("input_file_1");
        input2_path = getProperty("input_file_2");
        reference_path = getProperty("input_reference");
        outputPrefix = getProperty("output_prefix");
        outputDir = getProperty("output_dir");
        
        finalOutputDir = outputPrefix + outputDir + "/" ;
        }
    	
        catch (Exception e)
        {
               Logger.getLogger(WorkflowClient.class.getName()).log(Level.SEVERE, null, e);
               return(null);
        }
      
    	  
        // registers the first input file
        SqwFile file0 = this.createFile("file_in_1");
        file0.setSourcePath(input1_path);
        file0.setType("chemical/seq-na-fastq-gzip");
        file0.setIsInput(true);
        // registers the second input file 
        SqwFile file1 = this.createFile("file_in_2");
        file1.setSourcePath(input2_path);
        file1.setType("chemical/seq-na-fastq-gzip");
        file1.setIsInput(true);

        // registers an output file
        SqwFile file2 = this.createFile("file_out");
        file2.setSourcePath("finalOut.bam");
        file2.setType("application/bam");
        file2.setIsOutput(true);
        file2.setForceCopy(true);
        file2.setOutputPath(finalOutputDir);
       
      return this.getFiles();
    }
    
    @Override
    public void setupDirectory() {
        // creates a dir1 directory in the current working directory where the workflow runs
        this.addDirectory(finalOutputDir);
    }
    
    @Override
    public void buildWorkflow() {

        // a simple bash job to call mkdir

	Job job01 = this.getWorkflow().createBashJob("bwa_align1");
	job01.getCommand().addArgument("bwa aln "+reference_path+(" ") 
	+this.getFiles().get("file_in_1").getProvisionedPath()+("> aligned_1.sai"));
        job01.setMaxMemory("16000");

        Job job02 = this.getWorkflow().createBashJob("bwa_align2");
        job01.getCommand().addArgument("bwa aln "+ reference_path+(" ") 
        +this.getFiles().get("file_in_2").getProvisionedPath()+("> aligned_2.sai"));
        job02.setMaxMemory("16000");
        
        Job job03 = this.getWorkflow().createBashJob("bwa_sampe");
        job03.setCommand("bwa sampe " + reference_path 
           +(" aligned_1.sai")
           +(" aligned_2.sai ")
           + input1_path +(" ")
           + input2_path
           + (" > file_out.sam")); 
        job03.addParent(job01);
        job03.addParent(job02);
        job03.setMaxMemory("16000");
        //   bwa sampe /oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.6.2/hg19_random.fa test2.sai test3.sai /.mounts/labs/PDE/data/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.cleaned.fixed.bam.R1.fastq.gz /.mounts/labs/PDE/data/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.cleaned.fixed.bam.R2.fastq.gz > aln_ne.sam

        Job job04 = this.getWorkflow().createBashJob("sam_job");
        job04.getCommand().addArgument("samtools view -bS "
                + "file_out.sam > finalOut.bam");
        job04.addParent(job03);
        job04.setMaxMemory("16000");
        
    }

}
