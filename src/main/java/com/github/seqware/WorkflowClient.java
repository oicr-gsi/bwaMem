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
        String fileName = null;
        String o_file;
        
    @Override
    public Map<String, SqwFile> setupFiles() {
        
        try {
    	
        input1_path = getProperty("input_file_1");
        input2_path = getProperty("input_file_2");
        reference_path = getProperty("input_reference");
        outputPrefix = getProperty("output_prefix");
        outputDir = getProperty("output_dir");
        
        outputDir = outputDir.lastIndexOf("/") == (outputDir.length() - 1) ? 
                outputDir : outputDir + "/";

        finalOutputDir = outputPrefix + outputDir ;
       
        }
    	
        catch (Exception e)
        {
                //e.printStackTrace();
        }
      try {
    	  
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
        file2.setSourcePath(o_file);
        file2.setType("text/sam");
        file2.setIsOutput(true);
        file2.setForceCopy(true);
        file2.setOutputPath(finalOutputDir);
       
        
        // if output_file is set in the ini then use it to set the destination of this file
        //if (hasPropertyAndNotNull("output_file")) { file2.setOutputPath(getProperty("output_file")); }
        //return this.getFiles();

      } catch (Exception ex) {
        Logger.getLogger(WorkflowClient.class.getName()).log(Level.SEVERE, null, ex);
        return(null);
      }
      
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

	Job job01 = this.getWorkflow().createBashJob("bash_bwa aln");
	job01.getCommand().addArgument("bwa aln "+reference_path+(" ") 
	+this.getFiles().get("file_in_1").getProvisionedPath()+("> aligned_1.sai"));

        Job job02 = this.getWorkflow().createBashJob("bash_bwa aln");
        job01.getCommand().addArgument("bwa aln "+ reference_path+(" ") 
        +this.getFiles().get("file_in_2").getProvisionedPath()+("> aligned_2.sai"));

        
        Job job03 = this.getWorkflow().createBashJob("bash_bwa sampe");
        job03.setCommand("bwa sampe " + reference_path 
           +(" aligned_1.sai")
           +(" aligned_2.sai")
           + input1_path
           + input2_path
           + (" > o_file.sam"));
        
        job03.addParent(job01);
        job03.addParent(job02);
        
    }

}
