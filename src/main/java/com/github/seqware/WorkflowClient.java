


import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class WorkflowClient extends AbstractWorkflowDataModel {

    @Override
    public Map<String, SqwFile> setupFiles() {

      try {

        // register an input file
        SqwFile file0 = this.createFile("file_in_1");
        file0.setSourcePath(getProperty("input_file1"));
        file0.setType("chemical/seq-na-fastq-gzip");
        file0.setIsInput(true);

	SqwFile file1 = this.createFile("file_in_2");
	file1.setSourcePath("input_file2");
	file1.setType("chemical/seq-na-fastq-gzip");
	file1.setIsInput(true);

        // register an output file
        SqwFile file2 = this.createFile("file_out");
        file2.setSourcePath("dir1/output");
        file2.setType("text/plain");
        file2.setIsOutput(true);
        file2.setForceCopy(true);
        // if output_file is set in the ini then use it to set the destination of this file
        if (hasPropertyAndNotNull("output_file")) { file2.setOutputPath(getProperty("output_file")); }
        return this.getFiles();

      } catch (Exception ex) {
        Logger.getLogger(WorkflowClient.class.getName()).log(Level.SEVERE, null, ex);
        return(null);
      }
    }
    
    @Override
    public void setupDirectory() {
        // creates a dir1 directory in the current working directory where the workflow runs
        this.addDirectory("dir1");
    }
    
    @Override
    public void buildWorkflow() {

        // a simple bash job to call mkdir
        Job job00 = this.getWorkflow().createBashJob("bash_mkdir");
        job00.getCommand().addArgument("mkdir sai-files");

	Job job01 = this.getWorkflow().createBashJob("bash_bwa aln");
	job01.getCommand().addArgument("bwa aln "+getProperty("input_reference")+(" ") 
	+this.getFiles().get("file_in_1").getProvisionedPath()+("> aligned_1.sai"));

	
        Job job02 = this.getWorkflow().createBashJob("bash_bwa aln");
        job01.getCommand().addArgument("bwa aln "+getProperty("input_reference")+(" ") 
        +this.getFiles().get("file_in_2").getProvisionedPath()+("> aligned_2.sai"));

        
        // a simple bash job to copy an input to an output file
        Job job03 = this.getWorkflow().createBashJob("bash_bwa sampe");
        job11.setCommand("bwa sampe " + this.getFiles().get("file_in_0").getProvisionedPath() + " dir1/output");
        job11.addParent(job00);
               

    }

}
