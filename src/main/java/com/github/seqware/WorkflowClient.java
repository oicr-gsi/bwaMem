package com.github.seqware;
import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class WorkflowClient extends OicrWorkflow {

        String input1_path = null;
        String input2_path = null;
        String reference_path = null;
        String outputPrefix = null;
        String outputDir = null;
        String finalOutputDir = null;
        String outputFileName = null;
        
        
        //BWA parameters
        int readTrimming; //aln
        int numOfThreads; //aln 
        int pairingAccuracy; //aln
        int maxInsertSize; //sampe
        String readGroup;//sampe
        String bwa_aln_params;
        String bwa_sampe_params;
        String samtools_view_params;
        
        SqwFile file0;
        SqwFile file1;
        SqwFile file2;
        
    @Override
    public Map<String, SqwFile> setupFiles() {
        
        try {
    	
        input1_path = getProperty("input_file_1");
        input2_path = getProperty("input_file_2");
        reference_path = getProperty("input_reference");
        outputPrefix = getProperty("output_prefix");
        outputDir = getProperty("output_dir");        
        finalOutputDir = outputPrefix + outputDir + "/" ;
        
        if (getProperty("ouputFileName") != null){
            outputFileName = getProperty("outputFileName");
        }
        else {
            outputFileName = ((input1_path.substring(input1_path.lastIndexOf("/") + 1))
                            +(input2_path.substring(input2_path.lastIndexOf("/") +1))
                            +(".bam"));
        }
        
        }
    	
        catch (Exception e)
        {
               Logger.getLogger(WorkflowClient.class.getName()).log(Level.SEVERE, null, e);
               return(null);
        }

        // registers the first input file
        file0 = this.createFile("file_in_1");
        file0.setSourcePath(input1_path);
        String[] filepath = input1_path.split(".");
            if (filepath.length >=2){
                //for (int i = filepath.length; i > filepath.length -1; i--){
                    if (filepath[filepath.length].equals("gz") && filepath[filepath.length-1].equals("fastq")) {
                        file0.setType("chemical/seq-na-fastq-gzip");
                    } else if (filepath[filepath.length].equals("fastq")) {
                        file0.setType("chemical/seq-na-fastq");
                    } else {
                    System.exit(1);
     
                    }                 
                }
            
        file0.setIsInput(true);

        // registers the second input file 
        file1 = this.createFile("file_in_2");
        file1.setSourcePath(input2_path);
        String[] filepath2 = input2_path.split(".");
            if (filepath2.length >=2){
                    if (filepath[filepath.length].equals("gz") && filepath[filepath.length-1].equals("fastq")) {
                        file0.setType("chemical/seq-na-fastq-gzip");
                    } else if (filepath[filepath.length].equals("fastq")) {
                        file0.setType("chemical/seq-na-fastq");
                    } else {
                    System.exit(1);
     
                    }              
            }
        file1.setIsInput(true);

        //file2 = createOutputFile(outputFileName, "application/bam", setManualpath);
        
        // registers an output file
        file2 = this.createFile("file_out");
        file2.setSourcePath(outputFileName);
        file2.setType("application/bam");
        file2.setIsOutput(true);
        file2.setForceCopy(true);
        file2.setOutputPath(finalOutputDir);
       
      return this.getFiles();
    }
    
    @Override
    public void setupDirectory() {
        // creates the final output 
        this.addDirectory(finalOutputDir);
    }
    
    @Override
    public void buildWorkflow() {

	Job job01 = this.getWorkflow().createBashJob("bwa_align1");
	job01.getCommand().addArgument(this.getWorkflowBaseDir()+"/bin/bwa-0.6.2/bwa aln "
            +(this.parameters("aln")==null ? " " :this.parameters("aln"))  
            +reference_path+(" ") 
            +this.getFiles().get("file_in_1").getProvisionedPath()+(" > aligned_1.sai"));
        job01.setMaxMemory("16000");

        Job job02 = this.getWorkflow().createBashJob("bwa_align2");
        job02.getCommand().addArgument(this.getWorkflowBaseDir()+"/bin/bwa-0.6.2/bwa aln "
            +(this.parameters("aln")==null ? " " :this.parameters("aln"))
            + reference_path+(" ")
            +this.getFiles().get("file_in_2").getProvisionedPath()
            +(" > aligned_2.sai"));
        job02.setMaxMemory("16000");
        
        Job job03 = this.getWorkflow().createBashJob("bwa_sampe");
        job03.getCommand().addArgument(this.getWorkflowBaseDir()+"/bin/bwa-0.6.2/bwa sampe " 
           +(this.parameters("sampe")==null ? " " :this.parameters("sampe")) 
           +reference_path 
           +(" aligned_1.sai")
           +(" aligned_2.sai ")
           + input1_path +(" ")
           + input2_path
           + (" > file_out.sam")); 
        job03.addParent(job01);
        job03.addParent(job02);
        job03.setMaxMemory("16000");
        
        Job job04 = this.getWorkflow().createBashJob("samToBam_job");
        
        job04.getCommand().addArgument(this.getWorkflowBaseDir()+"/bin/samtools-0.1.19/samtools view -bS "
                +(this.parameters("view")==null ? " " :this.parameters("view")) 
                + "file_out.sam > "
                +outputFileName);
        job04.addParent(job03);
        job04.addFile(file2);
        job04.setMaxMemory("16000");
    }   
        
     public String parameters(String setup) {
         
         String paramCommand = null;
         StringBuilder a = new StringBuilder();

         try {  
            if (setup.equals("aln")){
             
                if (getProperty("readTrimming")!= null) {
                   readTrimming = Integer.parseInt(getProperty("readTrimming"));
                   a.append(" -q ");
                   a.append(readTrimming);
                   }

               if (getProperty("numOfThreads")!= null) {
                   numOfThreads = Integer.parseInt(getProperty("numOfThreads"));
                   a.append(" -t ");
                   a.append(numOfThreads);
                   }

               if (getProperty("pairingAccuracy") != null){
                   pairingAccuracy = Integer.parseInt(getProperty("pairingAccuracy"));
                   a.append(" -R ");
                   a.append(pairingAccuracy);
                   }
               if  (getProperty("bwa_aln_params") != null){
                   bwa_aln_params = getProperty("bwa_aln_params");
                   a.append(" ");
                   a.append(bwa_aln_params);
                   } 

                paramCommand = a.toString();
                return paramCommand ;  
                }
            
            if (setup.equals("sampe")){
            
               if (getProperty("maxInsertSize") != null){
                   maxInsertSize = Integer.parseInt(getProperty("maxInsertSize"));
                   a.append(" -a ");
                   a.append(maxInsertSize);
                   }
               
               if (getProperty("readGroup") != null){
                   a.append(" -r ");
                   a.append(readGroup); 
               }
 
               if  (getProperty("bwa_sampe_params") != null){
                   bwa_sampe_params = getProperty("bwa_sampe_params");
                   a.append(" ");
                   a.append(bwa_sampe_params);
                   } 
                paramCommand = a.toString();
                return paramCommand ;
            }
            
            if (setup.equals("view")){
                
                if  (getProperty("samtools_view_params") != null){
                   samtools_view_params = getProperty("samtools_view_params");
                   a.append(" ");
                   a.append(samtools_view_params);
                   } 
                paramCommand = a.toString();
                return paramCommand ;
            }
         }
         
         catch (Exception e){
            }
            return paramCommand ;
        }
 
    }


