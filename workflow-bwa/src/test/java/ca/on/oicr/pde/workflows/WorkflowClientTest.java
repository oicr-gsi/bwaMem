package ca.on.oicr.pde.workflows;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import org.apache.commons.io.FileUtils;
import org.testng.annotations.Test;

/**
 *
 * @author mlaszloffy
 */
public class WorkflowClientTest {

    public WorkflowClientTest() {

    }

    @Test
    public void testDefaultConfig() throws IOException {
        WorkflowClient w = new WorkflowClient();
        w.setConfigs(getDefaultConfig());
        buildWorkflowModel(w);
    }

    @Test
    public void testEmptyParams() throws IOException {
        Map<String, String> config = getDefaultConfig();
        config.put("input_file_1", "/tmp/abc.fastq.gz");
        config.put("input_file_2", "/tmp/abcd.fastq.gz");
        config.put("outputFileName", "1.bam");
        config.put("additionalPicardParams", "");
        config.put("maxInsertSize", "");
        config.put("readGroup", "");
        config.put("bwa_sampe_params", "");

        WorkflowClient w = new WorkflowClient();
        w.setConfigs(config);
        buildWorkflowModel(w);
    }

    private Map<String, String> getDefaultConfig() throws IOException {
        Properties p = new Properties();
        File defaults = new File("workflow/config/workflow.ini");
        p.load(FileUtils.openInputStream(defaults));
        return new HashMap<>((Map) p);
    }

    private void buildWorkflowModel(AbstractWorkflowDataModel workflowObject) {
        AbstractWorkflowDataModel.prepare(workflowObject);
        workflowObject.setupDirectory();
        workflowObject.setupFiles();
        workflowObject.setupWorkflow();
        workflowObject.setupEnvironment();
        workflowObject.buildWorkflow();
        workflowObject.wrapup();
    }

}
