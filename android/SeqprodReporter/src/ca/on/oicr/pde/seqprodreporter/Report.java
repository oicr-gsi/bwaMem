package ca.on.oicr.pde.seqprodreporter;

public class Report {
    private String rSampleName;
    private String rWorkflowName;
    private String rWorkflowVersion;
    private String rCreateTime;
    private String rLastmodTime;
    private String rProgress;

	
	public Report(String sname, String wname, String wversion, String ctime, String ltime) {
		rSampleName   = sname;
		rWorkflowName = wname;
		rWorkflowVersion = wversion;
		rCreateTime = ctime;
		rLastmodTime = ltime;
		rProgress = null;
	}


	public String getrSampleName() {
		return rSampleName;
	}


	public String getrLastmodTime() {
		return rLastmodTime;
	}


	public String getrWorkflowName() {
		return rWorkflowName;
	}


	public String getrWorkflowVersion() {
		return rWorkflowVersion;
	}


	public String getrCreateTime() {
		return rCreateTime;
	}


	public void setrProgress(String progress) {
		rProgress = progress;
	}
	
	public String getrProgress() {
		return rProgress;
	}

	
	public int progressValue() {
		int progress = 0;
		if (null == rProgress)
			return progress;
		
    	try {
			progress = Integer.parseInt(rProgress);
		} catch (NumberFormatException ne) {
			ne.printStackTrace();
		}
    	return progress;
	}

}
