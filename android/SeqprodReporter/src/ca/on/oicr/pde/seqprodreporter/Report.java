package ca.on.oicr.pde.seqprodreporter;

public class Report {
    private String rSampleName;
    private String rWorkflowName;
    private String rWorkflowVersion;
    private String rCreateTime;
    private String rLastmodTime;
    private String rProgress;
    //TODO we need a boolean that would flag the state of a Report (set it to true depending on rLastmodTime)
    //Later we can use this boolean flag for highlighting updated Reports
    private boolean rUpSinceLastTime;

	
	public Report(String sname, String wname, String wversion, 
			      String ctime, String ltime, boolean updated) {
		rSampleName   = sname;
		rWorkflowName = wname;
		rWorkflowVersion = wversion;
		rCreateTime = ctime;
		rLastmodTime = ltime;
		rProgress = null;
		rUpSinceLastTime = updated;
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

	public boolean getrUpdated() {
		return rUpSinceLastTime;
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
