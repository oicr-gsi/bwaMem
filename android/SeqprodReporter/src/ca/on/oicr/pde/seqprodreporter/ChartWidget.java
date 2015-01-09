package ca.on.oicr.pde.seqprodreporter;

import android.content.Context;
import android.graphics.Canvas;
import android.graphics.Color;
import android.graphics.Paint;
import android.graphics.Path;
import android.graphics.RectF;
import android.graphics.Typeface;
import android.graphics.Paint.Style;
import android.view.View;

public class ChartWidget extends View {
	private RectF bBox;
	private RectF inBox;
    private Paint mPaint;
    private Paint tPaint;
    private Path tArc;
    private String mName = "";
    private int[] angles;
    private int[] totals;
    private float PADDING = 150.0f;
	
	public ChartWidget(Context context) {
		super(context);

		this.bBox   = new RectF();
		this.inBox  = new RectF();
		this.mPaint = new Paint(Paint.ANTI_ALIAS_FLAG);
		this.tPaint = new Paint(Paint.ANTI_ALIAS_FLAG);
		this.tArc   = new Path();
        this.angles = new int []{0, 0, 0};
		
		this.mPaint.setColor(Color.WHITE);
		this.mPaint.setStyle(Style.FILL);
		this.tPaint.setColor(Color.DKGRAY);
		this.tPaint.setStyle(Style.FILL_AND_STROKE);
		this.tPaint.setTextSize(48.0f);
		this.tPaint.setTypeface(Typeface.DEFAULT_BOLD);
		
	}
	
	protected String getmName() {
		if (null == this.mName || this.mName.isEmpty()) {
			return WorkflowStatsActivity.ALL_WORKFLOWS;
		} else {
		    return mName;
		}
	}

	protected void setmName(String mName) {
		this.mName = mName;
	}
	
	@Override
	protected void onDraw(Canvas canvas) {
		canvas.drawColor(Color.WHITE);
			
		int xDiff = 0;
		int yDiff = 0;
		int inRadius = 0;
		int centerX = this.getWidth()/2;
		int centerY = this.getHeight()/2;
		// Setup boxes, font size
		if (centerX > centerY) {
			xDiff = centerX - centerY;
			PADDING = centerY * 0.3f;
			inRadius = (int)(centerX - PADDING)/3;
		} else {
			yDiff = centerY - centerX;
			PADDING = centerX * 0.3f;
			inRadius = (int)(centerY - PADDING)/3;
		}
		
		this.tPaint.setTextSize(PADDING/3.0f);
		this.bBox.set(PADDING + xDiff, 
				      PADDING + yDiff,
				      (float) this.getWidth()  - xDiff - PADDING, 
				      (float) this.getHeight() - yDiff - PADDING);
		this.inBox.set(centerX - inRadius, 
				       centerY - inRadius, 
				       (float) centerX + inRadius, 
				       (float) centerY + inRadius);

		int currentStart = 0;
		int currentStop  = 0;
		for (int t = 0; t < ReporterActivity.types.length; t++) {
			//cycle for drawing sectors
			currentStop = currentStart + this.angles[t];
			if (currentStop == currentStart) {
				continue;
			}
			//Sector
			this.mPaint.setColor(this.getTypeFillColor(t));
			canvas.drawArc(this.bBox, currentStart, this.angles[t], true, this.mPaint);
			//Text Label
			this.tArc.reset();
		    this.tArc.addArc(this.bBox, (currentStart + currentStop)/2, 20);
		    canvas.drawTextOnPath("" + this.totals[t],tArc, 0, -this.tPaint.getTextSize(), this.tPaint);
			currentStart = currentStop;
		}
        //Inner circle (doughnut hole)
	    this.mPaint.setColor(Color.WHITE);
		canvas.drawArc(this.inBox, 0, 360, true, this.mPaint);
		
		//Paint name
		if (!this.getmName().isEmpty()) {
			this.tPaint.setTextSize(PADDING/2.5f);
			canvas.drawText(this.getmName(), PADDING/2, PADDING/1.5f, this.tPaint);
			this.tPaint.setTextSize(PADDING/3.0f);
		}
	}

	@Override
	protected void onMeasure(int widthMeasureSpec, int heightMeasureSpec) {
		setMeasuredDimension(
				widthMeasureSpec,
				heightMeasureSpec);
	}
	
	//function called from outside, repaints
	public void updateTypeData (int[] totals, String name) {
		this.setmName(name);
		this.totals = totals;
		if (null != this.totals) {
		 this.calculateAngles();	
		 this.invalidate();
		}
	}
	
	//Set angles proportionally to totals
	private void calculateAngles() {
		int sum = 0;
		float gradePerCount = 0.0f;
		
		for (int i = 0; i < this.totals.length; i++)
			sum += this.totals[i];
		if (sum == 0) {
			this.angles = new int[]{0, 0, 0};
			return;
		} 
		
		gradePerCount = (float) 360/sum;
		int angleSum = 0;
		for (int t = 0; t < this.totals.length; t++) {
			if (this.totals[t] == 0) {
			  this.angles[t] = 0;
			} else {
			  this.angles[t] = Math.round(this.totals[t] * gradePerCount);
			}
			angleSum += this.angles[t];
		}
		
		//Make sure we have sum of angles 360 degrees
		if (angleSum < 360) {
			for (int t = this.totals.length -1; t >= 0 ; t--) {
				if (this.angles[t] > 0) {
					this.angles[t] += (360 - angleSum);
					break;
				}
			}
		}
	}

    //Get colour for a certain type
	private int getTypeFillColor(int workflowType){
		int fillColor;
		switch(workflowType){
		case 0:
			fillColor = getResources().getColor(R.color.green);
			break;
		case 1:
			fillColor = getResources().getColor(R.color.red);
			break;
		case 2: 
			fillColor = getResources().getColor(R.color.orange);
			break;
		default:
			fillColor = getResources().getColor(R.color.white);
			break;
	}
		return fillColor;
	}
}
