package ca.on.oicr.pde.seqprodreporter;

import android.annotation.SuppressLint;
import android.app.Activity;
import android.app.AlertDialog;
import android.app.Dialog;
import android.app.DialogFragment;
import android.content.DialogInterface;
import android.content.SharedPreferences;
import android.os.Bundle;
import android.text.Editable;
import android.text.TextWatcher;
import android.util.Log;
import android.view.LayoutInflater;
import android.view.View;
import android.view.WindowManager;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.EditText;
import android.widget.Spinner;

/**
 * This is a dialog that takes care of setting parameters such as 
 * Update Frequency, Time Range, Notification type and host URL 
 */
public class PrefDialogFragment extends DialogFragment implements AdapterView.OnItemSelectedListener {

	private String updateFrequency;
	private String timeRange;
	private String notificationType;
	private String hostUrl;
	
	public String getUpdateFrequency() {
		return updateFrequency;
	}

	public String getHostUrl() {
		return hostUrl;
	}
	
	public String getTimeRange() {
		return timeRange;
	}
	
	public String getNotificationType() {
		return notificationType;
	}
	
	// Use this instance of the interface to deliver action events
    PrefDialogListener mListener;
    
	public interface PrefDialogListener {
		public void onDialogPositiveClick(PrefDialogFragment dialog);
	}
	
	
	@SuppressLint("InflateParams")
	@Override
    public Dialog onCreateDialog(Bundle savedInstanceState) {
        // Build the dialog and set up the button click handlers
		AlertDialog.Builder builder = new AlertDialog.Builder(getActivity());
		LayoutInflater inflater = getActivity().getLayoutInflater();
		View rootView = inflater.inflate(R.layout.configure_dialog, null);  	
    	SharedPreferences sp = getActivity().getSharedPreferences(ReporterActivity.PREFERENCE_FILE, Activity.MODE_PRIVATE);
    	
	    builder.setView(rootView)
               .setPositiveButton(R.string.ok, new DialogInterface.OnClickListener() {
                   public void onClick(DialogInterface dialog, int id) {
                       mListener.onDialogPositiveClick(PrefDialogFragment.this);
                   }
               });
	    
	    Spinner upSpinner = (Spinner) rootView.findViewById(R.id.automatic_updates_settings);
	    ArrayAdapter<CharSequence> upAdapter = ArrayAdapter.createFromResource(this.getActivity(),
	            R.array.pref_automaticUpdates_entries, android.R.layout.simple_spinner_item);
	    upAdapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
	    upSpinner.setAdapter(upAdapter);
	    upSpinner.setSelection(this.findItemPosition(upSpinner,
	    		                                     sp.getString(ReporterActivity.PREF_SYNC_FREQ,
	    		                                    		      upSpinner.getItemAtPosition(0).toString())));
	    upSpinner.setOnItemSelectedListener(this);
	    
	    Spinner rangeSpinner = (Spinner) rootView.findViewById(R.id.summary_scope);
	    ArrayAdapter<CharSequence> rAdapter = ArrayAdapter.createFromResource(this.getActivity(),
	            R.array.pref_summaryScope_entries, android.R.layout.simple_spinner_item);
	    rAdapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
	    rangeSpinner.setAdapter(rAdapter);
	    rangeSpinner.setSelection(this.findItemPosition(rangeSpinner, 
	    		                                        sp.getString(ReporterActivity.PREF_SUMMARY_SCOPE, 
	    		                                                     rangeSpinner.getItemAtPosition(0).toString())));
	    rangeSpinner.setOnItemSelectedListener(this);
	    
	    Spinner alertSpinner = (Spinner) rootView.findViewById(R.id.notifications_settings);
	    ArrayAdapter<CharSequence> alrtAdapter = ArrayAdapter.createFromResource(this.getActivity(),
	            R.array.pref_notificationSettings_entries, android.R.layout.simple_spinner_item);
	    alrtAdapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
	    alertSpinner.setAdapter(alrtAdapter);
	    alertSpinner.setSelection(this.findItemPosition(alertSpinner,
	    		                                        sp.getString(ReporterActivity.PREF_NOTIFICATIONS,
	    		                                        		     alertSpinner.getItemAtPosition(0).toString())));
	    alertSpinner.setOnItemSelectedListener(this);
	    
	    EditText hostText = (EditText) rootView.findViewById(R.id.remote_server_url);
	    hostText.setText(sp.getString(ReporterActivity.PREF_HOSTNAME, null));
	    this.hostUrl = hostText.getText().toString();
	    
	    //Disable/Enable controls depending on host setting
	    final Spinner[] sps = {upSpinner,rangeSpinner,alertSpinner};
	    if (null != hostText.getText() && !hostText.getText().toString().isEmpty()) {
	    	this.toggleSpinners(sps, true);
	    } else {
	    	this.toggleSpinners(sps, false);
	    }
	    
    	hostText.addTextChangedListener(new TextWatcher() {

			@Override
			public void beforeTextChanged(CharSequence s, int start, int count,
					int after) {}

			@Override
			public void onTextChanged(CharSequence s, int start, int before,
					int count) {}

			@Override
			public void afterTextChanged(Editable s) {
				
				hostUrl = s.toString();
				if (hostUrl.isEmpty() && sps[0].isEnabled()) {
				  toggleSpinners(sps, false);
				} else if (!sps[0].isEnabled()) {
				  toggleSpinners(sps, true);
				}
			}});
    	
    	hostText.setActivated(false);
    	Dialog settingsDialog = builder.create();
    	settingsDialog.getWindow().setSoftInputMode(WindowManager.LayoutParams.SOFT_INPUT_STATE_HIDDEN);
    	return settingsDialog;
    }

	private void toggleSpinners(Spinner[] spinners, boolean enable) {
		
		for (Spinner sp : spinners) {
			sp.setEnabled(enable);
		}
	}
	
	private int findItemPosition (Spinner spinner, String entry) {
	 	for (int p = 0; p < spinner.getCount(); p++) {
	 		if (spinner.getItemAtPosition(p).toString().equals(entry)) {
	 		    return p;
	 		}
	 	}
		return 0;
	}
	
    @Override
    public void onAttach(Activity activity) {
        super.onAttach(activity);

        try {
            mListener = (PrefDialogListener) activity;
        } catch (ClassCastException e) {
            throw new ClassCastException(activity.toString()
                    + " must implement ConfigureDialogListener");
        }
    }

	@Override
	public void onItemSelected(AdapterView<?> parent, View view, int position,
			long id) {

		switch(parent.getId()) {
		case R.id.automatic_updates_settings:
			this.updateFrequency = parent.getItemAtPosition(position).toString();
			Log.d(ReporterActivity.TAG, "Update Range set to " + this.updateFrequency);
		break;
		case R.id.summary_scope:
			this.timeRange = parent.getItemAtPosition(position).toString();
			Log.d(ReporterActivity.TAG, "Time Range set to " + this.timeRange);
		break;
		case R.id.notifications_settings:
			this.notificationType = parent.getItemAtPosition(position).toString();
			Log.d(ReporterActivity.TAG, "Notification type set to " + this.notificationType);
        break;
		default:
		break;
		};
	}

	@Override
	public void onNothingSelected(AdapterView<?> parent) {
		// Auto-generated method stub
		
	}

}
