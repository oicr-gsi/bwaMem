package ca.on.oicr.pde.seqprodreporter;

import android.content.Context;
import android.content.Intent;
import android.os.Message;
import android.support.v7.widget.SearchView;
import android.util.AttributeSet;
import android.util.Log;

public class SeqprodSearchView extends SearchView {
	private String searchFilter;

	// Get search filter for this SearchView
	public String getSearchFilter() {
		return searchFilter;
	}

	// Set search filter for this SearchView
	public void setSearchFilter(String searchFilter) {
		this.searchFilter = searchFilter;
	}

	public SeqprodSearchView(Context context) {
		super(context);
		searchFilter = null;
	}

	public SeqprodSearchView(Context context, AttributeSet attrs) {
		super(context, attrs);
		searchFilter = null;
	}
	
	@Override
	public void onActionViewCollapsed() {
		super.onActionViewCollapsed();
		Log.d(ReporterActivity.TAG,"Unsetting search filter");
		this.setQuery("", false);
		this.setSearchFilter(null);
		Intent searchChange = new Intent(getContext(), ReporterActivity.class);
		searchChange.setAction(Intent.ACTION_SEARCH);
		//searchChange.set
		getContext().startActivity(searchChange);
		//Message msg =  ReporterActivity new Message();
		
	}

}