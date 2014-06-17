package ca.on.oicr.pde.seqprodreporter;

import android.content.Context;
import android.content.Intent;
import android.support.v7.widget.SearchView;
import android.util.Log;

public class SeqprodSearchView extends SearchView {
    private final String EMPTY_QUERY = "";
	
    // Default constructor
    public SeqprodSearchView(Context context) {
		super(context);
	}

	@Override
	public void onActionViewCollapsed() {
		super.onActionViewCollapsed();
		Log.d(ReporterActivity.TAG, "Unsetting search filter");
		this.setQuery(EMPTY_QUERY, false);
		Intent searchReset = new Intent(getContext(), ReporterActivity.class);
		searchReset.setAction(Intent.ACTION_SEARCH);
		getContext().startActivity(searchReset);
	}

}