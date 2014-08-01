package ca.on.oicr.pde.seqprodreporter;

import android.content.Context;
import android.content.Intent;
import android.support.v7.widget.SearchView;
import android.util.Log;

/**
 * 
 * @author pruzanov
 *
 * This class provides the following behavior:
 * On collapsing the search input box it will always reset the filter,
 * So if user wants to keep the search filter enabled, he needs to 
 * keep the widget open, otherwise it will reset the search query
 */
public class SeqprodSearchView extends SearchView {

	private String mSearchString;

	@Override
	public void onActionViewExpanded() {
		super.onActionViewExpanded();
		this.setQuery(getSearchString(), false);
	}

	private final static String EMPTY_QUERY = "";

	// Default constructor
	public SeqprodSearchView(Context context) {
		super(context);
	}

	@Override
	public void onActionViewCollapsed() {
		super.onActionViewCollapsed();
		Log.d(ReporterActivity.TAG, "Unsetting search filter");
		this.setQuery(EMPTY_QUERY, false);
		this.mSearchString = EMPTY_QUERY;
		Intent searchReset = new Intent(getContext(), ReporterActivity.class);
		searchReset.setAction(Intent.ACTION_SEARCH);
		getContext().startActivity(searchReset);

	}

	public String getSearchString() {
		return mSearchString;
	}

	public void setSearchString(String mSearchQuery) {
		this.mSearchString = mSearchQuery;
	}

}