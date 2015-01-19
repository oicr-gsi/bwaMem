/*
 * Copyright (C) 2013 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package ca.on.oicr.pde.seqprodreporter;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.support.v4.app.FragmentPagerAdapter;
import android.support.v4.view.ViewPager;
import android.text.format.Time;
import android.util.Log;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.view.ViewGroup.OnHierarchyChangeListener;
import ca.on.oicr.pde.seqprodreporter.ui.SlidingTabLayout;


/**
 * A test view that implements tab UI with sliding tabs (a-la Google tabs)
 */
public class SlidingTabsBasicFragment extends android.support.v4.app.Fragment {

    static final String LOG_TAG = "SlidingTabsBasicFragment";
    //static final String[] tabTypes = {"COMPLETED", " FAILED ", " PENDING"}; 
    
    private String searchFilter;
    private int mCurrentTabIndex;
    private int mSortIndex;

	private void setSortIndex(int index) {
		this.mSortIndex = index;
	}

	protected void setCurrentTabIndex(int index) {
		this.mCurrentTabIndex = index;
	}

	/**
     * A custom {@link ViewPager} title strip which looks much like Tabs present in Android v4.0 and
     * above, but is designed to give continuous feedback to the user when scrolling.
     */
    private SlidingTabLayout mSlidingTabLayout;

    /**
     * A {@link ViewPager} which will be used in conjunction with the {@link SlidingTabLayout} above.
     */
    private ViewPager mViewPager;
    
    /**
     * A {@link SamplePagerAdapter} which will be used in conjunction with the {@link ViewPager} above.
     */
    private SamplePagerAdapter mFragmentAdapter;

    /**
     * Inflates the {@link View} which will be displayed by this {@link Fragment}, from the app's
     * resources.
     */
    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container,
            Bundle savedInstanceState) {
        return inflater.inflate(R.layout.fragment_tab, container, false);
    }
    
    public static SlidingTabsBasicFragment newInstance(int selectedTab, int sortIndex, String searchFilter) {
    	SlidingTabsBasicFragment fragment = new SlidingTabsBasicFragment();
		fragment.setCurrentTabIndex(selectedTab);
		fragment.updateSearchFilter(searchFilter);
		fragment.setSortIndex(sortIndex);
		return fragment;
	}
    
    /**
     * Service method for updating search filter
     */
    public void updateSearchFilter(String filter) {
    	this.searchFilter = filter;
    	if (null != this.mViewPager) {
    		this.mFragmentAdapter.updateFilter(filter);
    		this.updateUI();
    	}
    }
    
    /**
     * Service method for updating sort index
     */
    public void updateSortIndex(int sortIndex) {
    	this.setSortIndex(sortIndex);   	
    	if (null != this.mViewPager) {
    		this.mFragmentAdapter.sortLists(sortIndex);
    		this.updateUI();
    	} 
    }
    
    protected void updateUI() {
    	this.mFragmentAdapter.notifyDataSetChanged();
    }
    
    protected Time getFirstUpdateTime() {
    	try {
    	    ReportListFragment f = (ReportListFragment) this.mFragmentAdapter.getItem(0);
    	    return f.getFirstUpdateTime();
    	} catch (Exception e) {
    		Log.e(ReporterActivity.TAG, "Error retrieving First Update Time");
    		return null;
    	}
    }

    // BEGIN_INCLUDE (fragment_onviewcreated)
    /**
     * This is called after the {@link #onCreateView(LayoutInflater, ViewGroup, Bundle)} has finished.
     * Here we can pick out the {@link View}s we need to configure from the content view.
     *
     * We set the {@link ViewPager}'s adapter to be an instance of {@link SamplePagerAdapter}. The
     * {@link SlidingTabLayout} is then given the {@link ViewPager} so that it can populate itself.
     *
     * @param view View created in {@link #onCreateView(LayoutInflater, ViewGroup, Bundle)}
     */
    @Override
    public void onViewCreated(View view, Bundle savedInstanceState) {

        
        mFragmentAdapter = new SamplePagerAdapter(getFragmentManager(), SlidingTabsBasicFragment.this.searchFilter);

        mViewPager = (ViewPager) view.findViewById(R.id.viewpager);
        mViewPager.setOffscreenPageLimit(ReporterActivity.types.length - 1);
        mViewPager.setAdapter(mFragmentAdapter);
        
        mSlidingTabLayout = (SlidingTabLayout) view.findViewById(R.id.sliding_tabs);
        mSlidingTabLayout.setViewPager(mViewPager);
        
        mViewPager.setOnHierarchyChangeListener(new OnHierarchyChangeListener() {

			@Override
			public void onChildViewAdded(View parent, View child) {
				if (mCurrentTabIndex != 0
						&& mViewPager.getChildCount() > mCurrentTabIndex
						&& mCurrentTabIndex != mViewPager
								.getCurrentItem())
					mViewPager.setCurrentItem(mCurrentTabIndex);
			}

			@Override
			public void onChildViewRemoved(View parent, View child) {
				// Do nothing, this is not supposed to happen
			}
		});
        
    }
    
    

    /**
     * The {@link android.support.v4.view.PagerAdapter} used to display pages in this sample.
     * The individual pages are simple and just display two lines of text. The important section of
     * this class is the {@link #getPageTitle(int)} method which controls what is displayed in the
     * {@link SlidingTabLayout}.
     */
    class SamplePagerAdapter extends FragmentPagerAdapter { 	
    	private List<ReportListFragment> fragments;
    	private String searchFilter;
    	
    	public void updateFilter (String filter) {
    		this.searchFilter = filter;
    		if (!fragments.isEmpty()) {
    			for (int i = 0; i < fragments.size(); i++) {
    				ReportListFragment tmp = fragments.get(i);
    				if (null != tmp)
    					tmp.setSearchFilter(this.searchFilter);
    			}
    		}
    	}
    	
    	public void sortLists (int sortIndex) {
			
    		if (!fragments.isEmpty()) {
    			for (int i = 0; i < fragments.size(); i++) {
    				ReportListFragment tmpf = fragments.get(i);
    				if (null != tmpf) {
    					tmpf.setSortIndex(sortIndex);
    					tmpf.sortFragment();
    				}
    			}
    		}
    	}
    	
        public SamplePagerAdapter(android.support.v4.app.FragmentManager fragmentManager, String filter) {
			super(fragmentManager);
			fragments = new ArrayList<ReportListFragment>();
			this.searchFilter = filter;
		}

        @Override
        public int getCount() {
            return 3;
        }

        @Override
        public CharSequence getPageTitle(int position) {
            return ReporterActivity.types[position].toUpperCase(Locale.CANADA);
        }
        

		@Override
		public Fragment getItem(int position) {

			if (position >= this.fragments.size() || this.fragments.size() == 0
					|| null == this.fragments.get(position)) {

				fragments.add(position, ReportListFragment.newInstance(position, this.searchFilter));
			}

			ReportListFragment fragment = fragments.get(position);
			fragment.setSearchFilter(this.searchFilter);
			if (mSortIndex >= 0) {
				fragment.setSortIndex(mSortIndex);
			}
            return fragments.get(position);
		}
		
		@Override
		public int getItemPosition(Object object) {
		    return POSITION_NONE;
		}

    }
}
