package org.apache.commons.math3.stat.descriptive.rank;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Random;

import org.junit.Test;

public class HighDynamicRangeQuantileTest {
	
	//@Ignore
	@Test
	public void testPerformanceRandom() {
	
		final double min = 434.;
		final double max = 4656.;
		final double numBins = 1000;
		
		final int M = 200;
		final int N = 5000000;
		
		final double values[] = new double[N];
		final Random random = new Random(0);
		for (int i = 0; i < N; ++i) {
			values[i] = min*Math.pow(max/min,random.nextDouble());
		}
		
		
		// Arrays.sort(values);
		
		final long  start = System.currentTimeMillis();
		for (int m = 0; m < M; ++m) {			
			final HighDynamicRangeQuantile quantile =  new HighDynamicRangeQuantile(min, max, Math.pow(max/min, 1./numBins));
			for (int i = 0; i < N; ++i) {
				quantile.add(values[i], 1);
			}
			// assertEquals(N, quantile.getCount());
		}
		final long end = System.currentTimeMillis();
		
		
		System.out.println("Avg time add operation unsorted data = " + ((end - start)*1e6)/(N*M) + "ns.");
		
	}
	
	private static final double searchTransitionToIndex(final int idx, final HighDynamicRangeQuantile quantile) {
		double low = Double.MIN_VALUE;
		double high = Double.MAX_VALUE;
		
		while(Math.nextUp(low)!=high) {
			double mid = (low+high)*0.5;
			if (mid>= high) {
				mid = Math.nextDown(high);
			}
			if (mid <= low) {
				mid = Math.nextUp(low);
			}
			final int i = quantile.getIndex(mid);
			if (i >= idx) {
				high = mid;
			}
			else {
				low = mid;
			}
		}
		double low2 = low;
		for(int i = 0; i < 100; ++i) {
			assertEquals(idx-1, quantile.getIndex(low2));
			low2 = Math.nextDown(low2);
		}
		
		double high2 = high;
		for(int i = 0; i < 100; ++i) {
			assertEquals(idx, quantile.getIndex(high2));
			high2 = Math.nextUp(high2);
		}
		return high;
	}

	@Test
	public void testIndex() {
		final double r = 1.000324;
		final double min = 434.;
		final double max = 4656.;
		final HighDynamicRangeQuantile quantile =  new HighDynamicRangeQuantile(min, max, r);
		
		double current = min;
		int idx = 0;
		while(current <= max) {
			idx += 1;
			final double last = current;
			current = searchTransitionToIndex(idx, quantile);
			
			//System.out.println((idx-1) + "   " + idx + "    " + last + "     " + current + "      " + current/last);
			assertTrue(current/last <= r);
		}			
	}
	
	
	

}
