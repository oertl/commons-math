package org.apache.commons.math4.stat.descriptive.rank;

import java.util.concurrent.ThreadLocalRandom;

public abstract class AbstractQuantileEstimator implements QuantileEstimator {

	@Override
	public void add(final double value) {
		add(value, 1);
	}

	@Override
	public void add(final double value, final long count) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double quantile(final double pValue) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double probability(final double from, final double to) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public long frequency(final double from, final double to) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double minimum() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double maximum() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public long count() {
		// TODO Auto-generated method stub
		return 0;
	}
	
	public void merge() {
		
	}

	public static final void sort(final double[] values, final long[] weights, final int size) {
		quickSort(values, weights, 0, size, ThreadLocalRandom.current());
	}
	
	private static void quickSort(final double[] values, final long[] weights, final int from, final int to, final ThreadLocalRandom random) {
		if (to - from <= 1) {
            // sorted by definition
            return;
        }
        final int p = partition(values, weights, from, to, random);
        quickSort(values, weights, from, p, random);
        quickSort(values, weights, p + 1, to, random);
		
	}

	private static int partition(final double[] values, final long[] weights, final int from, final int to, final ThreadLocalRandom random) {
		final int pivotIndex = from + random.nextInt(to - from);
        final double pivotValue = values[pivotIndex];
        swap(values, weights, pivotIndex, to - 1);
        int p = from;
        for (int i = from; i < to - 1; ++i) {
            if (values[i] < pivotValue) {
                swap(values, weights, i, p++);
            }
        }
        swap(values, weights, p, to - 1);
        return p;
	}

	private static final void swap(final double[] values, final long[] weights, final int i, final int j) {
		{
	       	final double tmp = values[i];
	        values[i] = values[j];
	        values[j] = tmp;			
		}
		{
	       	final long tmp = weights[i];
	       	weights[i] = weights[j];
	       	weights[j] = tmp;			
		}
    }
}
