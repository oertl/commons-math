package org.apache.commons.math3.stat.descriptive.rank;

import java.util.concurrent.ThreadLocalRandom;

public abstract class AbstractQuantileEstimator implements QuantileEstimator {

	public void add(double value) {
		add(value, 1);
	}

	public void add(double value, long count) {
		// TODO Auto-generated method stub
		
	}

	public double quantile(double pValue) {
		// TODO Auto-generated method stub
		return 0;
	}

	public double probability(double from, double to) {
		// TODO Auto-generated method stub
		return 0;
	}

	public long frequency(double from, double to) {
		// TODO Auto-generated method stub
		return 0;
	}

	public double minimum() {
		// TODO Auto-generated method stub
		return 0;
	}

	public double maximum() {
		// TODO Auto-generated method stub
		return 0;
	}

	public long count() {
		// TODO Auto-generated method stub
		return 0;
	}
	
	public void merge() {
		
	}

	public static final void sort(double[] values, long[] weights, int size) {
		quickSort(values, weights, 0, size, ThreadLocalRandom.current());
	}
	
	private static void quickSort(double[] values, long[] weights, int from, int to, ThreadLocalRandom random) {
		if (to - from <= 1) {
            // sorted by definition
            return;
        }
        final int p = partition(values, weights, from, to, random);
        quickSort(values, weights, from, p, random);
        quickSort(values, weights, p + 1, to, random);
		
	}

	private static int partition(double[] values, long[] weights, int from, int to, ThreadLocalRandom random) {
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

	private static final void swap(double[] values, long[] weights, int i, int j) {
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
