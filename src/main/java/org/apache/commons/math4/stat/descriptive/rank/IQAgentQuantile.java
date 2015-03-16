package org.apache.commons.math4.stat.descriptive.rank;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;

import org.apache.commons.math4.exception.OutOfRangeException;
import org.apache.commons.math4.stat.descriptive.rank.MinHeapUtils.MinDoubleIntHeap;
import org.apache.commons.math4.util.MathArrays;
import org.apache.commons.math4.util.MathArrays.OrderDirection;
import org.apache.commons.math4.util.MathUtils;


public class IQAgentQuantile {
	
	private final double[] buffer;
	
	private int bufferCounter;
	
	private final double[] pValues;
	
	private double[] quantiles;
	
	private long totalCount;	
	
	// TODO exclude 0. and 1. from pValues
	
	/**
	 * @param pValues array of double values the first must be 0., the last must be 1.
	 * @param bufferSize buffer size, must be at least 1
	 */
	public IQAgentQuantile(final double[] pValues, final int bufferSize) {
		
		// TODO improve argument checks
		MathUtils.checkNotNull(pValues);
		MathArrays.checkOrder(pValues, OrderDirection.INCREASING, true);
		if (pValues.length < 2) {
			throw new IllegalArgumentException();
		}
		if (pValues[0] != 0.) {
			throw new IllegalArgumentException();
		}
		if (pValues[pValues.length-1] != 1.) {
			throw new IllegalArgumentException();
		}
		
		if (bufferSize<1) {
			throw new IllegalArgumentException();
		}
		
		
		this.pValues = pValues;
		this.buffer = new double[bufferSize];
	}
	
	static final Histogram asHistogram(final double[] pValues, final double[] quantiles, final long count, final double scale) {
		
		final double extremeValueWeight = 0.5/count;
		
		return new Histogram() {
			
			@Override
			public int getNumberOfBins() {
				return pValues.length+1;
			}
			
			@Override
			public double getFrequency(final int idx) {
				if (idx > 0 && idx < pValues.length) {
					return Math.max(0., Math.min(1. - extremeValueWeight, pValues[idx]) - Math.max(extremeValueWeight, pValues[idx-1]))*scale;
				} else {					
					return extremeValueWeight*scale;
				}
			}
			
			@Override
			public double getBoundary(final int idx) {
				if (idx <= 0) {
					return quantiles[0];
				}
				else if (idx > quantiles.length) {
					return quantiles[quantiles.length-1];
				}
				else {
					return quantiles[idx-1];
				}
			}
		};
	}
	
	private void collapse() {
		
		if (bufferCounter == 0) {
			return;
		}
		
		totalCount += bufferCounter;
		
		final Collection<Histogram> histograms = new ArrayList<Histogram>(2);
		
		final double fractionValuesInBuffer = ((double)bufferCounter)/totalCount;
		
		// add histogram of buffer values
		
		Arrays.sort(buffer, 0, bufferCounter);
		
		double maximum = buffer[bufferCounter-1]; 

		histograms.add(sortedValuesAsHistogram(buffer, bufferCounter, fractionValuesInBuffer));
		
		if (quantiles != null) {
			maximum = Math.max(maximum, quantiles[quantiles.length-1]);
			histograms.add(asHistogram(pValues, quantiles, totalCount-bufferCounter, 1.-fractionValuesInBuffer));
			
		}
		quantiles = evaluateSumOfHistograms(histograms, pValues);
		
		// ensure that last quantile is equal to maximum
		quantiles[quantiles.length-1] = maximum;
		
		bufferCounter = 0;
	}

	public void add(final double value) {
		MathUtils.checkFinite(value);
		buffer[bufferCounter] = value;
		bufferCounter += 1;
		if (bufferCounter == buffer.length) {
			collapse();
		}
	}
	
	public double getQuantile(final double pValue) {

		if (pValue < 0. || pValue > 1.) {
			throw new OutOfRangeException(pValue, 0., 1.);
		}
		
		if (totalCount == 0) {
			return Double.NaN;
		}
		
		collapse();

		
		// TODO use binary search
		
		final Histogram histogram = asHistogram(pValues, quantiles, totalCount-bufferCounter, 1.);
		
		return evaluateSumOfHistograms(Collections.singleton(histogram), new double[]{pValue})[0];
		
	}
	
	/**
	 * A histogram.
	 */
	interface Histogram {
		
		/**
		 * Returns the number of bins
		 * 
		 * @return the number of bins
		 */
		int getNumberOfBins();
		
		/**
		 * @param idx 0 <= idx <= {@link #getNumberOfBins()}
		 * @return
		 */
		public double getBoundary(int idx);
		
		/**
		 * The frequency between {@code getBoundary(idx)} and {@code getBoundary(idx+1)}. 
		 * 
		 * @param idx 0 <= idx < {@link #getNumberOfBins()}
		 * @return
		 */
		public double getFrequency(int idx);
		
	}
	
	interface HistogramIterator {
		
		boolean advance();
		
		/**
		 * Returns the minimum value of the current interval.
		 * 
		 * @return the minimum value of the current interval
		 */
		double getMinimumValue();

		/**
		 * Returns the maximum value of the current interval.
		 * 
		 * @return the maximum value of the current interval
		 */
		double getMaximumValue();
		
		double getMinimumCumulativeFrequency();
		
		double getMaximumCumulativeFrequency();
		
		double getDensity();
	}

	
	final static class DynamicSum {
		
		private final double partialSums[];
		
		public DynamicSum(final int n) {
			partialSums = new double[(n-1) << 1];
		}
		
		public double update(final int idx, double value) {
			int i = partialSums.length-idx;	
			while(i!=0) {
				i -= 1;
				partialSums[i] = value;
				value += partialSums[i ^ 1];
				i = i >>> 1;
			}
			return value;
		}
	}
	
	final static class HistogramIterator1 implements HistogramIterator {
			
		private final Histogram[] histograms;
		private final MinDoubleIntHeap heap;
		private final DynamicSum densitySum;
		private final int[] counters;
		private double minCumulativeFrequency;
		private double maxCumulativeFrequency;
		private double minValue;
		private double maxValue;
		private double density;

		public HistogramIterator1(final Collection<? extends Histogram> histograms) {
			
			final int numHistograms = histograms.size();
			
			this.histograms = histograms.toArray(new Histogram[numHistograms]);
			this.counters = new int[numHistograms];
			this.densitySum = new DynamicSum(numHistograms);
			this.density = 0.;
			
			final double yBegin = 0.;
			
			final double[] minValues = new double[numHistograms];
			for(int histogramIdx = 0; histogramIdx < numHistograms; ++histogramIdx) {
				final Histogram histogram = this.histograms[histogramIdx];
				minValues[histogramIdx] = histogram.getBoundary(0); 
			}
			
			this.heap = MinHeapUtils.createMinDoubleIntHeap(minValues);
			
			this.minCumulativeFrequency = yBegin;
			this.maxCumulativeFrequency = yBegin;
			this.minValue = Double.NEGATIVE_INFINITY;			
			this.maxValue = heap.getMinValue();
		}
		
		@Override
		public boolean advance() {
			minValue = maxValue;
			minCumulativeFrequency = maxCumulativeFrequency;
			
			final int histogramIdx = heap.getMinIndex();
			
			final int pos = counters[histogramIdx] + 1;
			counters[histogramIdx] = pos;
			
			final Histogram histogram = histograms[histogramIdx];
			
			if (pos <= histogram.getNumberOfBins()) {			
				final double nextPointX = histogram.getBoundary(pos);
				final double deltaY = histogram.getFrequency(pos-1);
				if (minValue < nextPointX) {
					density = densitySum.update(histogramIdx, deltaY/(nextPointX-minValue)); // TODO handle potential overflow
					heap.update(nextPointX);
					maxValue = heap.getMinValue();
					maxCumulativeFrequency += density*(maxValue - minValue);
				}
				else {
					density = Double.POSITIVE_INFINITY;
					maxCumulativeFrequency += deltaY;
				}
				return true;
			}
			else {
				heap.update(Double.POSITIVE_INFINITY);
				density = densitySum.update(histogramIdx, 0.0);
				maxValue = heap.getMinValue();
				if (maxValue != Double.POSITIVE_INFINITY) {
					maxCumulativeFrequency += density*(maxValue - minValue);
					return true;
				}
				else {
					return false;
				}
			}
		}
		
		@Override
		public double getMinimumValue() {
			return minValue;
		}

		@Override
		public double getMaximumValue() {
			return maxValue;
		}
		
		@Override
		public double getDensity() {
			return density;
		}
		
		@Override
		public double getMinimumCumulativeFrequency() {
			return minCumulativeFrequency;
		}
		
		@Override
		public double getMaximumCumulativeFrequency() {
			return maxCumulativeFrequency;
		}

		@Override
		public String toString() {
			return "HistogramIterator1 [histograms="
					+ Arrays.toString(histograms) + ", heap=" + heap
					+ ", gradientSum=" + densitySum + ", counters="
					+ Arrays.toString(counters) + ", minY=" + minCumulativeFrequency + ", maxY="
					+ maxCumulativeFrequency + ", minX=" + minValue + ", maxX=" + maxValue
					+ ", gradient=" + density + "]";
		}
		
		
	}
	
	/**
	 * Returns a {@link Histogram} representing the cumulative distribution function of
	 * a given array of double values in ascending order. 
	 * The returned histogram is backed by the array.
	 * 
	 * @param values an array of double values in ascending order with length at least 1
	 * @param size defines the number of leading array elements that are used for the histogram
	 * @return a histogram backed by the given array
	 */
	static final Histogram sortedValuesAsHistogram(final double[] values, final int size, final double scale) {
		
		final double increment = scale/size;
		
		return new Histogram() {
			
			@Override
			public double getBoundary(final int idx) {
				return values[idx >>> 1];
			}
			
			@Override
			public int getNumberOfBins() {
				return (size<<1)-1;
			}

			@Override
			public double getFrequency(final int idx) {
				return ((idx & 1) == 0)?increment:0.;
			}
		};
	}
	
	
	static final Histogram cumulativeFrequenciesAsHistogram(final double[] cumulativeFrequencies, final double[] binBoundaries, final double scale) {
		
		final double factor = scale/cumulativeFrequencies[cumulativeFrequencies.length-1];
		
		return new Histogram() {

			@Override
			public int getNumberOfBins() {
				return cumulativeFrequencies.length;
			}

			@Override
			public double getBoundary(final int idx) {
				return binBoundaries[idx];
			}

			@Override
			public double getFrequency(final int idx) {
				double result = cumulativeFrequencies[idx];
				if (idx >0 ) {
					result -= cumulativeFrequencies[idx-1];
				}
				return result*factor; 
			}
		};
	}
	
	/**
	 * Get value between two given points using linear interpolation.
	 * 
	 * @param x1 x-value of point 1
	 * @param y1 y-value of point 1
	 * @param x2 x-value of point 2 (x1 < x2)
	 * @param y2 y-value of point 2 (y1 < y2)
	 * @param x x-value for which the corresponding y-value needs to interpolated 
	 * @return the interpolated value, is always in the range [y1, y2]
	 */
	static double interpolate(final double x1, final double y1, final double x2, final double y2, final double x) {
		final double alpha = (x - x1)/(x2 - x1);
		final double result = y1 + (y2 - y1)*alpha;
		return (result <= y2)?result:y2;
	}
	
	
	static final double[] evaluateSumOfHistograms(final Collection<? extends Histogram> histograms, final double[] cumulativeFrequencies) {
		
		final double[] values = new double[cumulativeFrequencies.length];
		
		final HistogramIterator iterator = new HistogramIterator1(histograms);

		int valueCounter = 0;
		
		while(valueCounter < cumulativeFrequencies.length && cumulativeFrequencies[valueCounter] <= iterator.getMinimumCumulativeFrequency()) {
			values[valueCounter] = iterator.getMaximumValue();
			valueCounter += 1;
		}
		
		while(valueCounter < cumulativeFrequencies.length && iterator.advance()) {
			
			if (iterator.getMinimumCumulativeFrequency() != iterator.getMaximumCumulativeFrequency()) {
				if (cumulativeFrequencies[valueCounter] == iterator.getMinimumCumulativeFrequency()) {
					values[valueCounter]+=iterator.getMinimumValue();
					values[valueCounter]*=0.5;
					valueCounter += 1;
				}
				
				while(valueCounter < cumulativeFrequencies.length && cumulativeFrequencies[valueCounter] < iterator.getMaximumCumulativeFrequency()) {
					values[valueCounter] = interpolate(iterator.getMinimumCumulativeFrequency(), iterator.getMinimumValue(), iterator.getMaximumCumulativeFrequency(), iterator.getMaximumValue(), cumulativeFrequencies[valueCounter]);
					valueCounter += 1;
				}
				if (valueCounter < cumulativeFrequencies.length && cumulativeFrequencies[valueCounter] == iterator.getMaximumCumulativeFrequency()) {
					values[valueCounter] = iterator.getMaximumValue();
				}
			}
		}
		
		while(valueCounter < cumulativeFrequencies.length) {
			values[valueCounter] = iterator.getMinimumValue();
			valueCounter += 1;
		}
	
		return values;
	}

	@Override
	public String toString() {
		return "IQAgentQuantile [buffer=" + Arrays.toString(buffer)
				+ ", bufferCounter=" + bufferCounter + ", pValues="
				+ Arrays.toString(pValues) + ", quantiles="
				+ Arrays.toString(quantiles) + ", totalCount=" + totalCount
				+ "]";
	}
	
	
}
