package org.apache.commons.math3.stat.descriptive.rank;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;

import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.MathUtils;
import org.apache.commons.math3.util.MathArrays.OrderDirection;

public class IQAgentQuantile {
	
	private final double[] buffer;
	
	private int bufferCounter;
	
	private final double[] pValues;
	
	private double[] quantiles;
	
	private long totalCount;	
	
	/**
	 * @param pValues array of double values the first must be 0., the last must be 1.
	 * @param bufferSize buffer size, must be at least 1
	 */
	public IQAgentQuantile(double[] pValues, int bufferSize) {
		
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
			
			public int getNumberOfBins() {
				return pValues.length+1;
			}
			
			public double getFrequency(int idx) {
				if (idx > 0 && idx < pValues.length) {
					return Math.max(0., Math.min(1. - extremeValueWeight, pValues[idx]) - Math.max(extremeValueWeight, pValues[idx-1]))*scale;
				} else {					
					return extremeValueWeight*scale;
				}
			}
			
			public double getBoundary(int idx) {
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
		
		final Collection<Histogram> histograms = new ArrayList<Histogram>(2);
		
		final double fractionValuesInBuffer = ((double)bufferCounter)/totalCount;
		
		// add histogram of buffer values
		
		// TODO
		double[] bufferCopy = Arrays.copyOf(buffer, bufferCounter);
		Arrays.sort(bufferCopy);
		histograms.add(sortedValuesAsHistogram(bufferCopy, fractionValuesInBuffer));
		
		if (quantiles != null) {
			histograms.add(asHistogram(pValues, quantiles, totalCount-bufferCounter, 1.-fractionValuesInBuffer));
		}
		quantiles = evaluateSumOfHistograms(histograms, pValues);
		
		bufferCounter = 0;
	}

	public void add(double value) {
		
		buffer[bufferCounter] = value;
		bufferCounter += 1;
		totalCount += 1;
		
		if (bufferCounter == buffer.length) {
			collapse();
		}
	}
	
	public double getQuantile(double pValue) {
	
		// TODO improve argument checking
		if (pValue < 0. || pValue > 1) {
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
		
		public DynamicSum(int n) {
			partialSums = new double[(n-1) << 1];
		}
		
		public double update(int idx, double value) {
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
	
	final static class MinDoubleIntHeap {
		
		private final double[] values;
		private final int[] indices;
		
		public MinDoubleIntHeap(double[] values) {
			int n = values.length;
			
			this.values = new double[n];
			this.indices = new int[n];
			for (int i = 0; i < n; ++i) {
				this.values[i] = Double.NEGATIVE_INFINITY;
			}
			for (int i = 0; i < n; ++i) {
				this.indices[0] = i;
				update(values[i]);
			}
		}
		
		public int getMinIndex() {
			return indices[0];
		}
		
		public double getMinValue() {
			return values[0];
		}
		
		public void update(double newValue) {
			
			final int n = values.length;
			
			final int updatedIndex = indices[0];
			
			int parentIdx = 0;
			
			while(true) {
				
				int leftChildIdx = (parentIdx << 1) | 1;			
				
				if (leftChildIdx < n-1) {
					
					final int rightChildIdx = leftChildIdx+1;
					
					final double leftChildValue = values[leftChildIdx];
					final double rightChildValue = values[rightChildIdx];
					
					final int minChildIdx;
					final double minChildValue;
					if (rightChildValue < leftChildValue) {
						minChildIdx = rightChildIdx;
						minChildValue = rightChildValue;
					}
					else {
						minChildIdx = leftChildIdx;
						minChildValue = leftChildValue;
					}
					if (minChildValue < newValue) {
						values[parentIdx] = minChildValue;
						indices[parentIdx] = indices[minChildIdx];
						parentIdx = minChildIdx;
						continue;
					}
				} else if (leftChildIdx == n-1) {
					final double leftChildValue = values[leftChildIdx];
					if (leftChildValue < newValue) {
						values[parentIdx] = leftChildValue;
						indices[parentIdx] = indices[leftChildIdx];
						parentIdx = leftChildIdx;
					}
				}
				break;
			}
			values[parentIdx] = newValue;
			indices[parentIdx] = updatedIndex;		
		}
	}
	
	final static class MinDoubleIntHeap2 {
		
		private final double[] values;
		private final int[] indices;
		
		public MinDoubleIntHeap2(double[] values) {
			int nt = values.length;
			
			int n = (Integer.highestOneBit(nt)<<1)-1;
			
			
			
			this.values = new double[n];
			this.indices = new int[n];
			for (int i = 0; i < nt; ++i) {
				this.values[i] = Double.NEGATIVE_INFINITY;
			}
			for (int i = nt; i < n; ++i) {
				this.values[i] = Double.NaN;
			}
			for (int i = 0; i < nt; ++i) {
				this.indices[0] = i;
				update(values[i]);
			}
		}
		
		public int getMinIndex() {
			return indices[0];
		}
		
		public double getMinValue() {
			return values[0];
		}
		
		public void update(double newValue) {
			
			final int n = values.length;
			
			final int updatedIndex = indices[0];
			
			int parentIdx = 0;
			
			for(int i = n; i!=1; i>>=1 ) {
				
				final int leftChildIdx = (parentIdx << 1) | 1;			
				final int rightChildIdx = leftChildIdx+1;
				
				final double leftChildValue = values[leftChildIdx];
				final double rightChildValue = values[rightChildIdx];
				
				final int minChild;
				final double minChildValue;
				if (rightChildValue < leftChildValue) {
					minChild = rightChildIdx;
					minChildValue = rightChildValue;
				}
				else {
					minChild = leftChildIdx;
					minChildValue = leftChildValue;
				}
				if (minChildValue < newValue) {
					values[parentIdx] = minChildValue;
					indices[parentIdx] = indices[minChild];
					parentIdx = minChild;
				}
				else {
					break;
				}
			}
			values[parentIdx] = newValue;
			indices[parentIdx] = updatedIndex;		
		}
	}
	
	final static class HistogramIterator1 implements HistogramIterator {
			
		private final Histogram[] histograms;
		private final MinDoubleIntHeap heap;
		private final DynamicSum gradientSum;
		private final int[] counters;
		private double minY;
		private double maxY;
		private double minX;
		private double maxX;
		private double gradient;

		public HistogramIterator1(final Collection<? extends Histogram> histograms) {
			
			final int numHistograms = histograms.size();
			
			this.histograms = histograms.toArray(new Histogram[numHistograms]);
			this.counters = new int[numHistograms];
			this.gradientSum = new DynamicSum(numHistograms);
			this.gradient = 0.;
			
			double yBegin = 0.;
			
			double[] minValues = new double[numHistograms];
			for(int histogramIdx = 0; histogramIdx < numHistograms; ++histogramIdx) {
				final Histogram histogram = this.histograms[histogramIdx];
				minValues[histogramIdx] = histogram.getBoundary(0); 
			}
			
			this.heap = new MinDoubleIntHeap(minValues);
			
			this.minY = yBegin;
			this.maxY = yBegin;
			this.minX = Double.NEGATIVE_INFINITY;			
			this.maxX = heap.getMinValue();
		}
		
		public boolean advance() {
			minX = maxX;
			minY = maxY;
			
			final int histogramIdx = heap.getMinIndex();
			
			final int pos = counters[histogramIdx] + 1;
			counters[histogramIdx] = pos;
			
			final Histogram histogram = histograms[histogramIdx];
			
			if (pos <= histogram.getNumberOfBins()) {			
				final double nextPointX = histogram.getBoundary(pos);
				final double deltaY = histogram.getFrequency(pos-1);
				if (minX < nextPointX) {
					gradient = gradientSum.update(histogramIdx, deltaY/(nextPointX-minX)); // TODO handle potential overflow
					heap.update(nextPointX);
					maxX = heap.getMinValue();
					maxY += gradient*(maxX - minX);
				}
				else {
					gradient = Double.POSITIVE_INFINITY;
					maxY += deltaY;
				}
				return true;
			}
			else {
				heap.update(Double.POSITIVE_INFINITY);
				gradient = gradientSum.update(histogramIdx, 0.0);
				maxX = heap.getMinValue();
				if (maxX != Double.POSITIVE_INFINITY) {
					maxY += gradient*(maxX - minX);
					return true;
				}
				else {
					return false;
				}
			}
		}
		
		public double getMinimumValue() {
			return minX;
		}

		public double getMaximumValue() {
			return maxX;
		}
		
		public double getDensity() {
			return gradient;
		}
		
		public double getMinimumCumulativeFrequency() {
			return minY;
		}
		
		public double getMaximumCumulativeFrequency() {
			return maxY;
		}

		@Override
		public String toString() {
			return "HistogramIterator1 [histograms="
					+ Arrays.toString(histograms) + ", heap=" + heap
					+ ", gradientSum=" + gradientSum + ", counters="
					+ Arrays.toString(counters) + ", minY=" + minY + ", maxY="
					+ maxY + ", minX=" + minX + ", maxX=" + maxX
					+ ", gradient=" + gradient + "]";
		}
		
		
	}
	
	/**
	 * Returns a {@link Histogram} representing the cumulative distribution function of
	 * a given array of double values in ascending order. 
	 * The returned histogram is backed by the array.
	 * 
	 * @param values an array of double values in ascending order with length at least 1
	 * @return a histogram backed by the given array
	 */
	static final Histogram sortedValuesAsHistogram(final double[] values, final double scale) {
		
		final double increment = scale/values.length;
		
		return new Histogram() {
			
			public double getBoundary(int idx) {
				return values[idx >>> 1];
			}
			
			public int getNumberOfBins() {
				return (values.length<<1)-1;
			}

			public double getFrequency(int idx) {
				return ((idx & 1) == 0)?increment:0.;
			}
		};
	}
	
	
	static final Histogram cumulativeFrequenciesAsHistogram(final double[] cumulativeFrequencies, final double[] binBoundaries, final double scale) {
		
		final double factor = scale/cumulativeFrequencies[cumulativeFrequencies.length-1];
		
		return new Histogram() {

			public int getNumberOfBins() {
				return cumulativeFrequencies.length;
			}

			public double getBoundary(int idx) {
				return binBoundaries[idx];
			}

			public double getFrequency(int idx) {
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
	 * @param y2 y-value of point 2
	 * @param x x-value for which the corresponding y-value needs to interpolated 
	 * @return the interpolated value, is always in the range [y1, y2]
	 */
	static double interpolate(double x1, double y1, double x2, double y2, double x) {
		double alpha = (x - x1)/(x2 - x1);
		double result = y1 + (y2 - y1)*alpha;
		return (result <= y2)?result:y2;
	}
	
	
	static final double[] evaluateSumOfHistograms(final Collection<? extends Histogram> histograms, double[] cumulativeFrequencies) {
		
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
