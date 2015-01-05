package org.apache.commons.math3.stat.descriptive.rank;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

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
		this.pValues = pValues;
		this.buffer = new double[bufferSize];
	}

	public void add(double value) {
		
		buffer[bufferCounter] = value;
		bufferCounter += 1;
		totalCount += 1;
				
		if (bufferCounter == buffer.length) {
			
			final Collection<Histogram> histograms = new ArrayList<Histogram>(2);
			
			final double fractionValuesInBuffer = ((double)bufferCounter)/totalCount;
			
			// add histogram of buffer values
			Arrays.sort(buffer);	
			histograms.add(sortedValuesAsHistogram(buffer, fractionValuesInBuffer));
			
			if (quantiles != null) {
			
				// add 
				int numberOfBins = pValues.length+1;
				
				final double[] binBoundaries = new double[numberOfBins+1];
				binBoundaries[0] = quantiles[0];
				for (int i = 1; i < binBoundaries.length-1; ++i) {
					binBoundaries[i] = quantiles[i-1];
				}
				binBoundaries[binBoundaries.length-1] = quantiles[quantiles.length-1];
				
				final double[] cumulativeFrequencies = new double[numberOfBins];
				final double cumulativeFrequencyMin = 0.5/totalCount;
				final double cumulativeFrequencyMax = 1.-cumulativeFrequencyMin;
				
				cumulativeFrequencies[0] = cumulativeFrequencyMin;
				for (int i = 1; i < cumulativeFrequencies.length-1; ++i) {
					double pValue = pValues[i-1]; 
					if (pValue < cumulativeFrequencyMin) {
						pValue = cumulativeFrequencyMin;
					} else if (pValue > cumulativeFrequencyMax) {
						pValue = cumulativeFrequencyMax;
					}
					cumulativeFrequencies[i] = pValue;  
				}
				cumulativeFrequencies[cumulativeFrequencies.length-1] = cumulativeFrequencyMax;
				histograms.add(cumulativeFrequenciesAsHistogram(cumulativeFrequencies, binBoundaries, 1.-fractionValuesInBuffer));
				
			}
			quantiles = evaluateSumOfHistograms(histograms, pValues);
			bufferCounter = 0;
		}
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
	
	final static class TwoHistogramsIterator implements HistogramIterator {
		
		private final Histogram histogram0;
		private final Histogram histogram1;
		private double next0;
		private double next1;
		private int nextIndex;
		private double gradient0;
		private double gradient1;
		private int counter0;
		private int counter1;
		private double minY;
		private double maxY;
		private double minX;
		private double maxX;
		private double gradient;

		public TwoHistogramsIterator(Histogram histogram0, Histogram histogram1) {
			
			this.histogram0 = histogram0;
			this.histogram1 = histogram1;
			this.counter0 = 0;
			this.counter1 = 0;
			this.gradient0 = 0.;
			this.gradient1 = 0.;
			this.gradient = 0.;
			
			double yBegin = 0.;
			
			this.next0 = histogram0.getBoundary(0);
			this.next1 = histogram1.getBoundary(0);
			if (next0 <= next1) {
				this.maxX = next0;
				this.nextIndex = 0;
			}
			else {
				this.maxX = next1;
				this.nextIndex = 1;
			}
			
			this.minY = yBegin;
			this.maxY = yBegin;
			this.minX = Double.NEGATIVE_INFINITY;			
		}
		
		public boolean advance() {
			minX = maxX;
			minY = maxY;
			
			final int pos;
			final Histogram histogram;
			if (nextIndex == 0) {
				counter0 += 1;
				pos = counter0;
				histogram = histogram0;
			}
			else {
				counter1 += 1;
				pos = counter1;
				histogram = histogram1;
				
			}
			
			if (pos <= histogram.getNumberOfBins()) {			
				final double nextPointX = histogram.getBoundary(pos);
				final double weightedDeltaY = histogram.getFrequency(pos-1);
				if (minX < nextPointX) {
					if (nextIndex==0) {
						gradient0 = weightedDeltaY/(nextPointX-minX);
						next0 = nextPointX;
					}
					else {
						gradient1 = weightedDeltaY/(nextPointX-minX);
						next1 = nextPointX;
					}
					gradient = gradient0 + gradient1;
					if (next0 <= next1) {
						maxX = next0;
						nextIndex = 0;
					}
					else {
						maxX = next1;
						nextIndex = 1;
					}
					maxY += gradient*(maxX - minX);
				}
				else {
					gradient = Double.POSITIVE_INFINITY;
					maxY += weightedDeltaY;
				}
				return true;
			}
			else {
				if (nextIndex==0) {
					gradient0 = 0.0;
					next0 = Double.POSITIVE_INFINITY;
					maxX = next1;
					nextIndex = 1;
				}
				else {
					gradient1 = 0.0;
					next1 = Double.POSITIVE_INFINITY;
					maxX = next0;
					nextIndex = 0;
				}
				gradient = Double.POSITIVE_INFINITY;
				return (maxX != Double.POSITIVE_INFINITY);
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
				final double weightedDeltaY = histogram.getFrequency(pos-1);
				if (minX < nextPointX) {
					gradient = gradientSum.update(histogramIdx, weightedDeltaY/(nextPointX-minX));
					heap.update(nextPointX);
					maxX = heap.getMinValue();
					maxY += gradient*(maxX - minX);
				}
				else {
					gradient = Double.POSITIVE_INFINITY;
					maxY += weightedDeltaY;
				}
				return true;
			}
			else {
				heap.update(Double.POSITIVE_INFINITY);
				gradient = gradientSum.update(histogramIdx, 0.0);
				maxX = heap.getMinValue();
				return (maxX != Double.POSITIVE_INFINITY);
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
				if (cumulativeFrequencies[valueCounter] != iterator.getMinimumCumulativeFrequency()) {
					while(valueCounter < cumulativeFrequencies.length && cumulativeFrequencies[valueCounter] < iterator.getMaximumCumulativeFrequency()) {
						values[valueCounter] = interpolate(iterator.getMinimumCumulativeFrequency(), iterator.getMinimumValue(), iterator.getMaximumCumulativeFrequency(), iterator.getMaximumValue(), cumulativeFrequencies[valueCounter]);
						valueCounter += 1;
					}
					if (valueCounter < cumulativeFrequencies.length && cumulativeFrequencies[valueCounter] == iterator.getMaximumCumulativeFrequency()) {
						values[valueCounter] = iterator.getMaximumValue();
					}
				}
				else {
					values[valueCounter]+=iterator.getMinimumValue();
					values[valueCounter]*=0.5;
					valueCounter += 1;
				}
			}
		}
		
		while(valueCounter < cumulativeFrequencies.length) {
			values[valueCounter] = iterator.getMinimumValue();
			valueCounter += 1;
		}
	
		return values;
	}
}
