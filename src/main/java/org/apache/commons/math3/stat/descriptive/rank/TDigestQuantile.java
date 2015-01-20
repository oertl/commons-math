package org.apache.commons.math3.stat.descriptive.rank;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NonMonotonicSequenceException;
import org.apache.commons.math3.exception.NotPositiveException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.stat.descriptive.rank.MinHeapUtils.MinDoubleIntHeap;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.MathUtils;
import org.apache.commons.math3.util.MathArrays.OrderDirection;


// TODO introduce abstract class for TDigestQuantile and IQAgentQuantile
public class TDigestQuantile {
	
	/**
	 * A partition strategy.
	 * <p>
	 * Implementations of this interface must be immutable.
	 */
	public interface PartitionStrategy {
		
		/**
		 * Calculates a partition for given sorted centroids.
		 * 
		 * @param sortedCentroids
		 * @param partitionIndices an array of size {@code centoids.size()} which is used to store the calculated partition indices
		 * @return the number of partitions
		 */
		int partition(SortedCentroids sortedCentroids, int[] partitionIndices);
	}
	
	
	public static abstract class AbstractPartitionStrategy implements PartitionStrategy {

		protected interface Partitioner {
			/**
			 * Returns {@code true} if the partition criterion is violated, meaning that
			 * the range [startIdx, endIdx[ need to be  split. 
			 * 
			 * @param startIdx
			 * @param endIdx
			 * @return
			 */
			boolean partitionCriterionViolated(int startIdx, int endIdx);
		}
		
		/**
		 * Creates a {@link Partitioner} for which the corresponding
		 * {@link Partitioner#partitionCriterionViolated(int startIdx, int endIdx)}
		 * must be defined for all {@code 0 <= startIdx < endIdx < sortedCentroids.size()}
		 * 
		 * @param sortedCentroids
		 * @return
		 */
		protected abstract Partitioner createPartitioner(SortedCentroids sortedCentroids);
		
		
		public final int partition(SortedCentroids sortedCentroids, int[] partitionIndices) {

			final Partitioner partitioner = createPartitioner(sortedCentroids);
			final int size = sortedCentroids.size();
			
			// partition
	        int partitionCounter = 0;
	        int startIdx = 0;
	        for (int endIdx = 1; endIdx < size; ++endIdx) {
	        	if (partitioner.partitionCriterionViolated(startIdx, endIdx)) {
	        		partitionIndices[partitionCounter] = endIdx;
	        		partitionCounter += 1;
	        		startIdx = endIdx;
	        	}
	        }
	        partitionIndices[partitionCounter] = size;
        	partitionCounter += 1;
	        return partitionCounter;
		}
	}
	
	/**
	 * Partition strategy that limits the integral/mean over a centroid of 1/(4*q*(1-q)).
	 */
	public final static class PartitionStrategy1 extends AbstractPartitionStrategy {

		private final double limit;
		
		public PartitionStrategy1(double delta) {
			if (delta < 0.) {
				throw new NotPositiveException(delta);
			}
			this.limit = FastMath.expm1(4.*delta);
		}
		
		@Override
		protected Partitioner createPartitioner(final SortedCentroids sortedCentroids) {
			return new Partitioner() {
				
				private final long w = sortedCentroids.getTotalWeight();
				private final double lim = Math.min(limit/w, Double.MAX_VALUE);
				
				public boolean partitionCriterionViolated(int startIdx, int endIdx) {					
					long q1 = (startIdx>0)?sortedCentroids.getAccumulatedWeight(startIdx-1):0L;
					long q2 = sortedCentroids.getAccumulatedWeight(endIdx);
					return (q2-q1) > lim*q1*(w-q2);
				}
			};
		}
	}
	
	/**
	 * Partition strategy that limits the integral/mean over a centroid of 1/(4*q*(1-q)).
	 * This strategy demonstrates the partition problem as partition problem of sums.
	 */
	public final static class PartitionStrategy2 implements PartitionStrategy {
		
		private final double delta;
			
		public PartitionStrategy2(double delta) {
			
			if (delta < 0.) {
				throw new NotPositiveException(delta);
			}
			
			this.delta = delta;
		}

		public int partition(SortedCentroids sortedCentroids, int[] partitionIndices) {
			
			MathUtils.checkNotNull(sortedCentroids);
			MathUtils.checkNotNull(partitionIndices);

			final int size = sortedCentroids.size();

			if (size != partitionIndices.length) {
				throw new DimensionMismatchException(partitionIndices.length, size);
			}
			
	        if (size == 0) {
	        	return 0;
	        }

			final long totalWeight = sortedCentroids.getTotalWeight();
			
			// calculate fill levels
			final double[] fillLevels = new double[size];
	        {
		        double lastValue = Double.NEGATIVE_INFINITY;
		        for(int i = 0; i < size; ++i) {
		        	final long accumulatedWeight = sortedCentroids.getAccumulatedWeight(i);
		        	final double nextValue = 0.25d*Math.log(((double)accumulatedWeight)/(totalWeight-accumulatedWeight));	
		        	fillLevels[i] = nextValue - lastValue;
		        	lastValue = nextValue;
		        }
	        }

	        // partition
	        int partitionCounter = 0;
	        double currentFillLevel = fillLevels[0];
	        for (int i = 1; i < size; ++i) {
	        	final double newFillLevel = fillLevels[i] + currentFillLevel;
	        	if (newFillLevel > delta) {
	        		partitionIndices[partitionCounter] = i;
	        		partitionCounter += 1;
	        		currentFillLevel = fillLevels[i];	        				        		
	        	}
	        	else {
	        		currentFillLevel = newFillLevel;
	        	}
	        }
	        partitionIndices[partitionCounter] = size;
        	partitionCounter += 1;
	        return partitionCounter;
		}	
	}
	
	/**
	 * Partition strategy that limits 1/(4*q*(1-q)) where q is the mean of a centroid.
	 * This strategy is similar to the criterion used in the original approach by Ted Dunning.
	 */
	public final static class PartitionStrategy3 extends AbstractPartitionStrategy {

		private final double delta;
		
		public PartitionStrategy3(double delta) {
			if (delta < 0.) {
				throw new NotPositiveException(delta);
			}
			this.delta = delta;
		}
		
		@Override
		protected Partitioner createPartitioner(final SortedCentroids sortedCentroids) {
			return new Partitioner() {
				
				private final long w2;
				private final double lim;
				{
					long w = sortedCentroids.getTotalWeight();
					w2 = 2*w;
					lim = Math.min(delta/w, Double.MAX_VALUE);
				}
				
				
				public boolean partitionCriterionViolated(int startIdx, int endIdx) {					
					long q1 = (startIdx>0)?sortedCentroids.getAccumulatedWeight(startIdx-1):0L;
					long q2 = sortedCentroids.getAccumulatedWeight(endIdx);
					long q = q1 +q2;
					return (q2-q1) > lim*q*(w2-q);
				}
			};
		}
	}
	
	/**
	 * Partition strategy that limits the integral/mean over a centroid of 1/sqrt(4*q*(1-q)).
	 * TODO Still needs to be optimized!
	 */
	public final static class PartitionStrategy4 extends AbstractPartitionStrategy {

		private final double delta;
		
		public PartitionStrategy4(double delta) {
			if (delta < 0.) {
				throw new NotPositiveException(delta);
			}
			this.delta = delta;
		}
		
		@Override
		protected Partitioner createPartitioner(final SortedCentroids sortedCentroids) {
			return new Partitioner() {
				
				private final double w = 2./sortedCentroids.getTotalWeight();
				private final double lim = 2.*delta;
				
				public boolean partitionCriterionViolated(int startIdx, int endIdx) {
					long q1 = (startIdx>0)?sortedCentroids.getAccumulatedWeight(startIdx-1):0L;
					long q2 = sortedCentroids.getAccumulatedWeight(endIdx);
					return FastMath.asin(w*q2)-FastMath.asin(w*q1) > lim; // TODO use addition law for asin to avoid asin evaluations
				}
			};
		}
	}

	private PartitionStrategy partitionStrategy;
	
	private long[] accumulatedCentroidWeights;
	private double[] centroidMeans;	
	private final double[] buffer;
	private int bufferCounter;
	private double minimum;
	private double maximum;
	
	
	// TODO introduce buffer strategy, choose buffer size relative to number of centroids
	public TDigestQuantile(PartitionStrategy partitionStrategy, int bufferSize) {
		
		assert bufferSize > 0;
		
		this.bufferCounter = 0;
		this.buffer = new double[bufferSize];
		this.centroidMeans = new double[0];
		this.accumulatedCentroidWeights = new long[0];
		this.partitionStrategy = partitionStrategy;
		this.minimum = Double.POSITIVE_INFINITY;
		this.maximum = Double.NEGATIVE_INFINITY;
	}

	public void add(double value) {
		MathUtils.checkFinite(value);
		buffer[bufferCounter] = value;
		bufferCounter += 1;
		if (bufferCounter == buffer.length) {
			collapse();
		}
	}
		
	public interface SortedCentroids {
		
		int size();
		
		long getWeight(int index);
		
		long getAccumulatedWeight(int index);
		
		double getMean(int index);
		
		long getTotalWeight();
		
	}
	
	static SortedCentroids asSortedCentroids(final double[] means, final long[] accumulatedWeights) {
		
		MathUtils.checkNotNull(means);
		MathUtils.checkNotNull(accumulatedWeights);
		
		int size = means.length;
		
		if (size != accumulatedWeights.length) {
			throw new DimensionMismatchException(accumulatedWeights.length, size);
		}
		
		if (size > 0) {	
			MathUtils.checkFinite(means[0]);
			MathArrays.checkOrder(means, OrderDirection.INCREASING, false);	// -0.0d and 0.0d are treated as equal		
			MathUtils.checkFinite(means[size-1]);
			
			if (accumulatedWeights[0] < 1) {
				throw new NotStrictlyPositiveException(accumulatedWeights[0]);
			}

			for (int i = 1; i < size; ++i) {
				if (accumulatedWeights[i-1] >= accumulatedWeights[i]) {
					throw new NonMonotonicSequenceException(accumulatedWeights[i-1], accumulatedWeights[i], i, OrderDirection.INCREASING, true);
				}
			}	
		}
		
		return asSortedCentroidsUnsafe(means, accumulatedWeights);
	}
	
	private static SortedCentroids asSortedCentroidsUnsafe(final double[] means, final long[] accumulatedWeights) {
				
		final int size = means.length;
		
		return new SortedCentroids() {
			
			public int size() {
				return size;
			}
			
			public long getWeight(int index) {
				return (index>0)?(accumulatedWeights[index]-accumulatedWeights[index-1]):accumulatedWeights[0];
			}
			
			public double getMean(int index) {
				return means[index];
			}

			public long getAccumulatedWeight(int index) {
				return accumulatedWeights[index];
			}

			public long getTotalWeight() {
				return (accumulatedWeights.length>0)?accumulatedWeights[accumulatedWeights.length-1]:0;
			}
		};
		
	}
	
	static SortedCentroids asSortedCentroids(final double[] values, final int size) {
		
		MathUtils.checkNotNull(values);
		
		if (size < 0) {
			throw new NotPositiveException(size);
		}
		
		if (size > values.length) {
			throw new NumberIsTooLargeException(size, values.length, true);
		}

		if (size > 0) {	
			MathUtils.checkFinite(values[0]);
			for (int i = 1; i < size; ++i) {
				if (values[i-1] > values[i]) {
					throw new NonMonotonicSequenceException(values[i-1], values[i], i, OrderDirection.INCREASING, false);
				}
			}	
			MathUtils.checkFinite(values[size-1]);
		}
		
		return asSortedCentroidsUnsafe(values, size);
	}
		
	private static SortedCentroids asSortedCentroidsUnsafe(final double[] values, final int size) {
		
		return new SortedCentroids() {
			
			public int size() {
				return size;
			}
			
			public long getWeight(int index) {
				return 1;
			}
			
			public double getMean(int index) {
				return values[index];
			}

			public long getAccumulatedWeight(int index) {
				return index+1;
			}

			public long getTotalWeight() {
				return size;
			}
		};
		
	}
	
	SortedCentroids merge(List<? extends SortedCentroids> sortedCentroids) {
		
		final int numSortedCentroids = sortedCentroids.size();
		double[] firstMeans = new double[numSortedCentroids];
		
		int finalSize = 0;
		for(int i = 0; i < numSortedCentroids; ++i) {
			final SortedCentroids centroids = sortedCentroids.get(i);
			firstMeans[i] = (centroids.size()>0)?centroids.getMean(0):Double.POSITIVE_INFINITY;
			finalSize += centroids.size();
		}
		final int[] centroidCounters = new int[numSortedCentroids];
		final MinDoubleIntHeap heap = MinHeapUtils.createMinDoubleIntHeap(firstMeans);
		
		final double[] mergedMeans = new double[finalSize];
		final long[] mergedAccumulatedWeights = new long[finalSize];
		
		int mergedCounter = 0;
		long mergedAccumulatedWeight = 0;
		
		while(heap.getMinValue()!=Double.POSITIVE_INFINITY) {
			int index = heap.getMinIndex();
			int centroidCounter = centroidCounters[index];
			SortedCentroids centroids = sortedCentroids.get(index);
			mergedAccumulatedWeight += centroids.getWeight(centroidCounter);
			mergedMeans[mergedCounter] = heap.getMinValue();
			mergedAccumulatedWeights[mergedCounter] = mergedAccumulatedWeight;
			mergedCounter += 1;
			centroidCounter += 1;
			heap.update((centroidCounter < centroids.size())?centroids.getMean(centroidCounter):Double.POSITIVE_INFINITY);
			centroidCounters[index] = centroidCounter;
		}
		return asSortedCentroidsUnsafe(mergedMeans, mergedAccumulatedWeights);
	}
		
	private void collapse() {
		
		if (bufferCounter == 0) {
			return;
		}
		
		Arrays.sort(buffer, 0, bufferCounter);
		
		// update minimum
		if (buffer[0] < minimum) {
			minimum = buffer[0];
		}
		
		// update maximum
		if (buffer[buffer.length-1] > maximum) {
			maximum = buffer[buffer.length-1];
		}

		final SortedCentroids mergedSortedCentroids = merge(
				Arrays.asList(
					asSortedCentroidsUnsafe(centroidMeans, accumulatedCentroidWeights),
					asSortedCentroidsUnsafe(buffer, bufferCounter)));
		
        // partition centroids        
        final int[] partitionIndices = new int[mergedSortedCentroids.size()];
        final int numPartitions = partitionStrategy.partition(mergedSortedCentroids, partitionIndices);
        
        // calculate new centroids
        accumulatedCentroidWeights = new long[numPartitions];
        centroidMeans = new double[numPartitions];
        calculateMergedCentroids(mergedSortedCentroids, partitionIndices, numPartitions, centroidMeans, accumulatedCentroidWeights);
        
        // reset buffer counter
		bufferCounter = 0;
	}
	
	static void calculateMergedCentroids(
			final SortedCentroids centroids, 
			final int[] partitionIndices,
			final int numPartitions, 
			final double[] mergedCentroidMeans,
			final long[] mergedAccumulatedCentroidWeights) {
		
		int centroidCounter = 0;
        
        long accumulatedWeight = 0l;
        
        for (int i = 0; i < numPartitions; ++i) {
        	
        	final int nextIndex = partitionIndices[i]; // always greater than 0
        	
        	mergedAccumulatedCentroidWeights[i] = centroids.getAccumulatedWeight(nextIndex-1);
        	
        	final long newWeight = centroids.getAccumulatedWeight(nextIndex-1) - accumulatedWeight;
        	accumulatedWeight = centroids.getAccumulatedWeight(nextIndex-1);
        	
        	final double firstMean = centroids.getMean(centroidCounter);
        	double mean = firstMean;
        	++centroidCounter;
        	while(centroidCounter < nextIndex) {
        		final long centroidWeight = centroids.getAccumulatedWeight(centroidCounter) - centroids.getAccumulatedWeight(centroidCounter-1);
        		mean += (((double)centroidWeight)/newWeight) *  (centroids.getMean(centroidCounter) - firstMean);
        		centroidCounter+=1;
        	}
        	
        	if (mean > centroids.getMean(centroidCounter-1)) {
        		mean = centroids.getMean(centroidCounter-1);
        	}
        	
        	mergedCentroidMeans[i] = mean;
        }
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
	static double interpolate(double x1, double y1, double x2, double y2, double x) { // TODO could be shared with IQAgentQuantile
		double alpha = (x - x1)/(x2 - x1);
		double result = y1 + (y2 - y1)*alpha;
		return (result <= y2)?result:y2;
	}


	// TODO improve method, especially the 4 different cases at its end 
	public double getQuantile(double pValue) {
		
		// TODO improve argument checking
		if (pValue < 0. || pValue > 1) {
			throw new OutOfRangeException(pValue, 0., 1.);
		}

		collapse();

		int numCentroids = centroidMeans.length;
		
		if (numCentroids == 0) {
			return Double.NaN;
		}

		final long totalCount = accumulatedCentroidWeights[accumulatedCentroidWeights.length-1];
	
		final double desiredAccumulatedCentroidWeight = pValue * totalCount;
		
		if (desiredAccumulatedCentroidWeight <= 0.5) {
			return minimum;
		}
		if (desiredAccumulatedCentroidWeight >= totalCount-0.5) {
			return maximum;
		}
		
		
		int minIndex = -1;
		int maxIndex = numCentroids;
		
		while (maxIndex-minIndex > 1) {
			
			int midIndex = (maxIndex+minIndex) >>> 1;
		
			long midWeight = accumulatedCentroidWeights[midIndex];
			if (midWeight > desiredAccumulatedCentroidWeight) {
				maxIndex = midIndex;	
			}
			else {
				minIndex = midIndex;
			}
		}
		
		// assert minIndex == maxIndex+1
		// (minIndex>=0)?accumulatedCentroidWeights[minIndex]:0 <= desiredAccumulatedCentroidWeight < accumulatedCentroidWeights[maxIndex]
		
		final long lowerAccumulatedCentroidWeight = (minIndex>=0)?accumulatedCentroidWeights[minIndex]:0;
		final long upperAccumulatedCentroidWeight = accumulatedCentroidWeights[maxIndex];
		final double centerAccumulatedCentroidWeight = lowerAccumulatedCentroidWeight + (upperAccumulatedCentroidWeight-lowerAccumulatedCentroidWeight)*0.5;
		
		if (desiredAccumulatedCentroidWeight < centerAccumulatedCentroidWeight)  {
			if (minIndex >= 0) {
				// interpolate between centroid[minIndex] and centroid[minIndex+1]
				final long lowerAccumulatedCentroidWeight2 = (minIndex>0)?accumulatedCentroidWeights[minIndex-1]:0;
				
				double upperWeight = centerAccumulatedCentroidWeight;
				double lowerWeight = (lowerAccumulatedCentroidWeight + lowerAccumulatedCentroidWeight2)*0.5;
				
				double upperValue = centroidMeans[minIndex+1];
				double lowerValue = centroidMeans[minIndex];
				
				return interpolate(lowerWeight, lowerValue, upperWeight, upperValue, desiredAccumulatedCentroidWeight);
			}
			else {
				// interpolate between minimum and 1st centroid

				double upperWeight = centerAccumulatedCentroidWeight;
				double lowerWeight = 0.5;
				
				double upperValue = centroidMeans[minIndex+1];
				double lowerValue = minimum;
				
				return interpolate(lowerWeight, lowerValue, upperWeight, upperValue, desiredAccumulatedCentroidWeight);
			}
		} else {
			if (maxIndex+1 < numCentroids) {
				// interpolate between centroid[maxIndex] and centroid [maxIndex+1]
				
				final long upperAccumulatedCentroidWeight2 = accumulatedCentroidWeights[maxIndex+1];
				
				double upperWeight = (upperAccumulatedCentroidWeight + upperAccumulatedCentroidWeight2) * 0.5;
				double lowerWeight = centerAccumulatedCentroidWeight;
				
				double upperValue = centroidMeans[maxIndex+1];
				double lowerValue = centroidMeans[maxIndex];
				
				return interpolate(lowerWeight, lowerValue, upperWeight, upperValue, desiredAccumulatedCentroidWeight);
			}
			else {
				// interpolate between last centroid and maximum

				double upperWeight = accumulatedCentroidWeights[numCentroids-1]-0.5;
				double lowerWeight = centerAccumulatedCentroidWeight;
				double upperValue = maximum;
				double lowerValue = centroidMeans[numCentroids-1];
				
				return interpolate(lowerWeight, lowerValue, upperWeight, upperValue, desiredAccumulatedCentroidWeight);
			}
		}
	}
}
