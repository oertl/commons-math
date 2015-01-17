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
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.MathUtils;
import org.apache.commons.math3.util.MathArrays.OrderDirection;

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
	
	public static class DefaultPartitionStrategy implements PartitionStrategy {
		
		private final double delta;
			
		public DefaultPartitionStrategy(double delta) {
			
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
			MathArrays.checkOrder(means, OrderDirection.INCREASING, false);			
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
		return asSortedCentroids(mergedMeans, mergedAccumulatedWeights);
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

	// TODO
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
		if (numCentroids == 1) {
			return centroidMeans[0];
		}
		
		final long totalCount = accumulatedCentroidWeights[accumulatedCentroidWeights.length-1];
		
		final double index = pValue * (totalCount - 1);
		
		int minIndex = 0;
		int maxIndex = numCentroids;
		
		while (maxIndex-minIndex > 1) {
			
			int midIndex = (maxIndex+minIndex) >>> 1;
		
			long midWeight = accumulatedCentroidWeights[midIndex];
			if (midWeight > index) {
				maxIndex = midIndex;	
			}
			else {
				minIndex = midIndex;
			}
		}
		/*
        final double delta = nextIndex - previousIndex;
        final double previousWeight = (nextIndex - index) / delta;
        final double nextWeight = (index - previousIndex) / delta;
        return previousMean * previousWeight + nextMean * nextWeight;*/
		
		/*
		
		
        final int centroidCount = 
        if (centroidCount == 0) {
            return Double.NaN;
        } else if (centroidCount == 1) {
            return centroids().iterator().next().mean();
        }

		
		*/
		
		// TODO not implemented yet
		return 0.;
		
		
		
		
		
	}
	

}
