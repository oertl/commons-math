package org.apache.commons.math3.stat.descriptive.rank;

import java.util.Arrays;

import org.apache.commons.math3.exception.NotPositiveException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.util.MathUtils;


// TODO introduce abstract class for TDigestQuantile and IQAgentQuantile
public class TDigestQuantile {
	
	private static final int INITIAL_BUFFER_SIZE = 20;
	private static final double RELATIVE_BUFFER_SIZE = 1.;
	
	
	public interface Partitioner {
		
		// TODO document access in increasing manner
		double accumulatedCentroidSizeLimit(double pValue);
	}
	
	/**
	 * A partition strategy.
	 * <p>
	 * Implementations of this interface must be immutable.
	 */
	public interface PartitionStrategy {

		/**
		 * Creates a {@link Partitioner} for which the corresponding
		 * {@link Partitioner#partitionCriterionViolated(int startIdx, int endIdx)}
		 * must be defined for all {@code 0 <= startIdx < endIdx < sortedCentroids.size()}
		 * 
		 * @param sortedCentroids
		 * @return
		 */
		Partitioner createPartitioner();

	}
	
	
	/**
	 * Partition strategy that limits the integral/mean over a centroid of 1/(4*q*(1-q)).
	 */
	public final static class PartitionStrategy1 implements PartitionStrategy {

		private final double compression;
		
		public PartitionStrategy1(double compression) {
			if (compression < 0.) {
				throw new NotPositiveException(compression);
			}
			this.compression = compression;
		}
	
		public Partitioner createPartitioner() {
			return new Partitioner() {
				public double accumulatedCentroidSizeLimit(double pValue) {
					return 0.25d*Math.log(pValue/(1.-pValue))*compression;
				}
			};
		}
	}
	
	public final static class PartitionStrategy2 implements PartitionStrategy {

		private final double compression;
		
		public PartitionStrategy2(double compression) {
			if (compression < 0.) {
				throw new NotPositiveException(compression);
			}
			this.compression = compression;
		}
	
		public Partitioner createPartitioner() {
			return new Partitioner() {
				public double accumulatedCentroidSizeLimit(double pValue) {
					return pValue*compression;
				}
			};
		}
	}

	
	
	/**
	 * Calculates a partition for given sorted centroids.
	 * 
	 * @param sortedCentroids
	 * @param partitionIndices an array of size {@code centoids.size()} which is used to store the calculated partition indices
	 * @return the number of partitions
	 */
	// TOOD javadoc
	static final int partition(long[] accumulatedWeights, double[] means, int size, PartitionStrategy strategy) {
		
		if (size <= 0) {
			return 0;
		}

		final double totalWeight = accumulatedWeights[size-1];
		final Partitioner partitioner = strategy.createPartitioner();
		
		// partition
        int partitionCounter = 0;
        long startAccumulatedWeights = 0;
        double startSizeLimit = partitioner.accumulatedCentroidSizeLimit(0.);
        long prevEndAccumulatedWeight = accumulatedWeights[0];
        double prevEndSizeLimit = partitioner.accumulatedCentroidSizeLimit(prevEndAccumulatedWeight/totalWeight);
        double firstMean = means[0];
        double mean = firstMean;
        for (int endIdx = 1; endIdx < size; ++endIdx) {
        	long endAccumulatedWeights = accumulatedWeights[endIdx];
        	double endMean = means[endIdx];
        	double endSizeLimit = partitioner.accumulatedCentroidSizeLimit(endAccumulatedWeights/totalWeight);
        	
        	if (!(endSizeLimit-startSizeLimit <= 1.)) {
        		accumulatedWeights[partitionCounter] = prevEndAccumulatedWeight;
        		means[partitionCounter] = mean;
        		firstMean = endMean;
        		mean = firstMean;
        		partitionCounter += 1;
        		startAccumulatedWeights = prevEndAccumulatedWeight;
        		startSizeLimit = prevEndSizeLimit;
        	}
        	else {
        		mean += (((double)(endAccumulatedWeights - prevEndAccumulatedWeight))/(endAccumulatedWeights-startAccumulatedWeights)) * (endMean - mean);
        		if (mean > endMean) {
        			mean = endMean;
        		}
        	}
        	prevEndAccumulatedWeight = endAccumulatedWeights;
        	prevEndSizeLimit = endSizeLimit;
        }
        accumulatedWeights[partitionCounter] = prevEndAccumulatedWeight;
        means[partitionCounter] = mean;
    	partitionCounter += 1;
    	
        return partitionCounter;
	}

	
	/**
	 * Partition strategy that limits the integral/mean over a centroid of 1/sqrt(4*q*(1-q)).
	 * TODO Still needs to be optimized!
	 */
	public final static class PartitionStrategy4 implements PartitionStrategy {

		private final double compression;
		
		public PartitionStrategy4(double compression) {
			if (compression < 0.) {
				throw new NotPositiveException(compression);
			}
			this.compression = compression;
		}
		
		public Partitioner createPartitioner() {
			return new Partitioner() {
				public double accumulatedCentroidSizeLimit(double pValue) {
					return 0.5*Math.asin(2.*pValue-1.)*compression;
				}
			};
		}
	}

	private final PartitionStrategy partitionStrategy;
	
	private long[] accumulatedCentroidWeights;
	private double[] centroidMeans;	
	private double[] bufferValues;
	private long[] bufferWeights;
	private int bufferCounter;
	private int centroidCounter;
	private double minimum;
	private double maximum;
	
	
	public TDigestQuantile(PartitionStrategy partitionStrategy) {
		
		this.bufferCounter = 0;
		
		this.centroidMeans = new double[0];
		this.accumulatedCentroidWeights = new long[0];
		this.partitionStrategy = partitionStrategy;
		this.minimum = Double.POSITIVE_INFINITY;
		this.maximum = Double.NEGATIVE_INFINITY;
		int initialBufferSize = INITIAL_BUFFER_SIZE;
		this.bufferValues = new double[initialBufferSize];
		this.bufferWeights = new long[initialBufferSize];
	}

	public void add(double value) {
		add(value, 1);
	}
	
	public void add(double value, long weight) {
		MathUtils.checkFinite(value);
		bufferValues[bufferCounter] = value;
		bufferWeights[bufferCounter] = weight;
		bufferCounter += 1;
		if (bufferCounter == bufferValues.length) {
			collapse();
		}
	}

	void ensureBufferCapacity(int size) {
		int currentSize = bufferValues.length;
		if (currentSize < size) {
			int newSize = size + (size >>> 1);
			bufferValues = new double[newSize];
			bufferWeights = new long[newSize];
		}
	}
		
	static final void merge(
			double[] means1,
			long[] accumulatedWeights1,
			int offset1,
			int size1,
			double[] means2,
			long[] accumulatedWeights2,
			int offset2,
			int size2,
			double[] mergedMeans,
			long[] mergedAccumulatedWeights,
			int mergedOffset) {
		
		int index1 = offset1;
		int index2 = offset2;
		double mean1 = (0 < size1)?means1[index1]:Double.POSITIVE_INFINITY;
		double mean2 = (0 < size2)?means2[index2]:Double.POSITIVE_INFINITY;
		long accumulatedWeight1 = 0;
		long accumulatedWeight2 = 0;
		int finalMergedIndex = size1 + size2 + mergedOffset;
		int finalIndex1 = size1 + offset1;
		int finalIndex2 = size2 + offset2;
		for(int mergedIndex = mergedOffset; mergedIndex < finalMergedIndex; ++mergedIndex) {
			if (mean1 <= mean2) {
				mergedMeans[mergedIndex] = mean1;
				accumulatedWeight1 = accumulatedWeights1[index1];
				index1 += 1;
				mean1 = (index1 < finalIndex1)?means1[index1]:Double.POSITIVE_INFINITY;
			}
			else {
				mergedMeans[mergedIndex] = mean2;
				accumulatedWeight2 = accumulatedWeights2[index2];
				index2 += 1;
				mean2 = (index2 < finalIndex2)?means2[index2]:Double.POSITIVE_INFINITY;
			}
			mergedAccumulatedWeights[mergedIndex] = accumulatedWeight1 + accumulatedWeight2;
		}
	}
	
	private final void ensureCentroidsCapacity(int minimumSize) {
		int currentSize = centroidMeans.length;
		if (currentSize < minimumSize) {
			int newSize = minimumSize + (minimumSize >>> 1); // TODO improve
			accumulatedCentroidWeights = Arrays.copyOf(accumulatedCentroidWeights, newSize);
			centroidMeans = Arrays.copyOf(centroidMeans, newSize);
		}
	}
	
	static final void mergeReverse(
			double[] means1,
			long[] accumulatedWeights1,
			int offset1,
			int size1,
			double[] means2,
			long[] accumulatedWeights2,
			int offset2,
			int size2,
			double[] mergedMeans,
			long[] mergedAccumulatedWeights,
			int mergedOffset) {
		
		int index1 = offset1 + size1 - 1;
		double mean1;
		long accumulatedWeight1;
		if (index1 >= offset1) {
			accumulatedWeight1 = accumulatedWeights1[index1];
			mean1 = means1[index1];					
		}
		else {
			accumulatedWeight1 = 0L;
			mean1 = Double.NEGATIVE_INFINITY;
		}
		
		int index2 = offset2 + size2 - 1;
		double mean2;
		long accumulatedWeight2;
		if (index2 >= offset2) {
			accumulatedWeight2 = accumulatedWeights2[index2];
			mean2 = means2[index2];					
		}
		else {
			accumulatedWeight2 = 0L;
			mean2 = Double.NEGATIVE_INFINITY;
		}
		
		for (int mergedIndex = mergedOffset + size1 + size2 - 1; mergedIndex >= mergedOffset; mergedIndex -= 1) {
			mergedAccumulatedWeights[mergedIndex] = accumulatedWeight1 + accumulatedWeight2;			
			if (mean1 >= mean2) {
				mergedMeans[mergedIndex] = mean1;
				index1 -= 1;
				if (index1 >= offset1) {
					accumulatedWeight1 = accumulatedWeights1[index1];
					mean1 = means1[index1];					
				}
				else {
					accumulatedWeight1 = 0L;
					mean1 = Double.NEGATIVE_INFINITY;
				}
			}
			else {
				mergedMeans[mergedIndex] = mean2;
				index2 -= 1;
				if (index2 >= offset2) {
					accumulatedWeight2 = accumulatedWeights2[index2];
					mean2 = means2[index2];					
				}
				else {
					accumulatedWeight2 = 0L;
					mean2 = Double.NEGATIVE_INFINITY;
				}
			}
		}
	}
	
	private void collapse() {
		
		if (bufferCounter == 0) {
			return;
		}
		
		AbstractQuantileEstimator.sort(bufferValues, bufferWeights, bufferCounter);
		
		// update minimum
		if (bufferValues[0] < minimum) {
			minimum = bufferValues[0];
		}
		
		// update maximum
		if (bufferValues[bufferCounter-1] > maximum) {
			maximum = bufferValues[bufferCounter-1];
		}
		
		// accumulate buffer weights
		for(int i = 1; i < bufferCounter; ++i) {
			bufferWeights[i] += bufferWeights[i-1];
		}

		ensureCentroidsCapacity(bufferCounter + centroidCounter);
		
		mergeReverse(
				centroidMeans, 
				accumulatedCentroidWeights, 
				0, 
				centroidCounter, 
				bufferValues, 
				bufferWeights,
				0,
				bufferCounter, 
				centroidMeans, 
				accumulatedCentroidWeights, 
				0);
		
		centroidCounter += bufferCounter;
		
		centroidCounter = partition(accumulatedCentroidWeights, centroidMeans, centroidCounter, partitionStrategy);
        
        // reset buffer counter
		bufferCounter = 0;
		ensureBufferCapacity(Math.max(INITIAL_BUFFER_SIZE, (int)(centroidCounter*RELATIVE_BUFFER_SIZE)));
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

	public double getQuantile(double pValue) {
		
		// TODO improve argument checking
		if (pValue < 0. || pValue > 1) {
			throw new OutOfRangeException(pValue, 0., 1.);
		}

		collapse();
		
		if (centroidCounter == 0) {
			return Double.NaN;
		}

		final long totalCount = accumulatedCentroidWeights[centroidCounter-1];
	
		final double desiredAccumulatedCentroidWeight = pValue * totalCount;
		
		if (desiredAccumulatedCentroidWeight <= 0.5) {
			return minimum;
		}
		if (desiredAccumulatedCentroidWeight >= totalCount-0.5) {
			return maximum;
		}
		
		
		int minIndex = -1;
		int maxIndex = centroidCounter;
		
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

		final long lowerAccumulatedCentroidWeight = (minIndex>=0)?accumulatedCentroidWeights[minIndex]:0;
		final long upperAccumulatedCentroidWeight = accumulatedCentroidWeights[maxIndex];
		final double centerAccumulatedCentroidWeight = lowerAccumulatedCentroidWeight + (upperAccumulatedCentroidWeight-lowerAccumulatedCentroidWeight)*0.5;
		
		final double upperWeight;
		final double lowerWeight;
		final double upperValue;
		final double lowerValue;
		if (desiredAccumulatedCentroidWeight < centerAccumulatedCentroidWeight)  {
			upperWeight = centerAccumulatedCentroidWeight;
			upperValue = centroidMeans[minIndex+1];
			if (minIndex >= 0) {
				// interpolate between centroid[minIndex] and centroid[minIndex+1]
				final long lowerAccumulatedCentroidWeight2 = (minIndex>0)?accumulatedCentroidWeights[minIndex-1]:0;
				lowerWeight = (lowerAccumulatedCentroidWeight + lowerAccumulatedCentroidWeight2)*0.5;
				lowerValue = centroidMeans[minIndex];
			}
			else {
				// interpolate between minimum and 1st centroid
				lowerWeight = 0.5;
				lowerValue = minimum;
			}
		} else {
			lowerWeight = centerAccumulatedCentroidWeight;
			lowerValue = centroidMeans[maxIndex];
			if (maxIndex+1 < centroidCounter) {
				// interpolate between centroid[maxIndex] and centroid [maxIndex+1]
				final long upperAccumulatedCentroidWeight2 = accumulatedCentroidWeights[maxIndex+1];
				upperWeight = (upperAccumulatedCentroidWeight + upperAccumulatedCentroidWeight2) * 0.5;				
				upperValue = centroidMeans[maxIndex+1];
			}
			else {
				// interpolate between last centroid and maximum
				upperWeight = accumulatedCentroidWeights[centroidCounter-1]-0.5;
				upperValue = maximum;
			}
		}
		return interpolate(lowerWeight, lowerValue, upperWeight, upperValue, desiredAccumulatedCentroidWeight);
	}
}
