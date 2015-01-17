package org.apache.commons.math3.stat.descriptive.rank;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertArrayEquals;

import java.util.Arrays;

import org.apache.commons.math3.stat.descriptive.rank.TDigestQuantile.DefaultPartitionStrategy;
import org.apache.commons.math3.stat.descriptive.rank.TDigestQuantile.PartitionStrategy;
import org.apache.commons.math3.stat.descriptive.rank.TDigestQuantile.SortedCentroids;
import org.junit.Ignore;
import org.junit.Test;

import com.sun.javafx.css.CalculatedValue;

public class TDigestQuantileTest {
	
	@Test
	public void testSortedValuesAsSortedCentroids() {
		
		double[] values = {-6.6,-5.1,-2.3, 3.4, 7.1, 6.9};
		
		SortedCentroids centroids = TDigestQuantile.asSortedCentroids(values, 5);
		
		assertEquals(5, centroids.size());
		
		assertEquals(1, centroids.getAccumulatedWeight(0));
		assertEquals(2, centroids.getAccumulatedWeight(1));
		assertEquals(3, centroids.getAccumulatedWeight(2));
		assertEquals(4, centroids.getAccumulatedWeight(3));
		assertEquals(5, centroids.getAccumulatedWeight(4));
		
		assertEquals(1, centroids.getWeight(0));
		assertEquals(1, centroids.getWeight(1));
		assertEquals(1, centroids.getWeight(2));
		assertEquals(1, centroids.getWeight(3));
		assertEquals(1, centroids.getWeight(4));
		
		assertEquals(5, centroids.getTotalWeight());
		
		assertEquals(-6.6, centroids.getMean(0), 0.0);
		assertEquals(-5.1, centroids.getMean(1), 0.0);
		assertEquals(-2.3, centroids.getMean(2), 0.0);
		assertEquals( 3.4, centroids.getMean(3), 0.0);
		assertEquals( 7.1, centroids.getMean(4), 0.0);
	}
	
	@Test
	public void testArraysAsSortedCentroids() {
		
		double[] means = {-1.1, 2.3, 5.6, 6.3};
		long[] accumulatedWeights = {21, 34, 45, 97};
		
		SortedCentroids centroids = TDigestQuantile.asSortedCentroids(means, accumulatedWeights);
		
		assertEquals(4, centroids.size());
		
		assertEquals(21, centroids.getAccumulatedWeight(0));
		assertEquals(34, centroids.getAccumulatedWeight(1));
		assertEquals(45, centroids.getAccumulatedWeight(2));
		assertEquals(97, centroids.getAccumulatedWeight(3));
		
		assertEquals(21, centroids.getWeight(0));
		assertEquals(13, centroids.getWeight(1));
		assertEquals(11, centroids.getWeight(2));
		assertEquals(52, centroids.getWeight(3));
		
		assertEquals(97, centroids.getTotalWeight());
		
		assertEquals(-1.1, centroids.getMean(0), 0.0);
		assertEquals( 2.3, centroids.getMean(1), 0.0);
		assertEquals( 5.6, centroids.getMean(2), 0.0);
		assertEquals( 6.3, centroids.getMean(3), 0.0);
		
	}
	
	private static final void testPartition(SortedCentroids centroids, PartitionStrategy strategy, int[] expectedPartitionIndices) {
		final int[] partitionIndices = new int[centroids.size()];
		final int numPartitions = strategy.partition(centroids, partitionIndices);		
		assertArrayEquals(expectedPartitionIndices, Arrays.copyOf(partitionIndices, numPartitions));		
	}
	
	@Test
	public void testDefaultStrategy() {

		final int N = 20;
		
		double means[] = new double[N];
		long weights[] = new long[N];
		long accumulatedWeights[] = new long[N];
		long accumulatedWeight = 0L;
		for (int i = 0; i < N; ++i) {
			means[i] = i;
			weights[i] = 1;
			accumulatedWeight += weights[i];
			accumulatedWeights[i] = accumulatedWeight;
		}
		final SortedCentroids centroids = TDigestQuantile.asSortedCentroids(means, accumulatedWeights);
		
		testPartition(centroids, new DefaultPartitionStrategy(0.05), new int[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20});

		testPartition(centroids, new DefaultPartitionStrategy(0.2), new int[] {1, 2, 3, 5, 8, 11, 14, 16, 17, 18, 19, 20});

		testPartition(centroids, new DefaultPartitionStrategy(1.0), new int[] {1, 14, 19, 20});
		
		testPartition(centroids, new DefaultPartitionStrategy(5.0), new int[] {1, 19, 20});
		
		testPartition(centroids, new DefaultPartitionStrategy(1000.0), new int[] {1, 19, 20});
	}
	
	@Ignore
	@Test
	public void testPerformance() {
		
		final int N = 500000000;
		final TDigestQuantile quantile =  new TDigestQuantile(new DefaultPartitionStrategy(1e-1), 1000);
		
		for (int i = 0; i < N; ++i) {
			quantile.add(i);
		}
	}
	
	@Test
	public void testCalculateMergedCentroids() {
		
		final double values[] = {-2.6, 3.7, 4.5, 6.7, 11.3, 14.7};
		
		final SortedCentroids centroids = TDigestQuantile.asSortedCentroids(values, values.length);
		final int[] partitionIndices = {2, 5, 6};
		final int numPartitions = partitionIndices.length;
		final double[] mergedCentroidMeans = new double[numPartitions];
		final long[] mergedAccumulatedCentroidWeights = new long[numPartitions];
		
		TDigestQuantile.calculateMergedCentroids(centroids, partitionIndices, numPartitions, mergedCentroidMeans, mergedAccumulatedCentroidWeights);
				
		assertArrayEquals(new long[]{2, 5, 6}, mergedAccumulatedCentroidWeights);
		assertArrayEquals(new double[]{(-2.6 + 3.7)/2., (4.5 + 6.7 + 11.3)/3., 14.7}, mergedCentroidMeans, 1e-8);
	}
	
	
}
