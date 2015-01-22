package org.apache.commons.math3.stat.descriptive.rank;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertArrayEquals;

import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.rank.TDigestQuantile.PartitionStrategy1;
import org.apache.commons.math3.stat.descriptive.rank.TDigestQuantile.PartitionStrategy2;
import org.apache.commons.math3.stat.descriptive.rank.TDigestQuantile.PartitionStrategy3;
import org.apache.commons.math3.stat.descriptive.rank.TDigestQuantile.PartitionStrategy4;
import org.apache.commons.math3.stat.descriptive.rank.TDigestQuantile.PartitionStrategy;
import org.apache.commons.math3.stat.descriptive.rank.TDigestQuantile.SimpleBufferStrategy;
import org.apache.commons.math3.stat.descriptive.rank.TDigestQuantile.SortedCentroids;
import org.junit.Ignore;
import org.junit.Test;

public class TDigestQuantileTest {
	
	private static final PartitionStrategy NEVER_MERGE_PARTITION_STRATEGY = new PartitionStrategy() {
		public int partition(SortedCentroids sortedCentroids, int[] partitionIndices) {
			int size = sortedCentroids.size();
			for (int i = 0; i < size; ++i) {
				partitionIndices[i] = i+1;
			}
			return size;
		}
	};
	
	private static final PartitionStrategy ALWAYS_MERGE_PARTITION_STRATEGY = new PartitionStrategy() {
		public int partition(SortedCentroids sortedCentroids, int[] partitionIndices) {
			int size = sortedCentroids.size();
			if (size > 0) {
				partitionIndices[0] = size;
				return 1;
			}
			else {
				return 0;
			}
		}
	};
	
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
		// System.out.println(Arrays.toString(Arrays.copyOf(partitionIndices, numPartitions)));
		assertArrayEquals(expectedPartitionIndices, Arrays.copyOf(partitionIndices, numPartitions));		
	}
	
	@Test
	public void testPartitionStrategies() {

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
		
		testPartition(centroids, new PartitionStrategy1(0.05), new int[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20});
		testPartition(centroids, new PartitionStrategy1(0.2), new int[] {1, 2, 3, 5, 8, 11, 14, 16, 17, 18, 19, 20});
		testPartition(centroids, new PartitionStrategy1(1.0), new int[] {1, 14, 19, 20});
		testPartition(centroids, new PartitionStrategy1(5.0), new int[] {1, 19, 20});
		testPartition(centroids, new PartitionStrategy1(1000.0), new int[] {1, 19, 20});

		testPartition(centroids, new PartitionStrategy2(0.05), new int[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20});
		testPartition(centroids, new PartitionStrategy2(0.2), new int[] {1, 2, 3, 5, 8, 11, 14, 16, 17, 18, 19, 20});
		testPartition(centroids, new PartitionStrategy2(1.0), new int[] {1, 14, 19, 20});
		testPartition(centroids, new PartitionStrategy2(5.0), new int[] {1, 19, 20});
		testPartition(centroids, new PartitionStrategy2(1000.0), new int[] {1, 19, 20});

		testPartition(centroids, new PartitionStrategy3(0.05), new int[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20});
		testPartition(centroids, new PartitionStrategy3(0.2), new int[] {1, 2, 4, 7, 10, 13, 16, 18, 19, 20});
		testPartition(centroids, new PartitionStrategy3(1.0), new int[] {20});
		testPartition(centroids, new PartitionStrategy3(5.0), new int[] {20});
		testPartition(centroids, new PartitionStrategy3(1000.0), new int[] {20});

		testPartition(centroids, new PartitionStrategy4(0.05), new int[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20});
		testPartition(centroids, new PartitionStrategy4(0.2), new int[] {1, 3, 6, 9, 12, 15, 17, 19, 20});
		testPartition(centroids, new PartitionStrategy4(1.0), new int[] {20});
		testPartition(centroids, new PartitionStrategy4(5.0), new int[] {20});
		testPartition(centroids, new PartitionStrategy4(1000.0), new int[] {20});
	}
	
	@Ignore
	@Test
	public void testPerformanceRandom() {
		
		final int M = 100;
		final int N = 1000000;
		
		final double values[] = new double[N];
		Random random = new Random(0);
		for (int i = 0; i < N; ++i) {
			values[i] = random.nextDouble();
		}
		
		long  start = System.currentTimeMillis();
		for (int m = 0; m < M; ++m) {			
			final TDigestQuantile quantile =  new TDigestQuantile(new PartitionStrategy4(1e-3), new SimpleBufferStrategy(2.0, 10));
			for (int i = 0; i < N; ++i) {
				quantile.add(values[i]);
			}
		}
		long end = System.currentTimeMillis(); 		
		
		System.out.println("Avg time add operation unsorted data = " + ((end - start)*1e6)/(N*M) + "ns.");
		
	}

	@Ignore
	@Test
	public void testPerformanceSequential() {
		
		final int M = 100;
		final int N = 1000000;
		
		final double values[] = new double[N];
		Random random = new Random(0);
		for (int i = 0; i < N; ++i) {
			values[i] = random.nextDouble();
		}
		Arrays.sort(values);
		
		long  start = System.currentTimeMillis();
		for (int m = 0; m < M; ++m) {			
			final TDigestQuantile quantile =  new TDigestQuantile(new PartitionStrategy4(1e-3), new SimpleBufferStrategy(2.0, 10));
			for (int i = 0; i < N; ++i) {
				quantile.add(values[i]);
			}
		}
		long end = System.currentTimeMillis(); 		
		
		System.out.println("Avg time add operation sorted data = " + ((end - start)*1e6)/(N*M) + "ns.");
		
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
	
	@Test
	public void testGetQuantileNotMerged() {
		
		// TODO improve test
		
		TDigestQuantile tDigestQuantile = new TDigestQuantile(NEVER_MERGE_PARTITION_STRATEGY, new SimpleBufferStrategy(0., 1));
		
		double min = 3.4;
		double max = 17.8;
		int N_ADD = 100;
		
		for (int i = 0; i < N_ADD; ++i) {
			tDigestQuantile.add(min+i*(max-min)/(N_ADD-1));
		}
	
		double eps = 1e-8;
		
		assertEquals(min, tDigestQuantile.getQuantile(0.0), eps);
		assertEquals((min+max)*0.5, tDigestQuantile.getQuantile(0.5), eps);
		assertEquals(max, tDigestQuantile.getQuantile(1.0), eps);
		
		int N_TEST = 23452;
		
		for (int i = 0; i < N_TEST; ++i) {
			assertEquals(min + i*(max-min)/(N_TEST-1), tDigestQuantile.getQuantile(0.5/N_ADD + ((double)i/(N_TEST-1))*((double)(N_ADD-1))/N_ADD), eps);
		}
		
		for (int i = 0; i < N_TEST; ++i) {
			assertEquals(min, tDigestQuantile.getQuantile((double)i/(N_TEST-1)*0.5/N_ADD), eps);
		}
		for (int i = 0; i < N_TEST; ++i) {
			assertEquals(max, tDigestQuantile.getQuantile(1.-0.5/N_ADD + (double)i/(N_TEST-1)*0.5/N_ADD), eps);
		}
	}
	
	@Test
	public void testGetQuantileAllMerged() {
		
		// TODO improve test
		
		TDigestQuantile tDigestQuantile = new TDigestQuantile(ALWAYS_MERGE_PARTITION_STRATEGY, new SimpleBufferStrategy(0., 1));
		
		double min = 3.4;
		double max = 17.8;
		int N_ADD = 100;
		
		for (int i = 0; i < N_ADD; ++i) {
			tDigestQuantile.add(min+i*(max-min)/(N_ADD-1));
		}
	
		double eps = 1e-8;
		
		assertEquals(min, tDigestQuantile.getQuantile(0.0), eps);
		assertEquals((min+max)*0.5, tDigestQuantile.getQuantile(0.5), eps);
		assertEquals(max, tDigestQuantile.getQuantile(1.0), eps);
		
		int N_TEST = 23452;
		
		for (int i = 0; i < N_TEST; ++i) {
			assertEquals(min + i*(max-min)/(N_TEST-1), tDigestQuantile.getQuantile(0.5/N_ADD + ((double)i/(N_TEST-1))*((double)(N_ADD-1))/N_ADD), eps);
		}
		
		for (int i = 0; i < N_TEST; ++i) {
			assertEquals(min, tDigestQuantile.getQuantile((double)i/(N_TEST-1)*0.5/N_ADD), eps);
		}
		for (int i = 0; i < N_TEST; ++i) {
			assertEquals(max, tDigestQuantile.getQuantile(1.-0.5/N_ADD + (double)i/(N_TEST-1)*0.5/N_ADD), eps);
		}
	}
	
	@Test
	public void testAddManyValues() {
		
		Random random = new Random(0);
		int numValues = 100000;
		
		double[] values = new double[numValues];
		for (int i = 0; i < numValues; ++i) {
			values[i] = random.nextDouble();
		}
		
		final double[] pValues = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
		
		TDigestQuantile quantileEstimator = new TDigestQuantile(new PartitionStrategy1(1e-2), new SimpleBufferStrategy(0., 100));
		
		for (double value : values) {
			quantileEstimator.add(value);
		}
		
		for (double pValue : pValues) {
			assertEquals(pValue, quantileEstimator.getQuantile(pValue), 0.01);
		}
	}
	
	
}
