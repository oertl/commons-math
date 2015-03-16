package org.apache.commons.math4.stat.descriptive.rank;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math4.stat.descriptive.rank.TDigestQuantile.PartitionStrategy;
import org.apache.commons.math4.stat.descriptive.rank.TDigestQuantile.PartitionStrategy1;
import org.apache.commons.math4.stat.descriptive.rank.TDigestQuantile.PartitionStrategy2;
import org.apache.commons.math4.stat.descriptive.rank.TDigestQuantile.Partitioner;
import org.junit.Test;

public class TDigestQuantileTest {
	
	private static final PartitionStrategy NEVER_MERGE_PARTITION_STRATEGY = new PartitionStrategy() {
		public Partitioner createPartitioner() {
			return new Partitioner() {
				public double accumulatedCentroidSizeLimit(final double pValue) {
					return Double.NaN;
				}
			};
		}
	};
	
	private static final PartitionStrategy ALWAYS_MERGE_PARTITION_STRATEGY = new PartitionStrategy() {
		public Partitioner createPartitioner() {
			return new Partitioner() {
				public double accumulatedCentroidSizeLimit(final double pValue) {
					return 0;
				}
			};
		}
	};
	
	@Test
	public void testAlwaysMergePartitionStrategy() {
		final int N = 143;
		
		final long[] weights = new long[N];
		final double[] means = new double[N];
		for (int i = 0; i < N; ++i) {
			weights[i] = i+1;
			means[i] = i+1.;
		}
		
		final int numCentroids = TDigestQuantile.partition(weights, means, N, ALWAYS_MERGE_PARTITION_STRATEGY);
		
		assertEquals(1, numCentroids);
		assertEquals((N+1.)*0.5, means[0], 1e-8);
		assertEquals(N, weights[0]);
		
	}
	
	@Test
	public void testNeverMergePartitionStrategy() {
		final int N = 143;
		
		final long[] weights = new long[N];
		final double[] means = new double[N];
		final long[] weightsExpected = new long[N];
		final double[] meansExpected = new double[N];
		for (int i = 0; i < N; ++i) {
			weights[i] = i+1;
			means[i] = i+1.;
			weightsExpected[i] = i+1;
			meansExpected[i] = i+1.;
		}
		
		final int numCentroids = TDigestQuantile.partition(weights, means, N, NEVER_MERGE_PARTITION_STRATEGY);
		
		assertEquals(N, numCentroids);
		assertArrayEquals(meansExpected, means, 1e-8);
		assertArrayEquals(weightsExpected, weights);
	}
	
	/*private static final void testPartition(ExtendedSortedCentroids centroids, PartitionStrategy strategy, int[] expectedPartitionIndices) {
		final int[] partitionIndices = new int[centroids.size()];
		final int numPartitions = strategy.partition(centroids, partitionIndices);
		// System.out.println(Arrays.toString(Arrays.copyOf(partitionIndices, numPartitions)));
		assertArrayEquals(expectedPartitionIndices, Arrays.copyOf(partitionIndices, numPartitions));		
	}*/
/*	
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
		final ExtendedSortedCentroids centroids = TDigestQuantile.asSortedCentroids(means, accumulatedWeights);
		
		testPartition(centroids, new PartitionStrategy1(0.05), new int[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20});
		testPartition(centroids, new PartitionStrategy1(0.2), new int[] {1, 2, 3, 5, 8, 11, 14, 16, 17, 18, 19, 20});
		testPartition(centroids, new PartitionStrategy1(1.0), new int[] {1, 14, 19, 20});
		testPartition(centroids, new PartitionStrategy1(5.0), new int[] {1, 19, 20});
		testPartition(centroids, new PartitionStrategy1(1000.0), new int[] {1, 19, 20});

//		testPartition(centroids, new PartitionStrategy2(0.05), new int[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20});
//		testPartition(centroids, new PartitionStrategy2(0.2), new int[] {1, 2, 3, 5, 8, 11, 14, 16, 17, 18, 19, 20});
//		testPartition(centroids, new PartitionStrategy2(1.0), new int[] {1, 14, 19, 20});
//		testPartition(centroids, new PartitionStrategy2(5.0), new int[] {1, 19, 20});
//		testPartition(centroids, new PartitionStrategy2(1000.0), new int[] {1, 19, 20});

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
	}*/
	
	/*@Test
	public void testPartitionStrategies2() {

		final int N = 10;
		
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
		final ExtendedSortedCentroids centroids = TDigestQuantile.asSortedCentroids(means, accumulatedWeights);
		
		testPartition(centroids, new PartitionStrategy1(0.6), new int[] {1, 5, 9, 10});

		testPartition(centroids, new PartitionStrategy3(0.6), new int[] {3, 8, 10});
	}*/
	
	@Test
	public void testPerformanceRandom() {
		
		final int M = 100;
		final int N = 1000000;
		
		final double values[] = new double[N];
		final Random random = new Random(0);
		for (int i = 0; i < N; ++i) {
			values[i] = random.nextDouble();
		}
		
		final long  start = System.currentTimeMillis();
		for (int m = 0; m < M; ++m) {			
			final TDigestQuantile quantile =  new TDigestQuantile(new PartitionStrategy2(1e3));
			for (int i = 0; i < N; ++i) {
				quantile.add(values[i]);
			}
		}
		final long end = System.currentTimeMillis(); 		
		
		System.out.println("Avg time add operation unsorted data = " + ((end - start)*1e6)/(N*M) + "ns.");
		
	}

	@Test
	public void testPerformanceSequential() {
		
		final int M = 100;
		final int N = 1000000;
		
		final double values[] = new double[N];
		final Random random = new Random(0);
		for (int i = 0; i < N; ++i) {
			values[i] = random.nextDouble();
		}
		Arrays.sort(values);
		
		final long  start = System.currentTimeMillis();
		for (int m = 0; m < M; ++m) {			
			final TDigestQuantile quantile =  new TDigestQuantile(new PartitionStrategy2(1e3));
			for (int i = 0; i < N; ++i) {
				quantile.add(values[i]);
			}
		}
		final long end = System.currentTimeMillis(); 		
		
		System.out.println("Avg time add operation sorted data = " + ((end - start)*1e6)/(N*M) + "ns.");
		
	}

	@Test
	public void testGetQuantileNotMerged() {
		
		// TODO improve test
		
		final TDigestQuantile tDigestQuantile = new TDigestQuantile(NEVER_MERGE_PARTITION_STRATEGY);
		
		final double min = 3.4;
		final double max = 17.8;
		final int N_ADD = 100;
		
		for (int i = 0; i < N_ADD; ++i) {
			tDigestQuantile.add(min+i*(max-min)/(N_ADD-1));
		}
	
		final double eps = 1e-8;
		
		assertEquals(min, tDigestQuantile.getQuantile(0.0), eps);
		assertEquals((min+max)*0.5, tDigestQuantile.getQuantile(0.5), eps);
		assertEquals(max, tDigestQuantile.getQuantile(1.0), eps);
		
		final int N_TEST = 23452;
		
		for (int i = 0; i < N_TEST; ++i) {
			assertEquals(min + i*(max-min)/(N_TEST-1), tDigestQuantile.getQuantile(0.5/N_ADD + ((double)i/(N_TEST-1))*(N_ADD-1)/N_ADD), eps);
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
		
		final TDigestQuantile tDigestQuantile = new TDigestQuantile(ALWAYS_MERGE_PARTITION_STRATEGY);
		
		final double min = 3.4;
		final double max = 17.8;
		final int N_ADD = 100;
		
		for (int i = 0; i < N_ADD; ++i) {
			tDigestQuantile.add(min+i*(max-min)/(N_ADD-1));
		}
	
		final double eps = 1e-8;
		
		assertEquals(min, tDigestQuantile.getQuantile(0.0), eps);
		assertEquals((min+max)*0.5, tDigestQuantile.getQuantile(0.5), eps);
		assertEquals(max, tDigestQuantile.getQuantile(1.0), eps);
		
		final int N_TEST = 23452;
		
		for (int i = 0; i < N_TEST; ++i) {
			assertEquals(min + i*(max-min)/(N_TEST-1), tDigestQuantile.getQuantile(0.5/N_ADD + ((double)i/(N_TEST-1))*(N_ADD-1)/N_ADD), eps);
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
		
		final Random random = new Random(0);
		final int numValues = 100000;
		
		final double[] values = new double[numValues];
		for (int i = 0; i < numValues; ++i) {
			values[i] = random.nextDouble();
		}
		
		final double[] pValues = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
		
		final TDigestQuantile quantileEstimator = new TDigestQuantile(new PartitionStrategy1(1e2));
		
		for (final double value : values) {
			quantileEstimator.add(value);
		}
		
		for (final double pValue : pValues) {
			assertEquals(pValue, quantileEstimator.getQuantile(pValue), 0.01);
		}
	}
	
	@Test
	public void testMerge() {
		
		final long[] accumulatedWeights1 = {11, 7, 3, 2, 3, 5, 8, 12};
		final double[] means1 = {1.0, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.1};
		final int offset1 = 3;
		final int size1 = 4;
		final double[] means2 = {1.1, 1.65, 1.75, 3.5};
		final long[] accumulatedWeights2 = {89, 7, 10, 16};
		final int size2 = 2;
		final int offset2 = 1;
		final double[] mergedMeans = {-1., -2., -3., -4., -5., -6., -7., -8., -9., -10.};
		final long[] mergedAccumulatedWeights = {-1, -2, -3, -4, -5, -6, -7, -8, -9, -10};
		final int mergedOffset = 1;
		
		final double[] expectedMeans = {-1., 1.6 , 1.65, 1.7, 1.75, 1.8, 1.9, -8., -9., -10.};
		final long[] expectedAccumulatedWeights = {-1, 2, 9, 10, 13, 15, 18, -8, -9, -10};
		
		TDigestQuantile.merge(means1, accumulatedWeights1, offset1, size1, means2, accumulatedWeights2, offset2, size2, mergedMeans, mergedAccumulatedWeights, mergedOffset);
		
		assertArrayEquals(expectedMeans, mergedMeans, 0.0d);
		assertArrayEquals(expectedAccumulatedWeights, mergedAccumulatedWeights);
	}
	
	@Test
	public void testMergeReverse() {
		
		final long[] accumulatedWeights1 = {11, 7, 3, 2, 3, 5, 8, 12};
		final double[] means1 = {1.0, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.1};
		final int offset1 = 3;
		final int size1 = 4;
		final double[] means2 = {1.1, 1.65, 1.75, 3.5};
		final long[] accumulatedWeights2 = {89, 7, 10, 16};
		final int size2 = 2;
		final int offset2 = 1;
		final double[] mergedMeans = {-1., -2., -3., -4., -5., -6., -7., -8., -9., -10.};
		final long[] mergedAccumulatedWeights = {-1, -2, -3, -4, -5, -6, -7, -8, -9, -10};
		final int mergedOffset = 1;
		
		final double[] expectedMeans = {-1., 1.6 , 1.65, 1.7, 1.75, 1.8, 1.9, -8., -9., -10.};
		final long[] expectedAccumulatedWeights = {-1, 2, 9, 10, 13, 15, 18, -8, -9, -10};
		
		TDigestQuantile.mergeReverse(means1, accumulatedWeights1, offset1, size1, means2, accumulatedWeights2, offset2, size2, mergedMeans, mergedAccumulatedWeights, mergedOffset);
		
		assertArrayEquals(expectedMeans, mergedMeans, 0.0d);
		assertArrayEquals(expectedAccumulatedWeights, mergedAccumulatedWeights);
	}
	
}
