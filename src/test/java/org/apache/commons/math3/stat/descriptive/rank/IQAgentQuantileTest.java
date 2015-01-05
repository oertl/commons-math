package org.apache.commons.math3.stat.descriptive.rank;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertArrayEquals;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.rank.IQAgentQuantile.DynamicSum;
import org.apache.commons.math3.stat.descriptive.rank.IQAgentQuantile.Histogram;
import org.apache.commons.math3.stat.descriptive.rank.IQAgentQuantile.HistogramIterator;
import org.apache.commons.math3.stat.descriptive.rank.IQAgentQuantile.HistogramIterator1;
import org.apache.commons.math3.stat.descriptive.rank.IQAgentQuantile.MinDoubleIntHeap;
import org.junit.Ignore;
import org.junit.Test;

public class IQAgentQuantileTest {
	
	@Ignore
	@Test
	public void testPerformance() {
		
		double[] pValues = {0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 1.0};
		
		final int N = 500000000;
		final IQAgentQuantile quantile =  new IQAgentQuantile(pValues, 30);
		
		for (int i = 0; i < N; ++i) {
			quantile.add(i);
		}
	}
	
	@Test
	public void testDynamicSum() {
		DynamicSum dynamicSum = new DynamicSum(7);
		assertEquals(1d, dynamicSum.update(0, 1),0);
		assertEquals(3d, dynamicSum.update(1, 2),0);
		assertEquals(6d, dynamicSum.update(2, 3),0);
		assertEquals(10d, dynamicSum.update(3, 4),0);
		assertEquals(15d, dynamicSum.update(4, 5),0);
		assertEquals(21d, dynamicSum.update(5, 6),0);
		assertEquals(28d, dynamicSum.update(6, 7), 0);
		assertEquals(33d, dynamicSum.update(0, 6), 0);
		
	}
	
	@Test
	public void testMinDoubleIntHeap1() {
		
		MinDoubleIntHeap heap = new MinDoubleIntHeap(new double[]{2., 1., 3., 5., 4.});		
		assertEquals(1, heap.getMinIndex());
		assertEquals(1., heap.getMinValue(), 0.);
		
		heap.update(6.);
		assertEquals(0, heap.getMinIndex());
		assertEquals(2., heap.getMinValue(), 0.);

		heap.update(7.);
		assertEquals(2, heap.getMinIndex());
		assertEquals(3., heap.getMinValue(), 0.);

		heap.update(8.);
		assertEquals(4, heap.getMinIndex());
		assertEquals(4., heap.getMinValue(), 0.);
		
		heap.update(9.);
		assertEquals(3, heap.getMinIndex());
		assertEquals(5., heap.getMinValue(), 0.);
		
		heap.update(10.);
		assertEquals(1, heap.getMinIndex());
		assertEquals(6., heap.getMinValue(), 0.);
	}
	
	@Test
	public void testMinDoubleIntHeap2() {
		
		MinDoubleIntHeap heap = new MinDoubleIntHeap(new double[]{2., 1.});		
		assertEquals(1, heap.getMinIndex());
		assertEquals(1., heap.getMinValue(), 0.);
		
		heap.update(3.);
		assertEquals(0, heap.getMinIndex());
		assertEquals(2., heap.getMinValue(), 0.);

		heap.update(4.);
		assertEquals(1, heap.getMinIndex());
		assertEquals(3., heap.getMinValue(), 0.);
	}
	
	@Test
	public void testMinDoubleIntHeap3() {
		
		MinDoubleIntHeap heap = new MinDoubleIntHeap(new double[]{3., 1.});		
		assertEquals(1, heap.getMinIndex());
		assertEquals(1., heap.getMinValue(), 0.);
		
		heap.update(2.);
		assertEquals(1, heap.getMinIndex());
		assertEquals(2., heap.getMinValue(), 0.);
	}
	
	@Test
	public void testMinDoubleIntHeap4() {
		
		MinDoubleIntHeap heap = new MinDoubleIntHeap(new double[]{2., 1., 3., 5., 4.});		
		assertEquals(1, heap.getMinIndex());
		assertEquals(1., heap.getMinValue(), 0.);
		
		heap.update(1.5);
		assertEquals(1, heap.getMinIndex());
		assertEquals(1.5, heap.getMinValue(), 0.);
	}

	
/*	@Test
	public void testDynamicSumPerformance() {
		int N = 100;
		int CYCLES = 5000000;
		
		DynamicSum ds = new DynamicSum(N);
		
		for (int k = 0; k < CYCLES; ++k) {
			
			for (int l = 0; l < N; ++l) {
				
				ds.update(l,  k);
			
			}
		}
		
		//assertEquals((double)N*(CYCLES-1), ds.getSum(), 0.0);
		
		
		
		
	}*/
	
	
	
	@Test
	public void testHistogramIterator1() {
		
		final HistrogramIteratorSupplier iteratorSupplier = getHistogramIterator1Supplier();
		
		double binBoundaries[] = {1., 2., 3.};
		double cumulativeFrequencies[] = {1., 2.};
	
		List<Histogram> histograms = new ArrayList<Histogram>();
		
		histograms.add(IQAgentQuantile.cumulativeFrequenciesAsHistogram(cumulativeFrequencies, binBoundaries, 2.));
		histograms.add(IQAgentQuantile.cumulativeFrequenciesAsHistogram(cumulativeFrequencies, binBoundaries, 2.));
		histograms.add(IQAgentQuantile.cumulativeFrequenciesAsHistogram(cumulativeFrequencies, binBoundaries, 2.));
		
		HistogramIterator iterator = iteratorSupplier.get(histograms);
		
		assertEquals(Double.NEGATIVE_INFINITY, iterator.getMinimumValue(), 0.0);
		assertEquals(1., iterator.getMaximumValue(), 0.0);
		assertEquals(0., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(0., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(0., iterator.getDensity(), 0.0);
		
		assertTrue(iterator.advance());
		
		assertEquals(1., iterator.getMinimumValue(), 0.0);
		assertEquals(1., iterator.getMaximumValue(), 0.0);
		assertEquals(0., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(0., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(1., iterator.getDensity(), 0.0);
		
		assertTrue(iterator.advance());
		
		assertEquals(1., iterator.getMinimumValue(), 0.0);
		assertEquals(1., iterator.getMaximumValue(), 0.0);
		assertEquals(0., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(0., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(2., iterator.getDensity(), 0.0);
		
		assertTrue(iterator.advance());
		
		assertEquals(1., iterator.getMinimumValue(), 0.0);
		assertEquals(2., iterator.getMaximumValue(), 0.0);
		assertEquals(0., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(3., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(3., iterator.getDensity(), 0.0);
		
		assertTrue(iterator.advance());
		
		assertEquals(2., iterator.getMinimumValue(), 0.0);
		assertEquals(2., iterator.getMaximumValue(), 0.0);
		assertEquals(3., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(3., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(3., iterator.getDensity(), 0.0);
		
		assertTrue(iterator.advance());
		
		assertEquals(2., iterator.getMinimumValue(), 0.0);
		assertEquals(2., iterator.getMaximumValue(), 0.0);
		assertEquals(3., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(3., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(3., iterator.getDensity(), 0.0);

		assertTrue(iterator.advance());
		
		assertEquals(2., iterator.getMinimumValue(), 0.0);
		assertEquals(3., iterator.getMaximumValue(), 0.0);
		assertEquals(3., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(6., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(3., iterator.getDensity(), 0.0);

		assertTrue(iterator.advance());
		
		assertEquals(3., iterator.getMinimumValue(), 0.0);
		assertEquals(3., iterator.getMaximumValue(), 0.0);
		assertEquals(6., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(6., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(2., iterator.getDensity(), 0.0);

		assertTrue(iterator.advance());
		
		assertEquals(3., iterator.getMinimumValue(), 0.0);
		assertEquals(3., iterator.getMaximumValue(), 0.0);
		assertEquals(6., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(6., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(1., iterator.getDensity(), 0.0);

		assertFalse(iterator.advance());
		
		assertEquals(3., iterator.getMinimumValue(), 0.0);
		assertEquals(Double.POSITIVE_INFINITY, iterator.getMaximumValue(), 0.0);
		assertEquals(6., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(6., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(0., iterator.getDensity(), 0.0);
	}
	
	@Test
	public void testHistogramIterator2() {

		final HistrogramIteratorSupplier iteratorSupplier = getHistogramIterator1Supplier();

		Collection<Histogram> histograms = new ArrayList<Histogram>();
		
		histograms.add(IQAgentQuantile.sortedValuesAsHistogram(new double[]{2.7}, 3.));
		
		HistogramIterator iterator = iteratorSupplier.get(histograms);
		
		assertEquals(Double.NEGATIVE_INFINITY, iterator.getMinimumValue(), 0.0);
		assertEquals(2.7, iterator.getMaximumValue(), 0.0);
		assertEquals(0., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(0., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(0., iterator.getDensity(), 0.0);
		
		assertTrue(iterator.advance());
		
		assertEquals(2.7, iterator.getMinimumValue(), 0.0);
		assertEquals(2.7, iterator.getMaximumValue(), 0.0);
		assertEquals(0., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(3., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(Double.POSITIVE_INFINITY, iterator.getDensity(), 0.0);
		
		assertFalse(iterator.advance());
		
		assertEquals(2.7, iterator.getMinimumValue(), 0.0);
		assertEquals(Double.POSITIVE_INFINITY, iterator.getMaximumValue(), 0.0);
		assertEquals(3., iterator.getMinimumCumulativeFrequency(), 0.0);
		assertEquals(3., iterator.getMaximumCumulativeFrequency(), 0.0);
		assertEquals(0., iterator.getDensity(), 0.0);
	}
		
	@Test
	public void testSortedValuesAsHistogram() {
		
		final double[] values = new double[]{1., 4., 5.};
		
		final double scale = 7.;
		
		final Histogram histogram = IQAgentQuantile.sortedValuesAsHistogram(values, scale);
		
		assertEquals(5, histogram.getNumberOfBins());
		
		assertEquals(1., histogram.getBoundary(0), 0.);
		assertEquals(1., histogram.getBoundary(1), 0.);
		assertEquals(4., histogram.getBoundary(2), 0.);
		assertEquals(4., histogram.getBoundary(3), 0.);
		assertEquals(5., histogram.getBoundary(4), 0.);
		assertEquals(5., histogram.getBoundary(5), 0.);
		
		assertEquals(scale/values.length, histogram.getFrequency(0), 0.);
		assertEquals(0., histogram.getFrequency(1), 0.);
		assertEquals(scale/values.length, histogram.getFrequency(2), 0.);
		assertEquals(0., histogram.getFrequency(3), 0.);
		assertEquals(scale/values.length, histogram.getFrequency(4), 0.);
	}

	@Test
	public void testCumulativeFrequenciesAsHistogram() {
		
		final double[] cumulativeFrequencies = new double[]{3., 5.};
		final double[] binBoundaries = new double[]{1.,11.,12.};
		
		final double scale = 7.;
		
		final Histogram histogram = IQAgentQuantile.cumulativeFrequenciesAsHistogram(cumulativeFrequencies, binBoundaries, scale);
		
		assertEquals(2, histogram.getNumberOfBins());
		
		assertEquals(1., histogram.getBoundary(0), 0.);
		assertEquals(11., histogram.getBoundary(1), 0.);
		assertEquals(12., histogram.getBoundary(2), 0.);
		
		assertEquals((scale/cumulativeFrequencies[1])*(cumulativeFrequencies[0]), histogram.getFrequency(0), 0.);
		assertEquals((scale/cumulativeFrequencies[1])*(cumulativeFrequencies[1]-cumulativeFrequencies[0]), histogram.getFrequency(1), 0.);
	}

	
	@Test
	public void testInterpolate() {
		{
			final double result = IQAgentQuantile.interpolate(3., 4., 7., 11., 4.);
			assertEquals(7/4.+4., result, 1e-20);				
		}
	}
	
	
	private interface HistrogramIteratorSupplier {
		HistogramIterator get(Collection<? extends Histogram> histograms);
	}
	
	private static final HistrogramIteratorSupplier getHistogramIterator1Supplier() {
		return new HistrogramIteratorSupplier() {
			public HistogramIterator get(Collection<? extends Histogram> histograms) {
				return new HistogramIterator1(histograms);
			}
		};
	}
	
	@Test
	public void testEvaluateSumOfHistogramsForUniformDistribution() {
		
		final double min = -4.7;
		final double max = 5.7;
			
		final Histogram histogram = IQAgentQuantile.cumulativeFrequenciesAsHistogram(new double[]{1.}, new double[]{min, max}, 1.);
		
		int N=100;
		
		final double[] cumulativeFrequencies = new double[N];
		final double[] expectedValues = new double[N];
		for(int i = 0; i < N; ++i) {
			cumulativeFrequencies[i] = i/(N-1.);
			expectedValues[i] = i/(N-1.)*(max-min)+min;
		}
		
		double[] values = IQAgentQuantile.evaluateSumOfHistograms(Collections.singleton(histogram), cumulativeFrequencies);
		
		assertArrayEquals(expectedValues, values, 1e-30);
	}
	
	@Test
	public void testEvaluateSumOfHistograms1() {
		
		final Histogram histogram = IQAgentQuantile.cumulativeFrequenciesAsHistogram(new double[]{1., 2.}, new double[]{1., 2., 3.}, 2.);
		
		final double[] cumulativeFrequencies = {0.5};
		final double[] expectedValues = {1.5};
		double[] values = IQAgentQuantile.evaluateSumOfHistograms(Collections.singleton(histogram), cumulativeFrequencies);
		
		assertArrayEquals(expectedValues, values, 1e-30);
	}
	
	@Test
	public void testEvaluateSumOfHistograms2() {
		
		final Histogram histogram = IQAgentQuantile.cumulativeFrequenciesAsHistogram(new double[]{1., 2.}, new double[]{1., 2., 3.}, 2.);
		
		final double[] cumulativeFrequencies = {-1.};
		final double[] expectedValues = {1.};
		double[] values = IQAgentQuantile.evaluateSumOfHistograms(Collections.singleton(histogram), cumulativeFrequencies);
		
		assertArrayEquals(expectedValues, values, 1e-30);
	}


	
	@Test
	public void testEvaluateSumOfHistogramsForTwoValues() {
		
		final double min = -4.7;
		final double max = 5.7;
			
		final Histogram histogram = IQAgentQuantile.sortedValuesAsHistogram(new double[]{min,max}, 2.);
		
		final double[] cumulativeFrequencies = {-1., 0., 1., 2., 3.};
		final double[] expectedValues = {min, min, (min+max)*0.5, max, max};
		
		double[] values = IQAgentQuantile.evaluateSumOfHistograms(Collections.singleton(histogram), cumulativeFrequencies);
		
		assertArrayEquals(expectedValues, values, 1e-30);
	}
	
	@Test
	public void testEvaluateSumOfHistogramsForTwoValues2() {
		
		final double min = -4.7;
		final double max = 5.7;
			
		final Histogram histogram = IQAgentQuantile.sortedValuesAsHistogram(new double[]{min,max}, 2.);
		
		final double[] cumulativeFrequencies = {-1., 0., 1.};
		final double[] expectedValues = {min, min, (min+max)*0.5};
		
		double[] values = IQAgentQuantile.evaluateSumOfHistograms(Collections.singleton(histogram), cumulativeFrequencies);
		
		assertArrayEquals(expectedValues, values, 1e-30);
	}
	
	

	
	/*
	//@Ignore
	@Test
	public void testHistogramIteratorPerformance() {
		int N_FUNC = 2;
		int N_ITERATIONS = 5000000;
		int N_POINTS = 30;
		
		
		Histogram[] functions = new Histogram[N_FUNC];
		for (int i = 0; i < N_FUNC; ++i) {
			Random random = new Random();
			
			double x[] = new double[N_POINTS];
			double y[] = new double[N_POINTS];
			for (int j=0; j < N_POINTS; ++j ) {
				x[j] = random.nextDouble();
				y[j] = random.nextDouble();
			}
			Arrays.sort(x);
			Arrays.sort(y);
			functions[i] = PiecewiseLinearFunctions.createFromPoints(x, y);
		}
				
		double sum = 0.;
		for (int i = 0; i < N_ITERATIONS; ++i) {
		
			PiecewiseLinearDistributionIterator iterator = getIterator(functions);
			
			do {
				sum += iterator.getDensity();
			} while(iterator.advance());
		}
		
		
		assertTrue(sum >=0.);
	}*/
	

}
