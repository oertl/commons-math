package org.apache.commons.math3.stat.descriptive.rank;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertArrayEquals;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.rank.IQAgentQuantile.DynamicSum;
import org.apache.commons.math3.stat.descriptive.rank.IQAgentQuantile.Histogram;
import org.apache.commons.math3.stat.descriptive.rank.IQAgentQuantile.HistogramIterator;
import org.apache.commons.math3.stat.descriptive.rank.IQAgentQuantile.HistogramIterator1;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.MathArrays.OrderDirection;
import org.junit.Ignore;
import org.junit.Test;

public class IQAgentQuantileTest {
	
	@Ignore
	@Test
	public void testPerformance() {
		int numPValues = 30;
		double[] pValues = new double[numPValues];
		for (int i = 0; i < numPValues; ++i) {
			pValues[i] = i/(numPValues-1.);
		}
			
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
		
		histograms.add(IQAgentQuantile.sortedValuesAsHistogram(new double[]{2.7}, 1, 3.));
		
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
	public void testHistogramIterator3() {

		final double[] pValues = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
		
		final double[] values1 = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
		Arrays.sort(values1);
		final Histogram histogram1 = IQAgentQuantile.sortedValuesAsHistogram(values1, values1.length, 10./170.);
		
		double[] quantiles = {0.0, 0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.0};
		final Histogram histogram2 = IQAgentQuantile.asHistogram(pValues, quantiles, 160, 160./170.);
		HistogramIterator iterator = new IQAgentQuantile.HistogramIterator1(Arrays.asList(histogram1, histogram2));
		
		double sum1 = 0.;
		double sum2 = 0.;
		
		while(iterator.advance()) {
			sum1 += iterator.getMaximumCumulativeFrequency() - iterator.getMinimumCumulativeFrequency();
			if (iterator.getDensity()!=Double.POSITIVE_INFINITY) {				
				sum2 += (iterator.getMaximumValue() - iterator.getMinimumValue())*iterator.getDensity();
			}
			else {
				sum2 += iterator.getMaximumCumulativeFrequency() - iterator.getMinimumCumulativeFrequency();
			}
		}
		
		assertEquals(1., sum1, 1e-8);
		assertEquals(1., sum2, 1e-8);
		
	}
		
	@Test
	public void testSortedValuesAsHistogram() {
		
		final double[] values = new double[]{1., 4., 5.};
		
		final double scale = 7.;
		
		final Histogram histogram = IQAgentQuantile.sortedValuesAsHistogram(values, values.length, scale);
		
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
	public void testEvaluateSumOfHistograms3() {
		
		final double[] values = {0.24053641567148587, 0.3332183994766498, 0.3851891847407185, 0.5504370051176339, 0.5975452777972018, 0.6374174253501083, 0.730967787376657, 0.8791825178724801, 0.9412491794821144, 0.984841540199809};
		
		final Histogram histogram = IQAgentQuantile.sortedValuesAsHistogram(values, values.length, 1.);
		
		final double[] cumulativeFrequencies = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

		double[] result = IQAgentQuantile.evaluateSumOfHistograms(Collections.singleton(histogram), cumulativeFrequencies);

		MathArrays.checkOrder(result, OrderDirection.INCREASING, false);	
	}
	
	/*@Test
	public void testEvaluateSumOfHistograms4() {
		

		final double[] pValues = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
		
		
		//final double[] values1 = {0.19990322644759395, 0.23665037197274819, 0.27401999592272896, 0.15905061152187727, 0.3724993797577332, 0.2157097180385359, 0.06491450914260377, 0.609041308147579, 0.43370181810972497, 0.18927685453270326};
		final double[] values1 = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
		Arrays.sort(values1);
		final Histogram histogram1 = IQAgentQuantile.sortedValuesAsHistogram(values1, 10./170.);
		
		//final double[] quantiles2 = {8.749017200373466E-4, 0.08386514655696836, 0.18491298087597574, 0.286650928807526, 0.40768478187592994, 0.5312325427273474, 0.6287885299580311, 0.7180895578947304, 0.792670973240814, 0.9460722563317701, 0.9949331294605093};
		
		double[] quantiles = {0.0, 0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.0};
		for (int i = 0; i < 100; ++i) {
			final Histogram histogram2 = IQAgentQuantile.asHistogram(pValues, quantiles, 160, 160./170.);			
			quantiles = IQAgentQuantile.evaluateSumOfHistograms(Arrays.asList(histogram1, histogram2), pValues);
		}
		
		
		System.out.println(Arrays.toString(quantiles));
		
	}*/
	
	
	@Test
	public void testEvaluateSumOfHistrogramsStability() {
		

		final double[] pValues = {0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.51, 0.6, 0.7, 0.8, 0.9, 1.0};
		
		double[] quantiles = {0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.51, 0.6, 0.7, 0.8, 0.9, 1.0};
		
		for (int i = 0; i< 1000;++i) {
			Histogram histogram = IQAgentQuantile.asHistogram(pValues, quantiles, 10, 1.);
			quantiles = IQAgentQuantile.evaluateSumOfHistograms(Collections.singleton(histogram), pValues);
		}
		
		assertArrayEquals(pValues, quantiles, 1e-30);
	}
	/*
	IQAgentQuantile [buffer=[0.19990322644759395, 0.23665037197274819, 0.27401999592272896, 0.15905061152187727, 0.3724993797577332, 0.2157097180385359, 0.06491450914260377, 0.609041308147579, 0.43370181810972497, 0.7180895578947304], bufferCounter=9, pValues=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], quantiles=[8.749017200373466E-4, 0.08386514655696836, 0.18491298087597574, 0.286650928807526, 0.40768478187592994, 0.5312325427273474, 0.6287885299580311, 0.7180895578947304, 0.792670973240814, 0.9460722563317701, 0.9949331294605093], totalCount=169]*/
	/*IQAgentQuantile [buffer=[0.19990322644759395, 0.23665037197274819, 0.27401999592272896, 0.15905061152187727, 0.3724993797577332, 0.2157097180385359, 0.06491450914260377, 0.609041308147579, 0.43370181810972497, 0.18927685453270326], bufferCounter=0, pValues=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], quantiles=[8.749017200373466E-4, 0.08386514655696836, 0.18491298087597574, 0.26757506357036037, 0.3774263186088291, 0.5003456025144932, 0.609041308147579, 0.7192022658818833, 0.8045472044624453, 0.9531293740131664, 0.9949331294605093], totalCount=170]*/
	

	@Test
	public void testEvaluateSumOfHistogramsForTwoValues() {
		
		final double min = -4.7;
		final double max = 5.7;
			
		final Histogram histogram = IQAgentQuantile.sortedValuesAsHistogram(new double[]{min,max}, 2, 2.);
		
		final double[] cumulativeFrequencies = {-1., 0., 1., 2., 3.};
		final double[] expectedValues = {min, min, (min+max)*0.5, max, max};
		
		double[] values = IQAgentQuantile.evaluateSumOfHistograms(Collections.singleton(histogram), cumulativeFrequencies);
		
		assertArrayEquals(expectedValues, values, 1e-30);
	}
	
	@Test
	public void testEvaluateSumOfHistogramsForTwoValues2() {
		
		final double min = -4.7;
		final double max = 5.7;
			
		final Histogram histogram = IQAgentQuantile.sortedValuesAsHistogram(new double[]{min,max}, 2, 2.);
		
		final double[] cumulativeFrequencies = {-1., 0., 1.};
		final double[] expectedValues = {min, min, (min+max)*0.5};
		
		double[] values = IQAgentQuantile.evaluateSumOfHistograms(Collections.singleton(histogram), cumulativeFrequencies);
		
		assertArrayEquals(expectedValues, values, 1e-30);
	}
	
	@Test
	public void testAddSingleValue() {
		
		IQAgentQuantile quantileEstimator = new IQAgentQuantile(new double[]{0., 1.}, 1);
		
		final double value = 0.5;
		
		quantileEstimator.add(value);
		
		assertEquals(quantileEstimator.getQuantile(0.), value, 0.0);
		assertEquals(quantileEstimator.getQuantile(0.1), value, 0.0);
		assertEquals(quantileEstimator.getQuantile(0.2), value, 0.0);
		assertEquals(quantileEstimator.getQuantile(0.3), value, 0.0);
		assertEquals(quantileEstimator.getQuantile(0.4), value, 0.0);
		assertEquals(quantileEstimator.getQuantile(0.5), value, 0.0);
		assertEquals(quantileEstimator.getQuantile(0.6), value, 0.0);
		assertEquals(quantileEstimator.getQuantile(0.7), value, 0.0);
		assertEquals(quantileEstimator.getQuantile(0.8), value, 0.0);
		assertEquals(quantileEstimator.getQuantile(0.9), value, 0.0);
		assertEquals(quantileEstimator.getQuantile(1.0), value, 0.0);
	}
	
	@Test
	public void testAddZeroValues() {
		
		final IQAgentQuantile quantileEstimator = new IQAgentQuantile(new double[]{0., 1.}, 1);
		
		assertEquals(quantileEstimator.getQuantile(0.0), Double.NaN, 0.0);
		assertEquals(quantileEstimator.getQuantile(0.5), Double.NaN, 0.0);
		assertEquals(quantileEstimator.getQuantile(1.0), Double.NaN, 0.0);
	}
	
	@Test
	public void testAddTwoValues() {
		
		final IQAgentQuantile quantileEstimator = new IQAgentQuantile(new double[]{0., 0.5, 1.}, 2);
		
		final double value1 = 545.4545;
		final double value2 = 765.1234;
		
		final double eps = 1e-8;
		
		quantileEstimator.add(value1);
		quantileEstimator.add(value2);
		
		assertEquals(value1, quantileEstimator.getQuantile(0.0), eps);
		assertEquals(value1, quantileEstimator.getQuantile(0.1), eps);
		assertEquals(value1, quantileEstimator.getQuantile(0.2), eps);
		assertEquals(value1 + 0.1*(value2 - value1), quantileEstimator.getQuantile(0.3), eps);
		assertEquals(value1 + 0.3*(value2 - value1), quantileEstimator.getQuantile(0.4), eps);
		assertEquals(value1 + 0.5*(value2 - value1), quantileEstimator.getQuantile(0.5), 1e-30);
		assertEquals(value1 + 0.7*(value2 - value1), quantileEstimator.getQuantile(0.6), eps);
		assertEquals(value1 + 0.9*(value2 - value1), quantileEstimator.getQuantile(0.7), eps);
		assertEquals(value2, quantileEstimator.getQuantile(0.8), eps);
		assertEquals(value2, quantileEstimator.getQuantile(0.9), eps);
		assertEquals(value2, quantileEstimator.getQuantile(1.0), eps);
	}
	
	// TODO test histogram iterator for histogram made of sorted values with identical values
	
	@Test
	public void testAddManyValues() {
		
		Random random = new Random(0);
		int numValues = 100000;
		
		double[] values = new double[numValues];
		for (int i = 0; i < numValues; ++i) {
			values[i] = random.nextDouble();
		}
		
		final double[] pValues = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
		
		IQAgentQuantile quantileEstimator = new IQAgentQuantile(pValues, 10);
		
		for (double value : values) {
			quantileEstimator.add(value);
		}
		
		for (double pValue : pValues) {
			assertEquals(pValue, quantileEstimator.getQuantile(pValue), 0.01);
		}
	}
	
	@Test
	public void testAsHistogram() {
		
		double[] pValues = {0., 0.3, 0.7, 0.8, 1.0};
		double[] quantiles = {-9., -3., 1., 7., 8.};
		
		long count = 4;
		double scale = 2.;
		
		double eps = 1e-8;
		
		final Histogram histogram = IQAgentQuantile.asHistogram(pValues, quantiles, count, scale);
		
		assertEquals(histogram.getNumberOfBins(), 6);
		
		assertEquals(-9., histogram.getBoundary(0), eps);
		assertEquals(-9., histogram.getBoundary(1), eps);
		assertEquals(-3., histogram.getBoundary(2), eps);
		assertEquals(1., histogram.getBoundary(3), eps);
		assertEquals(7., histogram.getBoundary(4), eps);
		assertEquals(8., histogram.getBoundary(5), eps);
		assertEquals(8., histogram.getBoundary(6), eps);
		
		assertEquals(0.25, histogram.getFrequency(0), eps);
		assertEquals(0.35, histogram.getFrequency(1), eps);
		assertEquals(0.8, histogram.getFrequency(2), eps);
		assertEquals(0.2, histogram.getFrequency(3), eps);
		assertEquals(0.15, histogram.getFrequency(4), eps);
		assertEquals(0.25, histogram.getFrequency(5), eps);	
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
