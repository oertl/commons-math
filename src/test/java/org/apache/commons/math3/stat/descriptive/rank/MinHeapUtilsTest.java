package org.apache.commons.math3.stat.descriptive.rank;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.apache.commons.math3.stat.descriptive.rank.MinHeapUtils.MinDoubleIntHeap;
import org.apache.commons.math3.stat.descriptive.rank.MinHeapUtils.MinDoubleIntHeapGeneral;
import org.apache.commons.math3.stat.descriptive.rank.MinHeapUtils.MinDoubleIntHeapSingleValue;
import org.apache.commons.math3.stat.descriptive.rank.MinHeapUtils.MinDoubleIntHeapTwoValues;
import org.junit.Test;

public class MinHeapUtilsTest {
	
	@Test
	public void testMinDoubleIntHeapGeneral1() {
		
		MinDoubleIntHeap heap = new MinDoubleIntHeapGeneral(new double[]{2., 1., 3., 5., 4.});		
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
	public void testMinDoubleIntHeapGeneral2() {
		
		MinDoubleIntHeap heap = new MinDoubleIntHeapGeneral(new double[]{2., 1.});		
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
	public void testMinDoubleIntHeapGeneral3() {
		
		MinDoubleIntHeap heap = new MinDoubleIntHeapGeneral(new double[]{3., 1.});		
		assertEquals(1, heap.getMinIndex());
		assertEquals(1., heap.getMinValue(), 0.);
		
		heap.update(2.);
		assertEquals(1, heap.getMinIndex());
		assertEquals(2., heap.getMinValue(), 0.);
	}
	
	@Test
	public void testMinDoubleIntHeapGeneral4() {
		
		MinDoubleIntHeap heap = new MinDoubleIntHeapGeneral(new double[]{2., 1., 3., 5., 4.});		
		assertEquals(1, heap.getMinIndex());
		assertEquals(1., heap.getMinValue(), 0.);
		
		heap.update(1.5);
		assertEquals(1, heap.getMinIndex());
		assertEquals(1.5, heap.getMinValue(), 0.);
	}
	
	@Test
	public void testMinDoubleIntHeapSingleValue1() {
		MinDoubleIntHeap heap = new MinDoubleIntHeapSingleValue(2.);
		assertEquals(0, heap.getMinIndex());
		assertEquals(2., heap.getMinValue(), 0.);
		
		heap.update(1.5);
		assertEquals(0, heap.getMinIndex());
		assertEquals(1.5, heap.getMinValue(), 0.);		
	}
	
	@Test
	public void testMinDoubleIntHeapTwoValues1() {
		MinDoubleIntHeap heap = new MinDoubleIntHeapTwoValues(2., 4.);
		assertEquals(0, heap.getMinIndex());
		assertEquals(2., heap.getMinValue(), 0.);
		
		heap.update(5.);
		assertEquals(1, heap.getMinIndex());
		assertEquals(4., heap.getMinValue(), 0.);		

		heap.update(4.5);
		assertEquals(1, heap.getMinIndex());
		assertEquals(4.5, heap.getMinValue(), 0.);		
	}

	@Test
	public void testMinDoubleIntHeapTwoValues2() {
		MinDoubleIntHeap heap = new MinDoubleIntHeapTwoValues(7., 4.);
		assertEquals(1, heap.getMinIndex());
		assertEquals(4., heap.getMinValue(), 0.);
		
		heap.update(8.);
		assertEquals(0, heap.getMinIndex());
		assertEquals(7., heap.getMinValue(), 0.);		

		heap.update(4.5);
		assertEquals(0, heap.getMinIndex());
		assertEquals(4.5, heap.getMinValue(), 0.);		
	}
	
	@Test
	public void testCreateMinDoubleIntHeap() {
		assertTrue(MinHeapUtils.createMinDoubleIntHeap(new double[0]) instanceof MinDoubleIntHeapGeneral);
		assertTrue(MinHeapUtils.createMinDoubleIntHeap(new double[1]) instanceof MinDoubleIntHeapSingleValue);
		assertTrue(MinHeapUtils.createMinDoubleIntHeap(new double[2]) instanceof MinDoubleIntHeapTwoValues);
		assertTrue(MinHeapUtils.createMinDoubleIntHeap(new double[3]) instanceof MinDoubleIntHeapGeneral);
		assertTrue(MinHeapUtils.createMinDoubleIntHeap(new double[4]) instanceof MinDoubleIntHeapGeneral);
		assertTrue(MinHeapUtils.createMinDoubleIntHeap(new double[5]) instanceof MinDoubleIntHeapGeneral);
	}

}
