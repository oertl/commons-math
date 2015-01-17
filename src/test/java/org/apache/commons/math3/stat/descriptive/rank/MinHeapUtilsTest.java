package org.apache.commons.math3.stat.descriptive.rank;

import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.stat.descriptive.rank.MinHeapUtils.MinDoubleIntHeap;
import org.junit.Test;

public class MinHeapUtilsTest {
	
	@Test
	public void testMinDoubleIntHeap1() {
		
		MinDoubleIntHeap heap = MinHeapUtils.createMinDoubleIntHeap(new double[]{2., 1., 3., 5., 4.});		
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
		
		MinDoubleIntHeap heap = MinHeapUtils.createMinDoubleIntHeap(new double[]{2., 1.});		
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
		
		MinDoubleIntHeap heap = MinHeapUtils.createMinDoubleIntHeap(new double[]{3., 1.});		
		assertEquals(1, heap.getMinIndex());
		assertEquals(1., heap.getMinValue(), 0.);
		
		heap.update(2.);
		assertEquals(1, heap.getMinIndex());
		assertEquals(2., heap.getMinValue(), 0.);
	}
	
	@Test
	public void testMinDoubleIntHeap4() {
		
		MinDoubleIntHeap heap = MinHeapUtils.createMinDoubleIntHeap(new double[]{2., 1., 3., 5., 4.});		
		assertEquals(1, heap.getMinIndex());
		assertEquals(1., heap.getMinValue(), 0.);
		
		heap.update(1.5);
		assertEquals(1, heap.getMinIndex());
		assertEquals(1.5, heap.getMinValue(), 0.);
	}

}
