package org.apache.commons.math3.stat.descriptive.rank;

final class MinHeapUtils {
	
	private MinHeapUtils() {}
	
	static final MinDoubleIntHeap createMinDoubleIntHeap(double[] initialValues) {
		switch(initialValues.length) {
		case 1:
			return new MinDoubleIntHeapSingleValue(initialValues[0]);
		case 2:
			return new MinDoubleIntHeapTwoValues(initialValues[0], initialValues[1]);
		default:
			return new MinDoubleIntHeapGeneral(initialValues);
		}
	}

	// TODO redefine interface, that first updates are guaranteed to return getMinIndex as sequential order
	interface MinDoubleIntHeap {
		
		int getMinIndex();
		
		double getMinValue();
		
		void update(double newValue);
	}
	
	final static class MinDoubleIntHeapGeneral implements MinDoubleIntHeap {
		
		private final double[] values;
		private final int[] indices;
		
		public MinDoubleIntHeapGeneral(double[] values) {
			int n = values.length;
			
			this.values = new double[n];
			this.indices = new int[n];
			for (int i = 0; i < n; ++i) {
				this.values[i] = Double.NEGATIVE_INFINITY;
			}
			for (int i = 0; i < n; ++i) {
				this.indices[0] = i;
				update(values[i]);
			}
		}
		
		public int getMinIndex() {
			return indices[0];
		}
		
		public double getMinValue() {
			return values[0];
		}
		
		public void update(double newValue) {
			
			final int n = values.length;
			
			final int updatedIndex = indices[0];
			
			int parentIdx = 0;
			
			while(true) {
				
				int leftChildIdx = (parentIdx << 1) | 1;			
				
				if (leftChildIdx < n-1) {
					
					final int rightChildIdx = leftChildIdx+1;
					
					final double leftChildValue = values[leftChildIdx];
					final double rightChildValue = values[rightChildIdx];
					
					final int minChildIdx;
					final double minChildValue;
					if (rightChildValue < leftChildValue) {
						minChildIdx = rightChildIdx;
						minChildValue = rightChildValue;
					}
					else {
						minChildIdx = leftChildIdx;
						minChildValue = leftChildValue;
					}
					if (minChildValue < newValue) {
						values[parentIdx] = minChildValue;
						indices[parentIdx] = indices[minChildIdx];
						parentIdx = minChildIdx;
						continue;
					}
				} else if (leftChildIdx == n-1) {
					final double leftChildValue = values[leftChildIdx];
					if (leftChildValue < newValue) {
						values[parentIdx] = leftChildValue;
						indices[parentIdx] = indices[leftChildIdx];
						parentIdx = leftChildIdx;
					}
				}
				break;
			}
			values[parentIdx] = newValue;
			indices[parentIdx] = updatedIndex;		
		}
	}
	
	final static class MinDoubleIntHeap2 implements MinDoubleIntHeap {
		
		private final double[] values;
		private final int[] indices;
		
		public MinDoubleIntHeap2(double[] values) {
			int nt = values.length;
			
			int n = (Integer.highestOneBit(nt)<<1)-1;
			
			
			
			this.values = new double[n];
			this.indices = new int[n];
			for (int i = 0; i < nt; ++i) {
				this.values[i] = Double.NEGATIVE_INFINITY;
			}
			for (int i = nt; i < n; ++i) {
				this.values[i] = Double.NaN;
			}
			for (int i = 0; i < nt; ++i) {
				this.indices[0] = i;
				update(values[i]);
			}
		}
		
		public int getMinIndex() {
			return indices[0];
		}
		
		public double getMinValue() {
			return values[0];
		}
		
		public void update(double newValue) {
			
			final int n = values.length;
			
			final int updatedIndex = indices[0];
			
			int parentIdx = 0;
			
			for(int i = n; i!=1; i>>=1 ) {
				
				final int leftChildIdx = (parentIdx << 1) | 1;			
				final int rightChildIdx = leftChildIdx+1;
				
				final double leftChildValue = values[leftChildIdx];
				final double rightChildValue = values[rightChildIdx];
				
				final int minChild;
				final double minChildValue;
				if (rightChildValue < leftChildValue) {
					minChild = rightChildIdx;
					minChildValue = rightChildValue;
				}
				else {
					minChild = leftChildIdx;
					minChildValue = leftChildValue;
				}
				if (minChildValue < newValue) {
					values[parentIdx] = minChildValue;
					indices[parentIdx] = indices[minChild];
					parentIdx = minChild;
				}
				else {
					break;
				}
			}
			values[parentIdx] = newValue;
			indices[parentIdx] = updatedIndex;		
		}
	}
	
	final static class MinDoubleIntHeapSingleValue implements MinDoubleIntHeap {
		
		private double value;
		
		public MinDoubleIntHeapSingleValue(double value) {
			this.value = value;
		}
		
		public int getMinIndex() {
			return 0;
		}
		
		public void update(double newValue) {
			value = newValue;
		}

		public double getMinValue() {
			return value;
		}
	}
	
	final static class MinDoubleIntHeapTwoValues implements MinDoubleIntHeap {
		
		private double value0;
		private double value1;
		private int minIndex;
		
		public MinDoubleIntHeapTwoValues(double value0, double value1) {
			this.value0 = value0;
			this.value1 = value1;
			minIndex = (value0 <= value1)?0:1;
		}
		
		public int getMinIndex() {
			return minIndex;
		}
		
		public double getMinValue() {
			return (minIndex==0)?value0:value1;
		}
		
		public void update(double newValue) {
			if (minIndex==0) {
				value0 = newValue;
			}
			else {
				value1 = newValue;
			}
			minIndex = (value0 <= value1)?0:1;
		}
	}

}
