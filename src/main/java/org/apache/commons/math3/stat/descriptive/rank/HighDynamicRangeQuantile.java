package org.apache.commons.math3.stat.descriptive.rank;


public class HighDynamicRangeQuantile {
	private final long[] counts;
	private double minimum = Double.POSITIVE_INFINITY;
	private double maximum = Double.NEGATIVE_INFINITY;
	private long underFlowCount = 0;
	private long overFlowCount = 0;	
	private final double factor;
	private final double offset;
	private final double minExpectedQuantileValue;
	
	public HighDynamicRangeQuantile(final double minExpectedQuantileValue, final double maxExpectedQuantileValue, final double maxRelativeError) {
		
		// TOOD range checks
		
		final int arraySize = getIndex(maxExpectedQuantileValue) + 1;
		counts = new long[arraySize];
		factor = 0.25/Math.log(maxRelativeError);
		offset = getIndexHelper(minExpectedQuantileValue);
		this.minExpectedQuantileValue = minExpectedQuantileValue;
	}
		
	int getIndex(final double value) {
		return (int) (getIndexHelper(value) - offset);
	}
	
	double getIndexHelper(final double value) {
		final long valueBits = Double.doubleToRawLongBits(value);
		final long exponent = (valueBits & 0x7ff0000000000000L) >> 51;
		final double exponent_mul_3 = exponent + (exponent >> 1);
		final double mantissa = Double.longBitsToDouble((valueBits & 0x800fffffffffffffL) | 0x3ff0000000000000L);
		return factor*((mantissa-1.)*(5.-mantissa)+exponent_mul_3);
	}

	public void add(final double value, final long count) {
		
		
		if (value >= minExpectedQuantileValue) {
			final int idx = getIndex(value);
			if (idx < counts.length) {
				counts[idx] += count;
			}
			else {
				overFlowCount += count;
			}
		}
		else {
			underFlowCount += count;
		}
		minimum = Math.min(minimum, value); 
		maximum = Math.max(maximum, value);
	}
	
	public long getCount() {
		long sum = underFlowCount + overFlowCount;
		for (final long c : counts) {
			sum += c;
		}
		return sum;
	}
	

}
