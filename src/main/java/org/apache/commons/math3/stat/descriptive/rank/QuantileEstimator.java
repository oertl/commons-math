package org.apache.commons.math3.stat.descriptive.rank;

public interface QuantileEstimator {
	
	void add(double value);
	
	void add(double value, long count);
	
	double quantile(double pValue);
	
	double probability(double from, double to);
	
	long frequency(double from, double to);
	
	double minimum();
	
	double maximum();
	
	long count();
	
}
