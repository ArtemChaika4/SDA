package edu.dnu.fpm.model;

public record LinearRegressionModel(RegressionParameter intercept, RegressionParameter slope,
                                    double pearsonCoefficient, double residualVariance) {

}
