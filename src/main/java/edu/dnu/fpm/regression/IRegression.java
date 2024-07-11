package edu.dnu.fpm.regression;

import javafx.util.Pair;

public interface IRegression {
    double getFunctionValue(double x);
    double getResidualVariance();
    double getFTest();
    double getDeterminationCoefficient();
    Pair<Double, Double> getFunctionBorders(double x);
    Pair<Double, Double> getFunctionValueBorders(double x);

}
