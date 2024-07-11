package edu.dnu.fpm.regression;

import edu.dnu.fpm.model.RegressionParameter;
import edu.dnu.fpm.service.DataService;
import javafx.util.Pair;

import java.util.List;

public class ExponentialRegression implements IRegression{
    private final int N;
    private final RegressionParameter a;
    private final RegressionParameter b;
    private final LinearRegression linearRegression;
    private final double fTest;
    private final double residualVariance;
    private final double determinationCoefficient;

    public ExponentialRegression(List<Double> xSample, List<Double> ySample){
        N = xSample.size();
        linearRegression = new LinearRegression(
                xSample,
                ySample.stream()
                        .map(Math::log)
                        .toList());
        RegressionParameter intercept = linearRegression.getModel().intercept();
        double aValue = Math.exp(intercept.getValue());
        double studentsQuantile = intercept.getQuantile();
        double as = (aValue - Math.exp(intercept.getIntervalLeftBorder())) / studentsQuantile;
        a = new RegressionParameter(aValue, as, studentsQuantile);
        b = linearRegression.getModel().slope();
        fTest = DataService.getFTest(xSample,ySample, this::getFunctionValue);
        double ys = DataService.getStandardDeviation(ySample, true);
        residualVariance = DataService.getResidualVariance(xSample, ySample, this::getFunctionValue);
        determinationCoefficient = (1 - (N - 2) * residualVariance / (N * ys * ys));
    }

    public RegressionParameter getA() {
        return a;
    }

    public RegressionParameter getB() {
        return b;
    }

    @Override
    public double getFunctionValue(double x) {
        return a.getValue() * Math.exp(b.getValue() * x);
    }

    @Override
    public double getResidualVariance() {
        return residualVariance;
    }

    @Override
    public double getFTest() {
        return fTest;
    }

    @Override
    public double getDeterminationCoefficient() {
        return determinationCoefficient;
    }

    @Override
    public Pair<Double, Double> getFunctionBorders(double x) {
        Pair<Double, Double> borders = linearRegression.getFunctionBorders(x);
        return new Pair<>(Math.exp(borders.getKey()), Math.exp(borders.getValue()));
    }

    @Override
    public Pair<Double, Double> getFunctionValueBorders(double x) {
        Pair<Double, Double> borders = linearRegression.getFunctionValueBorders(x);
        return new Pair<>(Math.exp(borders.getKey()), Math.exp(borders.getValue()));
    }
}
