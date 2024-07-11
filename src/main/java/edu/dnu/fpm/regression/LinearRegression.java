package edu.dnu.fpm.regression;

import edu.dnu.fpm.func.AnalysisFunctions;
import edu.dnu.fpm.model.LinearRegressionModel;
import edu.dnu.fpm.service.DataService;
import javafx.util.Pair;

import java.util.List;

public class LinearRegression implements IRegression{
    private final int N;
    private final LinearRegressionModel model;
    private final double xMean;
    private final double fTest;
    private final double studentsQuantile;

    public LinearRegression(List<Double> xSample, List<Double> ySample) {
        N = xSample.size();
        xMean = DataService.getMean(xSample);
        studentsQuantile = AnalysisFunctions.quantileStudentsDistribution(DataService.P, N - 2);
        model = DataService.getLinearRegression(xSample, ySample);
        fTest = DataService.getFTest(xSample, ySample, this::getFunctionValue);

    }

    public LinearRegressionModel getModel() {
        return model;
    }

    public double getLineRegressionDeviation(double x){
        return Math.sqrt(model.residualVariance() / N +
                Math.pow(model.slope().getStandardDeviation() * (x - xMean), 2));
    }

    @Override
    public double getFunctionValue(double x) {
        return model.intercept().getValue() + model.slope().getValue() * x;
    }

    @Override
    public double getResidualVariance() {
        return model.residualVariance();
    }

    @Override
    public double getFTest() {
        return fTest;
    }

    @Override
    public double getDeterminationCoefficient() {
        return Math.pow(model.pearsonCoefficient(), 2);
    }

    @Override
    public Pair<Double, Double> getFunctionBorders(double x) {
        double f = getFunctionValue(x);
        double ys = getLineRegressionDeviation(x);
        double v = studentsQuantile * ys;
        return new Pair<>(f - v, f + v);
    }

    @Override
    public Pair<Double, Double> getFunctionValueBorders(double x) {
        double f = getFunctionValue(x);
        double ys = getLineRegressionDeviation(x);
        double vs = Math.sqrt(ys * ys + model.residualVariance());
        double v = studentsQuantile * vs;
        return new Pair<>(f - v, f + v);
    }
}
