package edu.dnu.fpm.model;

public class RegressionParameter {
    private double value;
    private double standardDeviation;
    private double quantile;
    private String name;

    public RegressionParameter(){}
    public RegressionParameter(double value, double standardDeviation, double quantile) {
        this.value = value;
        this.standardDeviation = standardDeviation;
        this.quantile = quantile;
    }

    public double getIntervalLeftBorder(){
        return value - quantile * standardDeviation;
    }

    public double getIntervalRightBorder(){
        return value + quantile * standardDeviation;
    }

    public double getQuantile() {
        return quantile;
    }

    public double getValue() {
        return value;
    }

    public double getStandardDeviation() {
        return standardDeviation;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public double getTest() {
        return value / standardDeviation;
    }

}
