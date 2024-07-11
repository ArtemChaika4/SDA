package edu.dnu.fpm.model;


public class CorrelationCoefficient {
    private static final String BLANK = "—";
    private String name;
    private String value;
    private String interval;
    private String test;
    private String quantile;
    private String significance;
    private String result;
    private String standardDeviation;

    public CorrelationCoefficient(String name, RegressionParameter parameter){
        this(name, parameter.getValue(), parameter.getIntervalLeftBorder(), parameter.getIntervalRightBorder(), parameter.getTest(), parameter.getQuantile());
        this.standardDeviation = String.format("%.3f", parameter.getStandardDeviation());
    }

    public CorrelationCoefficient(String name, double value, double leftBorder, double rightBorder,
                                  double test, double quantile) {
        this(name, value, test, quantile);
        this.interval = String.format("[%.3f;%.3f]", leftBorder, rightBorder);
    }

    public String getStandardDeviation() {
        return standardDeviation;
    }

    public CorrelationCoefficient(String name, double value, double test, double quantile) {
        this.name = name;
        this.value = String.format("%.3f", value);
        this.interval = BLANK;
        this.test = String.format("%.3f", test);
        this.quantile = String.format("%.3f", quantile);
        boolean significance = Math.abs(test) > quantile;
        this.significance = significance ? "значущий" : "незначущий";
        this.result = significance ? "є" : "немає";
    }

    public String getName() {
        return name;
    }

    public String getValue() {
        return value;
    }

    public String getInterval() {
        return interval;
    }

    public String getTest() {
        return test;
    }

    public String getQuantile() {
        return quantile;
    }

    public String getSignificance() {
        return significance;
    }

    public String getResult() {
        return result;
    }
}
