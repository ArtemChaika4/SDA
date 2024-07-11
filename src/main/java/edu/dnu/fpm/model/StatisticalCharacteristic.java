package edu.dnu.fpm.model;

public class StatisticalCharacteristic {
    private static final String BLANK = "—";
    private final String name;
    private final String value;
    private final String standardDeviation;
    private final String interval;

    public StatisticalCharacteristic(String name, Double value, Double standardDeviation, Double leftBorder, Double rightBorder) {
        this.name = name;
        this.value = String.format("%.3f", value);
        this.standardDeviation = (standardDeviation == null) ? BLANK : String.format("%.3f", standardDeviation);
        this.interval = (leftBorder == null || rightBorder == null) ? BLANK :
                String.format("[%.3f;%.3f]", leftBorder, rightBorder);
    }

    public String getName() {
        return name;
    }

    public String getValue() {
        return value;
    }

    public String getStandardDeviation() {
        return standardDeviation;
    }

    public String getInterval() {
        return interval;
    }
}
