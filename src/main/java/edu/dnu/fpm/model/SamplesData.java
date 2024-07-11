package edu.dnu.fpm.model;

import java.util.ArrayList;
import java.util.List;

public class SamplesData {
    private String firstColumnName;
    private String secondColumnName;
    private List<Double> firstSample;
    private List<Double> secondSample;

    public SamplesData(){
        firstSample = new ArrayList<>();
        secondSample = new ArrayList<>();
    }

    public SamplesData(List<Double> firstSample, List<Double> secondSample) {
        this.firstSample = firstSample;
        this.secondSample = secondSample;
    }

    public SamplesData(String firstColumnName, String secondColumnName, List<Double> firstSample, List<Double> secondSample) {
        this.firstColumnName = firstColumnName;
        this.secondColumnName = secondColumnName;
        this.firstSample = firstSample;
        this.secondSample = secondSample;
    }

    public List<Double> getFirstSample() {
        return firstSample;
    }

    public List<Double> getSecondSample() {
        return secondSample;
    }

    public void setFirstSample(List<Double> firstSample) {
        this.firstSample = firstSample;
    }

    public String getFirstColumnName() {
        return firstColumnName;
    }

    public void setFirstColumnName(String firstColumnName) {
        this.firstColumnName = firstColumnName;
    }

    public String getSecondColumnName() {
        return secondColumnName;
    }

    public void setSecondColumnName(String secondColumnName) {
        this.secondColumnName = secondColumnName;
    }

    public void setSecondSample(List<Double> secondSample) {
        this.secondSample = secondSample;
    }
}
