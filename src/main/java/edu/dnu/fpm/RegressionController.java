package edu.dnu.fpm;

import edu.dnu.fpm.func.AnalysisFunctions;
import edu.dnu.fpm.model.*;
import edu.dnu.fpm.regression.ExponentialRegression;
import edu.dnu.fpm.regression.IRegression;
import edu.dnu.fpm.regression.LinearRegression;
import edu.dnu.fpm.service.DataService;
import javafx.fxml.Initializable;
import javafx.scene.control.*;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.util.Pair;

import java.net.URL;
import java.util.List;
import java.util.ResourceBundle;

public class RegressionController implements Initializable {
    public Label funcLabel;
    public TableView<CorrelationCoefficient> table;
    public TableColumn<CorrelationCoefficient, String> name;
    public TableColumn<CorrelationCoefficient, String> value;
    public TableColumn<CorrelationCoefficient, String> standardDeviation;
    public TableColumn<CorrelationCoefficient, String> interval;
    public TableColumn<CorrelationCoefficient, String> quantile;
    public TableColumn<CorrelationCoefficient, String> significance;
    public TableColumn<CorrelationCoefficient, String> test;
    public TextField residualVarianceField;
    public TextField determinationCoefficientField;
    public Label fTestValue;
    public Label fTestQuantile;
    public Label fTestResult;
    public TextField pointField;
    public Label pointValue;
    public Label pointInterval;
    public Button expButton;
    private List<Double> firstSample;
    private List<Double> secondSample;
    private IRegression regression;
    private LinearRegression linearRegression;
    private ExponentialRegression exponentialRegression;

    public void setOnLinear() {
        if(linearRegression == null){
            linearRegression = new LinearRegression(firstSample, secondSample);
        }
        regression = linearRegression;
        table.getItems().clear();
        table.getItems().add(new CorrelationCoefficient("a", linearRegression.getModel().intercept()));
        table.getItems().add(new CorrelationCoefficient("b", linearRegression.getModel().slope()));
        setRegressionCoefficients(linearRegression);
        funcLabel.setText("y = a + b * x");
    }

    public void setOnExponential() {
        if(exponentialRegression == null){
            exponentialRegression = new ExponentialRegression(firstSample, secondSample);
        }
        regression = exponentialRegression;
        table.getItems().clear();
        table.getItems().add(new CorrelationCoefficient("a", exponentialRegression.getA()));
        table.getItems().add(new CorrelationCoefficient("b", exponentialRegression.getB()));
        setRegressionCoefficients(exponentialRegression);
        funcLabel.setText("y = a * exp(b * x)");
    }

    private void setRegressionCoefficients(IRegression regression){
        residualVarianceField.setText(String.format("%.3f", regression.getResidualVariance()));
        determinationCoefficientField.setText(String.format("%.3f", regression.getDeterminationCoefficient() * 100));
        double fTest = regression.getFTest();
        double q = AnalysisFunctions.quantileFishersDistribution(1 - DataService.A, 1, firstSample.size() - 2);
        fTestValue.setText(String.format("%.3f", fTest));
        fTestQuantile.setText(String.format("%.3f", q));
        fTestResult.setText(fTest > q ? "значуща" : "незначуща");
        pointValue.setText("");
        pointInterval.setText("");
    }

    public void calculate() {
        double x;
        try {
            x = Double.parseDouble(pointField.getText().replace(',', '.'));
            pointField.setStyle("-fx-border-color: none;");
        }catch (NumberFormatException e){
            pointField.setStyle("-fx-border-color: red;");
            return;
        }
        double y = regression.getFunctionValue(x);
        Pair<Double, Double> p = regression.getFunctionValueBorders(x);
        pointValue.setText(String.format("%.3f", y));
        pointInterval.setText(String.format("[%.3f;%.3f]", p.getKey(), p.getValue()));
    }

    @Override
    public void initialize(URL url, ResourceBundle resourceBundle) {
        SamplesData samplesData = (SamplesData)resourceBundle.getObject("");
        firstSample = samplesData.getFirstSample();
        secondSample = samplesData.getSecondSample();

        name.setCellValueFactory(new PropertyValueFactory<>("name"));
        value.setCellValueFactory(new PropertyValueFactory<>("value"));
        standardDeviation.setCellValueFactory(new PropertyValueFactory<>("standardDeviation"));
        interval.setCellValueFactory(new PropertyValueFactory<>("interval"));
        test.setCellValueFactory(new PropertyValueFactory<>("test"));
        quantile.setCellValueFactory(new PropertyValueFactory<>("quantile"));
        significance.setCellValueFactory(new PropertyValueFactory<>("significance"));
        double yMin = secondSample.stream().mapToDouble(y -> y).min().orElse(0);
        if(yMin <= 0){
            expButton.setDisable(true);
        }
        setOnLinear();
    }
}
