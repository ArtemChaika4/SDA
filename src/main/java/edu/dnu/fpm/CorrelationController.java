package edu.dnu.fpm;

import com.jfoenix.controls.JFXToggleButton;
import edu.dnu.fpm.func.AnalysisFunctions;
import edu.dnu.fpm.model.CorrelationCoefficient;
import edu.dnu.fpm.model.SamplesData;
import edu.dnu.fpm.model.StatisticalCharacteristic;
import edu.dnu.fpm.regression.ExponentialRegression;
import edu.dnu.fpm.regression.IRegression;
import edu.dnu.fpm.regression.LinearRegression;
import edu.dnu.fpm.service.DataService;
import javafx.beans.property.SimpleStringProperty;
import javafx.collections.FXCollections;
import javafx.fxml.Initializable;
import javafx.geometry.Pos;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.*;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.layout.StackPane;
import javafx.scene.text.Text;
import javafx.scene.text.TextAlignment;
import javafx.util.Pair;

import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.ResourceBundle;

public class CorrelationController implements Initializable {
    public StackPane chartPane;
    public TableView<StatisticalCharacteristic> table;
    public TableColumn<StatisticalCharacteristic, String> name;
    public TableColumn<StatisticalCharacteristic, String> value;
    public TableColumn<StatisticalCharacteristic, String> standardDeviation;
    public TableColumn<StatisticalCharacteristic, String> interval;
    public Label skewnessLabel;
    public Label kurtosisLabel;
    public Label resultLabel;
    public TableView<CorrelationCoefficient> coefficientTable;
    public TableColumn<CorrelationCoefficient, String> coefficientName;
    public TableColumn<CorrelationCoefficient, String> coefficientValue;
    public TableColumn<CorrelationCoefficient, String> coefficientInterval;
    public TableColumn<CorrelationCoefficient, String> coefficientTest;
    public TableColumn<CorrelationCoefficient, String> coefficientQuantile;
    public TableColumn<CorrelationCoefficient, String> coefficientSignificance;
    public TableColumn<CorrelationCoefficient, String> coefficientResult;
    public TableView<?> ratioLineTable;
    public TableColumn<?, String> pearsonCoefficientValue;
    public TableColumn<?, String> correlationRatioValue;
    public TableColumn<?, String> ratioLineTest;
    public TableColumn<?, String> ratioLineQuantile;
    public TableColumn<?, String> ratioLineEquality;
    public TableColumn<?, String> ratioLineResult;
    public JFXToggleButton linToggle;
    public JFXToggleButton expToggle;
    private List<Double> firstSample;
    private List<Double> secondSample;
    private SamplesData samplesData;
    private List<XYChart.Series<Number, Number>> linearRegressionSeries;
    private List<XYChart.Series<Number, Number>> exponentialRegressionSeries;
    private LineChart<Number, Number> lineChart;

    @Override
    public void initialize(URL url, ResourceBundle resourceBundle) {
        samplesData = (SamplesData)resourceBundle.getObject("");
        firstSample = samplesData.getFirstSample();
        secondSample = samplesData.getSecondSample();
        setChartPane();
        name.setCellValueFactory(new PropertyValueFactory<>("name"));
        value.setCellValueFactory(new PropertyValueFactory<>("value"));
        standardDeviation.setCellValueFactory(new PropertyValueFactory<>("standardDeviation"));
        interval.setCellValueFactory(new PropertyValueFactory<>("interval"));
        setTableTextWrap(table);
        setFirstSample();
        coefficientName.setCellValueFactory(new PropertyValueFactory<>("name"));
        coefficientValue.setCellValueFactory(new PropertyValueFactory<>("value"));
        coefficientInterval.setCellValueFactory(new PropertyValueFactory<>("interval"));
        coefficientTest.setCellValueFactory(new PropertyValueFactory<>("test"));
        coefficientQuantile.setCellValueFactory(new PropertyValueFactory<>("quantile"));
        coefficientSignificance.setCellValueFactory(new PropertyValueFactory<>("significance"));
        coefficientResult.setCellValueFactory(new PropertyValueFactory<>("result"));
        setTableTextWrap(coefficientTable);
        setTableTextWrap(ratioLineTable);
        setCorrelationCoefficients();
        linearRegressionSeries = new ArrayList<>();
        exponentialRegressionSeries = new ArrayList<>();
        LinearRegression linearRegression = new LinearRegression(firstSample, secondSample);
        ExponentialRegression exponentialRegression = new ExponentialRegression(firstSample, secondSample);
        setRegressionSeries(linearRegression, linearRegressionSeries);
        setRegressionSeries(exponentialRegression, exponentialRegressionSeries);
        linearRegressionSeries.get(0).setName(String.format("y = %.3f + %.3fx",
                linearRegression.getModel().intercept().getValue(),
                linearRegression.getModel().slope().getValue()));
        exponentialRegressionSeries.get(0).setName(String.format("y = %.3f * exp(%.3fx)",
                exponentialRegression.getA().getValue(),
                exponentialRegression.getB().getValue()));
        double yMin = secondSample.stream().mapToDouble(y -> y).min().orElse(0);
        if(yMin <= 0){
            expToggle.setDisable(true);
        }
    }

    private <T> void setTableTextWrap(TableView<T> tableView){
        for (TableColumn<T, ?> tableColumn : tableView.getColumns()) {
            TableColumn<T, String> textColumn = (TableColumn<T, String>) tableColumn;
            textColumn.setCellFactory(tc -> {
                TableCell<T, String> cell = new TableCell<>();
                Text text = new Text();
                cell.setGraphic(text);
                cell.setPrefHeight(Control.USE_COMPUTED_SIZE);
                text.wrappingWidthProperty().bind(textColumn.widthProperty());
                text.textProperty().bind(cell.itemProperty());
                return cell;
            });
            makeHeaderWrappable(tableColumn);
        }
    }

    private void makeHeaderWrappable(TableColumn<?, ?> col) {
        Label label = new Label(col.getText());
        label.setStyle("-fx-padding: 5px;");
        label.setWrapText(true);
        label.setAlignment(Pos.CENTER);
        label.setTextAlignment(TextAlignment.CENTER);

        StackPane stack = new StackPane();
        stack.setMaxHeight(50);
        label.setMaxHeight(50);
        stack.getChildren().add(label);
        stack.prefWidthProperty().bind(col.widthProperty().subtract(5));
        label.prefWidthProperty().bind(stack.prefWidthProperty());
        col.setText(null);
        col.setGraphic(stack);
    }

    private void setChartPane(){
        final NumberAxis xAxis = new NumberAxis();
        final NumberAxis yAxis = new NumberAxis();
        lineChart = new LineChart<>(xAxis, yAxis);
        lineChart.setLegendVisible(true);
        lineChart.setAnimated(false);
        lineChart.getStylesheets().addAll(Objects.requireNonNull(
                getClass().getClassLoader().getResource("style/correlation-chart.css")).toExternalForm());
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        for (int i = 0; i < firstSample.size(); i++) {
            double x = firstSample.get(i);
            double y = secondSample.get(i);
            series.getData().add(new XYChart.Data<>(x, y));
        }
        series.setName(samplesData.getSecondColumnName() + " (" + samplesData.getFirstColumnName() + ")");
        lineChart.getData().add(series);
        chartPane.getChildren().add(lineChart);
    }

    private void setCorrelationCoefficients(){
        if(firstSample.isEmpty()){
            return;
        }
        final int N = firstSample.size();
        double pearsonCoefficient = DataService.getPearsonCorrelationCoefficient(firstSample, secondSample);
        double spearmanCoefficient = DataService.getSpearmanRankCorrelationCoefficient(firstSample, secondSample);
        double kendallCoefficient = DataService.getKendallRankCorrelationCoefficient(firstSample, secondSample);
        double ratioCoefficient = DataService.getCorrelationRatio(firstSample, secondSample);
        double k = DataService.getNumberOfClasses(N);

        List<CorrelationCoefficient> list = new ArrayList<>();
        list.add(new CorrelationCoefficient("Пірсона", pearsonCoefficient,
                DataService.getPearsonCorrelationIntervalBorder(pearsonCoefficient, N, true),
                DataService.getPearsonCorrelationIntervalBorder(pearsonCoefficient, N, false),
                DataService.getPearsonCorrelationCoefficientTest(pearsonCoefficient, N),
                AnalysisFunctions.quantileNormalDistribution(DataService.P)));
        list.add(new CorrelationCoefficient("Спірмена", spearmanCoefficient,
                DataService.getSpearmanCorrelationCoefficientTest(spearmanCoefficient, N),
                AnalysisFunctions.quantileStudentsDistribution(DataService.P, N - 2)));
        list.add(new CorrelationCoefficient("Кендалла", kendallCoefficient,
                DataService.getKendallRankCorrelationCoefficientTest(kendallCoefficient, N),
                AnalysisFunctions.quantileNormalDistribution(DataService.P)));
        list.add(new CorrelationCoefficient("Кореляційне відношення", ratioCoefficient,
                DataService.getCorrelationRatioTest(ratioCoefficient, N),
                AnalysisFunctions.quantileFishersDistribution(1 - DataService.A, k - 1, N - k)));
        coefficientTable.setItems(FXCollections.observableArrayList(list));

        if(!DataService.checkCorrelationRatio(ratioCoefficient, N)) {
            ratioLineTable.setDisable(true);
            return;
        }

        double test = DataService.getCorrelationRatioLineTest(ratioCoefficient, pearsonCoefficient, N);
        double pearsonRatioCoefficient = DataService.getPearsonCorrelationRatioCoefficient(firstSample, secondSample);
        double quantile = AnalysisFunctions.quantileFishersDistribution(1 - DataService.A, k - 2, N - k);
        pearsonCoefficientValue.setCellValueFactory(o -> new SimpleStringProperty(String.format("%.3f", pearsonRatioCoefficient)));
        correlationRatioValue.setCellValueFactory(o -> new SimpleStringProperty(String.format("%.3f", ratioCoefficient)));
        ratioLineTest.setCellValueFactory(o -> new SimpleStringProperty(String.format("%.3f", test)));
        ratioLineQuantile.setCellValueFactory(o -> new SimpleStringProperty(String.format("%.3f", quantile)));
        ratioLineEquality.setCellValueFactory(o -> new SimpleStringProperty(test <= quantile ? "рівні" : "не рівні"));
        ratioLineResult.setCellValueFactory(o -> new SimpleStringProperty(test <= quantile ? "лінійний" : "нелінійний"));
        ratioLineTable.getItems().add(null);
    }

    private void setSampleCharacteristics(List<Double> sample){
        if(sample.isEmpty()){
            return;
        }
        table.getItems().clear();
        table.getItems().addAll(DataService.getStaticalCharacteristics(sample));
        double tA = DataService.getSkewnessCoefficient(sample) /
                DataService.getSkewnessCoefficientDeviation(sample);
        double tE = DataService.getKurtosisCoefficient(sample) /
                DataService.getKurtosisCoefficientDeviation(sample);
        double t = AnalysisFunctions.quantileStudentsDistribution(DataService.P, sample.size() - 1);
        skewnessLabel.setText(Math.abs(tA) <= t ? "A = 0" : "A <> 0");
        kurtosisLabel.setText(Math.abs(tE) <= t ? "E = 0" : "E <> 0");
        resultLabel.setText(DataService.isNormalDistributionByCoefficients(sample) ?
                "Ідентифіковано нормальний розподіл" : "Нормальний розподіл НЕ ідентифіковано");
    }

    public void setFirstSample() {
        setSampleCharacteristics(firstSample);
    }

    public void setSecondSample() {
        setSampleCharacteristics(secondSample);
    }

    private void setRegressionSeries(IRegression regression, List<XYChart.Series<Number, Number>> seriesList){
        XYChart.Series<Number, Number> func = new XYChart.Series<>();
        XYChart.Series<Number, Number> funcUpperBorders = new XYChart.Series<>();
        XYChart.Series<Number, Number> funcLowerBorders = new XYChart.Series<>();
        XYChart.Series<Number, Number> funcValueUpperBorders = new XYChart.Series<>();
        XYChart.Series<Number, Number> funcValueLowerBorders = new XYChart.Series<>();
        for (double x : firstSample) {
            func.getData().add(new XYChart.Data<>(x, regression.getFunctionValue(x)));
            Pair<Double, Double> borders = regression.getFunctionBorders(x);
            funcLowerBorders.getData().add(new XYChart.Data<>(x, borders.getKey()));
            funcUpperBorders.getData().add(new XYChart.Data<>(x, borders.getValue()));
            borders = regression.getFunctionValueBorders(x);
            funcValueLowerBorders.getData().add(new XYChart.Data<>(x, borders.getKey()));
            funcValueUpperBorders.getData().add(new XYChart.Data<>(x, borders.getValue()));
        }
        seriesList.addAll(List.of(
                func, funcLowerBorders, funcUpperBorders, funcValueLowerBorders, funcValueUpperBorders));

    }

    public void onLinearRegression() {
        if(linToggle.isSelected()){
            expToggle.setSelected(false);
            lineChart.getData().removeAll(exponentialRegressionSeries);
            lineChart.getData().addAll(linearRegressionSeries);
        }else {
            lineChart.getData().removeAll(linearRegressionSeries);
        }
    }

    public void onExponentialRegression() {
        if(expToggle.isSelected()){
            linToggle.setSelected(false);
            lineChart.getData().removeAll(linearRegressionSeries);
            lineChart.getData().addAll(exponentialRegressionSeries);
        }else {
            lineChart.getData().removeAll(exponentialRegressionSeries);
        }
    }
}
