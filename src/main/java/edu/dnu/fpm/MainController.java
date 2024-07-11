package edu.dnu.fpm;


import edu.dnu.fpm.model.SamplesData;
import edu.dnu.fpm.service.DataService;
import javafx.beans.property.SimpleStringProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.fxml.FXMLLoader;
import javafx.fxml.Initializable;
import javafx.scene.Parent;
import javafx.scene.control.*;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.URL;
import java.util.*;

public class MainController implements Initializable {
    public ScrollPane content;
    public TableView<List<Double>> samplesTable;
    public TableColumn<List<Double>, String> firstSampleColumn;
    public TableColumn<List<Double>, String> secondSampleColumn;
    public TableView<List<Double>> table;
    public TableView<List<String>> correlationPMatrix;
    public TableView<List<String>> correlationRMatrix;
    public ComboBox<String> firstSampleBox;
    public ComboBox<String> secondSampleBox;
    public VBox parameters;
    private Parent correlation;
    private Parent regression;
    private SamplesData samplesData;
    private List<Double> firstSample;
    private List<Double> secondSample;


    public void loadData() {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Обрати файл");
        fileChooser.getExtensionFilters().add(
                new FileChooser.ExtensionFilter("Text Files", "*.txt","*.DAT"));
        File file = fileChooser.showOpenDialog(null);
        if(file != null){
            try {
                loadDataFromFileTable(file);
            } catch (Exception e) {
                new Alert(Alert.AlertType.ERROR, "Помилка при зчитуванні з файлу: " + file.getName())
                        .showAndWait();
            }
        }
    }

    private Map<String, List<Double>> map;

    private Map<String, List<Double>> getTableData(File file) throws FileNotFoundException {
        Map<String, List<Double>> newMap = new LinkedHashMap<>();
        Scanner fileScanner = new Scanner(file);
        String columns = fileScanner.nextLine();
        List<String> columnNames = new ArrayList<>();
        Scanner lineScanner = new Scanner(columns);
        while (lineScanner.hasNext()){
            columnNames.add(lineScanner.next());
        }
        for (String name : columnNames) {
            newMap.put(name, new ArrayList<>());
        }
        while (fileScanner.hasNextLine()){
            lineScanner = new Scanner(fileScanner.nextLine());
            int i = 0;
            while (lineScanner.hasNext()){
                double value = Double.parseDouble(lineScanner.next().replace(',', '.'));
                newMap.get(columnNames.get(i++)).add(value);
            }
        }
        return newMap;
    }

    private List<List<Double>> getTableRowsList(List<List<Double>> columns){
        if(columns.size() < 2){
            return null;
        }
        List<List<Double>> rows = new ArrayList<>();
        int rowSize = columns.get(0).size();
        for (int i = 0; i < rowSize; i++) {
            List<Double> row = new ArrayList<>();
            for (List<Double> column : columns) {
                row.add(column.get(i));
            }
            rows.add(row);
        }
        return rows;
    }

    private void setTables(List<String> columnNames){
        int size = columnNames.size();
        firstSampleBox.setItems(FXCollections.observableArrayList(columnNames));
        secondSampleBox.setItems(FXCollections.observableArrayList(columnNames));
        table.getColumns().clear();
        correlationPMatrix.getColumns().clear();
        correlationRMatrix.getColumns().clear();

        TableColumn<List<String>, String> matrixPName = new TableColumn<>();
        matrixPName.setCellValueFactory(o -> new SimpleStringProperty(o.getValue().get(0)));
        correlationPMatrix.getColumns().add(matrixPName);
        TableColumn<List<String>, String> matrixRName = new TableColumn<>();
        matrixRName.setCellValueFactory(o -> new SimpleStringProperty(o.getValue().get(0)));
        correlationRMatrix.getColumns().add(matrixRName);

        for (int i = 0; i < size; i++){
            int index = i;
            TableColumn<List<Double>, String> tableColumn = new TableColumn<>(columnNames.get(i));
            tableColumn.setCellValueFactory(o -> new SimpleStringProperty(o.getValue().get(index).toString()));
            table.getColumns().add(tableColumn);
            TableColumn<List<String>, String> matrixPColumn = new TableColumn<>(columnNames.get(i));
            matrixPColumn.setCellValueFactory(o -> new SimpleStringProperty(o.getValue().get(index + 1)));
            correlationPMatrix.getColumns().add(matrixPColumn);
            TableColumn<List<String>, String> matrixRColumn = new TableColumn<>(columnNames.get(i));
            matrixRColumn.setCellValueFactory(o -> new SimpleStringProperty(o.getValue().get(index + 1)));
            correlationRMatrix.getColumns().add(matrixRColumn);
        }
    }

    private void loadDataFromFileTable(File file) throws FileNotFoundException {
        map = getTableData(file);
        List<String> columnNames = map.keySet().stream().toList();
        int size = columnNames.size();
        setTables(columnNames);
        List<List<Double>> columns = map.values().stream().toList();
        ObservableList<List<Double>> tableList = FXCollections.observableArrayList(
                getTableRowsList(columns));
        table.setItems(tableList);
        ObservableList<List<String>> matrixPList = FXCollections.observableArrayList();
        ObservableList<List<String>> matrixRList = FXCollections.observableArrayList();
        for (int i = 0; i < size; i++) {
            List<String> matrixPValues = new ArrayList<>();
            List<String> matrixRValues = new ArrayList<>();
            matrixPValues.add(columnNames.get(i));
            matrixRValues.add(columnNames.get(i));
            for (int j = 0; j < size; j++) {
                double p = DataService.getPearsonCorrelationCoefficient(columns.get(i), columns.get(j));
                double r = DataService.getCorrelationRatio(columns.get(i), columns.get(j));
                matrixPValues.add(String.format("%.3f", p));
                matrixRValues.add(String.format("%.3f", r));
            }
            matrixPList.add(matrixPValues);
            matrixRList.add(matrixRValues);
        }
        correlationPMatrix.setItems(matrixPList);
        correlationRMatrix.setItems(matrixRList);
        samplesTable.getItems().clear();
        firstSampleColumn.setText("X");
        secondSampleColumn.setText("Y");
        content.setContent(parameters);
        firstSample.clear();
        secondSample.clear();
        correlation = null;
        regression = null;
    }

    @Override
    public void initialize(URL url, ResourceBundle resourceBundle) {
        firstSample = new ArrayList<>();
        secondSample = new ArrayList<>();
        firstSampleColumn.setCellValueFactory(o -> new SimpleStringProperty(o.getValue().get(0).toString()));
        secondSampleColumn.setCellValueFactory(o -> new SimpleStringProperty(o.getValue().get(1).toString()));
    }

    private Parent getNode(String name){
        FXMLLoader fxmlLoader = new FXMLLoader(MainController.class.getClassLoader().getResource(name));
        Parent root = null;
        try {
            fxmlLoader.setResources(new ResourceBundle() {
                @Override
                protected Object handleGetObject(String key) {
                    return samplesData;
                }
                @Override
                public Enumeration<String> getKeys() {
                    return null;
                }
            });
            root = fxmlLoader.load();
        }catch (Exception e){
            e.printStackTrace();
        }
        return root;
    }

    public void onCorrelation() {
        if(correlation == null){
            correlation = getNode("correlation.fxml");
        }
        content.setContent(correlation);
    }

    public void onRegression() {
        if(regression == null){
            regression = getNode("regression.fxml");
        }
        content.setContent(regression);
    }

    public void onParameters() {
        content.setContent(parameters);
    }

    public void setSamples() {
        String firstColumn = firstSampleBox.getValue();
        String secondColumn = secondSampleBox.getValue();
        if(firstColumn == null || secondColumn == null){
            return;
        }
        firstSample = map.get(firstColumn);
        secondSample = map.get(secondColumn);
        samplesData = new SamplesData(firstColumn, secondColumn,
                firstSample, secondSample);
        firstSampleColumn.setText(firstColumn);
        secondSampleColumn.setText(secondColumn);
        samplesTable.setItems(FXCollections.observableArrayList(
                getTableRowsList(List.of(firstSample, secondSample))));
        correlation = null;
        regression = null;
    }

}
