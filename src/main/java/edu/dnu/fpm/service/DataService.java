package edu.dnu.fpm.service;

import edu.dnu.fpm.func.AnalysisFunctions;
import edu.dnu.fpm.model.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class DataService {
    public static final double P = 0.975;
    public static final double A = 0.05;

    public static double getFTest(
            List<Double> firstSample, List<Double> secondSample, Function<Double, Double> func){
        double ym = getMean(secondSample);
        double rv = getResidualVariance(firstSample, secondSample, func);
        return firstSample.stream()
                .mapToDouble(x -> Math.pow(func.apply(x) - ym, 2))
                .sum() / rv;
    }

    public static LinearRegressionModel getLinearRegression(List<Double> firstSample, List<Double> secondSample){
        final int N = firstSample.size();
        double xm = getMean(firstSample);
        double ym = getMean(secondSample);
        double xs = getStandardDeviation(firstSample, true);
        double ys = getStandardDeviation(secondSample, true);
        double r = getPearsonCorrelationCoefficient(firstSample, secondSample);
        double q = AnalysisFunctions.quantileStudentsDistribution(P, N - 2);
        double b = r * ys / xs;
        double a = ym - b * xm;
        double rv = getResidualVariance(firstSample, secondSample, x -> a + b * x);
        double as = Math.sqrt(rv * (1 + Math.pow(xm / xs, 2)) / N);
        double bs = Math.sqrt(rv / (N * xs * xs));

        return new LinearRegressionModel(new RegressionParameter(a, as, q),
                new RegressionParameter(b, bs, q), r, rv);
    }


    public static double getResidualVariance
            (List<Double> firstSample, List<Double> secondSample, Function<Double, Double> func){
        final int N = firstSample.size();
        return IntStream
                .range(0, N)
                .mapToDouble(i -> Math.pow(secondSample.get(i) - func.apply(firstSample.get(i)), 2))
                .sum() / (N - 2);
    }

    public static double getPearsonCorrelationCoefficient
            (List<Double> firstSample, List<Double> secondSample){
        double productSampleMean = IntStream
                .range(0, firstSample.size())
                .mapToObj(i -> firstSample.get(i) * secondSample.get(i))
                .collect(Collectors.averagingDouble(x -> x));
        double firstSampleMean = getMean(firstSample);
        double secondSampleMean = getMean(secondSample);
        double firstSampleDeviation = getStandardDeviation(firstSample, true);
        double secondSampleDeviation = getStandardDeviation(secondSample, true);
        return (productSampleMean - firstSampleMean * secondSampleMean) /
                (firstSampleDeviation * secondSampleDeviation);
    }

    public static double getPearsonCorrelationCoefficientTest(double value, int size){
        if(Math.abs(value) >= 1){
            return Double.POSITIVE_INFINITY;
        }
        return value * Math.sqrt(size - 2) / Math.sqrt(1 - value * value);
    }

    public static boolean checkPearsonCorrelation(double value, int size){
        return Math.abs(getPearsonCorrelationCoefficientTest(value, size)) >
                AnalysisFunctions.quantileNormalDistribution(P);
    }

    public static double getPearsonCorrelationIntervalBorder(double value, int size, boolean lower){
        double u = AnalysisFunctions.quantileNormalDistribution(P);
        double v = value + value * (1 - value * value) / (2 * size);
        double d = u * (1 - value * value) / Math.sqrt(size - 1);
        return lower ? v - d : v + d;
    }

    public static Map<Double, List<Double>> getCorrelationRatioFunction
            (List<Double> firstSample, List<Double> secondSample){
        final int N = firstSample.size();
        int k = getNumberOfClasses(N);
        double min = getMin(firstSample);
        double max = getMax(firstSample);
        double h = (max - min) / k;
        double eps = 0.0000001;

        return Stream
                .iterate(0, l -> l < k, l -> l + 1)
                .collect(Collectors.toMap(
                        l -> (min + h * l + h / 2),
                        l -> IntStream
                                .range(0, N)
                                .filter(i -> {
                                    double left = min + h * l;
                                    double right = left + h;
                                    double value = firstSample.get(i);
                                    return (value >= left && value < right) ||
                                            (Math.abs(right - max) < eps && Math.abs(value - max) < eps);
                                })
                                .mapToObj(secondSample::get)
                                .toList()));
    }

    public static double getCorrelationRatio
            (List<Double> firstSample, List<Double> secondSample){
        Collection<List<Double>> y = getCorrelationRatioFunction(firstSample, secondSample).values();
        double yMean = getMean(secondSample);
        double s1 = y.stream()
                .mapToDouble(l -> l.size() * Math.pow(getMean(l) - yMean, 2))
                .sum();
        double s2 = y.stream()
                .mapToDouble(l -> l.stream()
                        .mapToDouble(yl -> Math.pow(yl - yMean, 2))
                        .sum())
                .sum();

        return Math.sqrt(s1 / s2);
    }

    public static double getCorrelationRatioTest(double value, int size){
        if(value >= 1){
            return Double.POSITIVE_INFINITY;
        }
        int k = getNumberOfClasses(size);
        double p = Math.pow(value, 2);
        return (p / (k - 1)) / ((1 - p) / (size - k));
    }

    public static boolean checkCorrelationRatio(double value, int size){
        int k = getNumberOfClasses(size);
        return getCorrelationRatioTest(value, size) >
                AnalysisFunctions.quantileFishersDistribution(1 - A, k - 1, size - k);
    }

    public static double getPearsonCorrelationRatioCoefficient
            (List<Double> firstSample, List<Double> secondSample){
        Map<Double, List<Double>> map = getCorrelationRatioFunction(firstSample, secondSample);
        List<Double> x = new ArrayList<>();
        List<Double> y = new ArrayList<>();
        map.forEach((key, val) -> {
            for (double el : val) {
                x.add(key);
                y.add(el);
            }
        });
        return getPearsonCorrelationCoefficient(x, y);
    }

    public static double getCorrelationRatioLineTest(double value, double pearsonCorrelation, int size){
        if(value >= 1){
            return Double.POSITIVE_INFINITY;
        }
        int k = getNumberOfClasses(size);
        double p = Math.pow(value, 2);
        double r = Math.pow(pearsonCorrelation, 2);
        return ((p - r) / (k - 2)) / ((1 - p) / (size - k));
    }

    public static boolean checkCorrelationLineRatio(double value, double pearsonCorrelation, int size){
        int k = getNumberOfClasses(size);
        return getCorrelationRatioLineTest(value, pearsonCorrelation, size) >
                AnalysisFunctions.quantileFishersDistribution(1 - A, k - 2, size - k);
    }

    public static double getSpearmanRankCorrelationCoefficient
            (List<Double> firstSample, List<Double> secondSample){
        Map<Double, Double> firstSampleRanks = getSampleRanks(firstSample);
        Map<Double, Double> secondSampleRanks = getSampleRanks(secondSample);
        List<Double> xr = new ArrayList<>();
        List<Double> yr = new ArrayList<>();
        for (int i = 0; i < firstSample.size(); i++) {
            double x = firstSample.get(i);
            double y = secondSample.get(i);
            xr.add(firstSampleRanks.get(x));
            yr.add(secondSampleRanks.get(y));
        }
        return getPearsonCorrelationCoefficient(xr, yr);
    }

    public static double getSpearmanCorrelationCoefficientTest(double value, int size){
        if(Math.abs(value) >= 1){
            return Double.POSITIVE_INFINITY;
        }
        return value * Math.sqrt(size - 2) / Math.sqrt(1 - value * value);
    }

    public static boolean checkSpearmanCorrelationCoefficient(double value, int size){
        return Math.abs(getSpearmanCorrelationCoefficientTest(value, size)) >
                AnalysisFunctions.quantileStudentsDistribution(P, size - 2);
    }

    public static double getKendallRankCorrelationCoefficient
            (List<Double> firstSample, List<Double> secondSample){
        final int N = firstSample.size();
        Map<Double, Rank> firstSampleRanks = getSampleElementRanks(firstSample);
        Map<Double, Rank> secondSampleRanks = getSampleElementRanks(secondSample);

        double s = 0;
        boolean contains = false;
        for (int i = 0; i < N - 1; i++) {
            double rxi = firstSampleRanks.get(firstSample.get(i)).value();
            double ryi = secondSampleRanks.get(secondSample.get(i)).value();
            for (int j = i + 1; j < N; j++) {
                double rxj = firstSampleRanks.get(firstSample.get(j)).value();
                double ryj = secondSampleRanks.get(secondSample.get(j)).value();
                if(Double.compare(rxi, rxj) == 0 || Double.compare(ryi, ryj) == 0){
                    contains = true;
                    continue;
                }
                s += (rxi < rxj && ryi < ryj) || (rxi > rxj && ryi > ryj) ? 1 : -1;
            }
        }
        if(!contains){
            return 2 * s / (N * (N - 1));
        }

        double c = firstSampleRanks.values().stream()
                .mapToInt(r -> r.count * (r.count - 1))
                .sum() / 2.0;
        double d = secondSampleRanks.values().stream()
                .mapToInt(r -> r.count * (r.count - 1))
                .sum() / 2.0;

        return s / Math.sqrt((0.5 * N * (N - 1) - c) * (0.5 * N * (N - 1) - d));
    }

    public static double getKendallRankCorrelationCoefficientTest(double value, int size){
        return value * Math.sqrt(9 * size * (size - 1)) / Math.sqrt(2 * (2 * size + 5));
    }

    public static boolean checkKendallCorrelationCoefficient(double value, int size) {
        return Math.abs(getPearsonCorrelationCoefficientTest(value, size)) >
                AnalysisFunctions.quantileNormalDistribution(P);
    }

    private record Rank(double value, int count) {}

    public static Map<Double, Rank> getSampleElementRanks(List<Double> sample){
        Double[] data = getSortedArray(sample);
        Map<Double, Rank> map = new HashMap<>();
        int count = 0;
        for (int i = 0; i < sample.size() - 1; i++) {
            count++;
            if(Double.compare(data[i], data[i + 1]) < 0){
                double rank = (2 * i - count + 3) / 2.0;
                map.put(data[i], new Rank(rank, count));
                count = 0;
            }
        }
        double rank = (2 * sample.size() - count) / 2.0;
        map.put(data[sample.size() - 1], new Rank(rank, count + 1));

        return map;
    }

    public static Map<Double, Double> getSampleRanks(List<Double> sample){
        Double[] data = getSortedArray(sample);
        Map<Double, Double> map = new HashMap<>();
        int count = 0;
        for (int i = 0; i < sample.size() - 1; i++) {
            count++;
            if(Double.compare(data[i], data[i + 1]) < 0){
                double rank = (2 * i - count + 3) / 2.0;
                map.put(data[i], rank);
                count = 0;
            }
        }
        double rank = (2 * sample.size() - count) / 2.0;
        map.put(data[sample.size() - 1], rank);

        return map;
    }

    public static double getSignedRankTest(List<Double> firstSample, List<Double> secondSample){
        List<Double> diffSample = getDiffSample(firstSample, secondSample);
        diffSample.removeIf(x -> x == 0);
        if(diffSample.isEmpty()){
            return 0;
        }
        Map<Double, Double> ranks = getSampleRanks(
                diffSample.stream()
                        .map(Math::abs)
                        .collect(Collectors.toList()));
        double t = 0;
        for (Double x : diffSample) {
            int a = x > 0 ? 1 : 0;
            t += a * ranks.get(Math.abs(x));
        }
        final int N = diffSample.size();
        final double E = N * (N + 1) / 4.0;
        final double D = N * (N + 1) * (2 * N + 1) / 24.0;

        return (t - E) / Math.sqrt(D);
    }


    public static double getRankSumTest(List<Double> firstSample, List<Double> secondSample){
        Map<Double, Double> ranks = getSampleRanks(
                getMergeSample(firstSample, secondSample));
        double w = 0;
        for (Double x : firstSample) {
            w += ranks.get(x);
        }
        final int N1 = firstSample.size();
        final int N2 = secondSample.size();
        final double E = N1 * (N1 + N2 + 1) / 2.0;
        final double D = N1 * N2 * (N1 + N2 + 1) / 12.0;

        return (w - E) / Math.sqrt(D);
    }

    public static List<Double> getMergeSample(List<Double> firstSample, List<Double> secondSample){
        List<Double> mergeSample = new ArrayList<>();
        mergeSample.addAll(firstSample);
        mergeSample.addAll(secondSample);
        return mergeSample;
    }

    public static List<Double> getDiffSample(List<Double> firstSample, List<Double> secondSample){
        List<Double> diffSample = new ArrayList<>();
        for (int i = 0; i < firstSample.size(); i++) {
            diffSample.add(firstSample.get(i) - secondSample.get(i));
        }
        return diffSample;
    }

    public static double getMeanTest(List<Double> firstSample, List<Double> secondSample, boolean depend){
        if(depend){
            List<Double> diffSample = getDiffSample(firstSample, secondSample);
            double m = getMean(diffSample);
            double s = getStandardDeviation(diffSample, false);
            if(s == 0){
                return 0;
            }
            return m * Math.sqrt(diffSample.size()) / s;
        }
        final int N1 = firstSample.size();
        final int N2 = secondSample.size();
        double d1 = Math.pow(getStandardDeviation(firstSample, false), 2);
        double d2 = Math.pow(getStandardDeviation(secondSample, false), 2);
        double m1 = getMean(firstSample);
        double m2 = getMean(secondSample);
        double d = ((N1 - 1) * d1 + (N2 - 1) * d2) / (N1 + N2 - 2);
        return checkVarianceEquals(firstSample, secondSample) ?
                (m1 - m2) / Math.sqrt(d / N1 + d / N2) :
                (m1 - m2) / Math.sqrt(d1 / N1 + d2 / N2);
    }

    public static double getMeanTestQuantile
            (List<Double> firstSample, List<Double> secondSample, boolean depend){
        double v = getMeanFreedomDegree(firstSample, secondSample, depend);
        return AnalysisFunctions.quantileStudentsDistribution(P, v);
    }

    public static double getMeanFreedomDegree
            (List<Double> firstSample, List<Double> secondSample, boolean depend){
        if(depend){
            return firstSample.size() - 1;
        }
        final int N1 = firstSample.size();
        final int N2 = secondSample.size();
        double d1 = Math.pow(getStandardDeviation(firstSample, false), 2);
        double d2 = Math.pow(getStandardDeviation(secondSample, false), 2);
        return checkVarianceEquals(firstSample, secondSample) ?
                (N1 + N2 - 2) : Math.pow(d1 / N1 + d2 / N2, 2) /
                (Math.pow(d1 / N1, 2) / (N1 - 1) + Math.pow(d2 / N2, 2) / (N2 - 1));
    }

    public static double getMeanPValue(List<Double> firstSample, List<Double> secondSample, boolean depend){
        double t = getMeanTest(firstSample, secondSample, depend);
        double v = getMeanFreedomDegree(firstSample, secondSample, depend);
        double studentsQuantile = AnalysisFunctions.CDF_S(Math.abs(t), v);

        return 2 * (1 - studentsQuantile);
    }

    public static boolean checkMeanEquals(List<Double> firstSample, List<Double> secondSample, boolean depend){
        double meanTest = getMeanTest(firstSample, secondSample, depend);
        double studentsQuantile = getMeanTestQuantile(firstSample,secondSample, depend);
        return Math.abs(meanTest) <= studentsQuantile;
    }

    public static double getVarianceTest(List<Double> firstSample, List<Double> secondSample){
        double firstVariance = Math.pow(getStandardDeviation(firstSample, false), 2);
        double secondVariance = Math.pow(getStandardDeviation(secondSample, false), 2);
        return firstVariance >= secondVariance ?
                firstVariance / secondVariance : secondVariance / firstVariance;
    }

    public static double getVarianceTestQuantile
            (List<Double> firstSample, List<Double> secondSample){
        double firstVariance = Math.pow(getStandardDeviation(firstSample, false), 2);
        double secondVariance = Math.pow(getStandardDeviation(secondSample, false), 2);
        double v1, v2;
        if(firstVariance >= secondVariance){
            v1 = firstSample.size() - 1;
            v2 = secondSample.size() - 1;
        } else {
            v1 = secondSample.size() - 1;
            v2 = firstSample.size() - 1;
        }
        return AnalysisFunctions.quantileFishersDistribution(1 - A, v1, v2);
    }

    public static double getVariancePValue(List<Double> firstSample, List<Double> secondSample){
        final int N1 = firstSample.size();
        final int N2 = secondSample.size();
        double firstVariance = Math.pow(getStandardDeviation(firstSample, false), 2);
        double secondVariance = Math.pow(getStandardDeviation(secondSample, false), 2);
        double f = firstVariance / secondVariance;
        double fishersQuantile = AnalysisFunctions.CDF_F(f, N1 - 1, N2 - 1);
        return f > 1 ? 2 * (1 - fishersQuantile) : 2 * fishersQuantile;
    }

    public static boolean checkVarianceEquals(List<Double> firstSample, List<Double> secondSample){
        double varianceTest = getVarianceTest(firstSample, secondSample);
        double fishersQuantile = getVarianceTestQuantile(firstSample, secondSample);
        return varianceTest <= fishersQuantile;
    }

    public static double getKurtosisCoefficient(List<Double> sample){
        final int N = sample.size();
        double sum = getDiffSum(sample, 4);
        double s = getStandardDeviation(sample,true);
        double biasCoefficient = sum / (N * Math.pow(s, 4)) - 3;
        if(s == 0){
            biasCoefficient = -3;
        }
        return (N * N - 1) * (biasCoefficient + 6.0 / (N + 1)) / ((N - 2) * (N - 3));
    }

    public static double getKurtosisCoefficientDeviation(List<Double> sample){
        final double N = sample.size();
        return Math.sqrt((24 * N * (N - 1) * (N - 1)) / ((N - 3) * (N - 2) * (N + 3) * (N + 5)));
    }

    public static double getSkewnessCoefficient(List<Double> sample){
        final int N = sample.size();
        double sum = getDiffSum(sample, 3);
        double s = getStandardDeviation(sample, true);
        double biasCoefficient = sum / (N * Math.pow(s, 3));
        if(s == 0){
            return 0;
        }
        return Math.sqrt(N * (N - 1)) * biasCoefficient / (N - 2);
    }

    public static double getSkewnessCoefficientDeviation(List<Double> sample){
        final double N = sample.size();
        return Math.sqrt((6 * N * (N - 1)) / ((N - 2) * (N + 1) * (N + 3)));
    }

    private static Double[] getSortedArray(List<Double> sample){
        Double[] tmp = sample.toArray(new Double[0]);
        Arrays.sort(tmp);
        return tmp;
    }

    public static double getMedian(List<Double> sample){
        Double[] data = getSortedArray(sample);
        int index = sample.size() / 2;
        return sample.size() % 2 == 0 ? (data[index - 1] + data[index]) / 2 : data[index];
    }

    public static double getMean(List<Double> sample) {
        return sample.stream().collect(Collectors.averagingDouble(x -> x));
    }

    public static double getMeanDeviation(List<Double> sample){
        double s = getStandardDeviation(sample, false);
        return s / Math.sqrt(sample.size());
    }

    public static double getStandardDeviation(List<Double> sample, boolean bias){
        final int N = sample.size();
        double div = bias ? N : N - 1;
        return Math.sqrt(getDiffSum(sample, 2) / div);
    }

    public static double getStandardDeviationDeviation(List<Double> sample){
        return getStandardDeviation(sample,false) / Math.sqrt(2 * sample.size());
    }

    public static double getMedianIntervalBorder(List<Double> sample, boolean isLeft){
        final int N = sample.size();
        double p = AnalysisFunctions.quantileNormalDistribution(P) * Math.sqrt(N) / 2;
        int index = (int) (isLeft ? N / 2 - p : N / 2 + p);
        return getSortedArray(sample)[index];
    }

    private static double getDiffSum(List<Double> sample, int p){
        double m = getMean(sample);
        double sum = 0;
        for (double x : sample) {
            sum += Math.pow(x - m, p);
        }
        return sum;
    }

    public static double getMin(List<Double> sample){
        return sample.stream()
                .min(Comparator.comparingDouble(x -> x))
                .orElse(0.0);
    }

    public static double getMax(List<Double> sample){
        return sample.stream()
                .max(Comparator.comparingDouble(x -> x))
                .orElse(0.0);
    }

    public static boolean isNormalDistributionByCoefficients(List<Double> sample){
        double tA = getSkewnessCoefficient(sample) / getSkewnessCoefficientDeviation(sample);
        double tE = getKurtosisCoefficient(sample) / getKurtosisCoefficientDeviation(sample);
        double t = AnalysisFunctions.quantileStudentsDistribution(P, sample.size() - 1);
        return Math.abs(tA) <= t && Math.abs(tE) <= t;
    }

    public static int getNumberOfClasses(int size){
        return (int)Math.ceil(1 + 1.44 * Math.log(size));
    }

    public static List<CorrelationCoefficient> getCorrelationCoefficients
            (List<Double> firstSample, List<Double> secondSample){
        final int N = firstSample.size();
        double pearsonCoefficient = getPearsonCorrelationCoefficient(firstSample, secondSample);
        double spearmanCoefficient = getSpearmanRankCorrelationCoefficient(firstSample, secondSample);
        double kendallCoefficient = getKendallRankCorrelationCoefficient(firstSample, secondSample);
        double ratioCoefficient = getCorrelationRatio(firstSample, secondSample);
        double k = getNumberOfClasses(N);

        List<CorrelationCoefficient> list = new ArrayList<>();
        list.add(new CorrelationCoefficient("Пірсона", pearsonCoefficient,
                getPearsonCorrelationIntervalBorder(pearsonCoefficient, N, true),
                getPearsonCorrelationIntervalBorder(pearsonCoefficient, N, false),
                getPearsonCorrelationCoefficientTest(pearsonCoefficient, N),
                AnalysisFunctions.quantileNormalDistribution(P)));
        list.add(new CorrelationCoefficient("Спірмена", spearmanCoefficient,
                getSpearmanCorrelationCoefficientTest(spearmanCoefficient, N),
                AnalysisFunctions.quantileStudentsDistribution(P, N - 2)));
        list.add(new CorrelationCoefficient("Кендалла", kendallCoefficient,
                getKendallRankCorrelationCoefficientTest(kendallCoefficient, N),
                AnalysisFunctions.quantileNormalDistribution(P)));
        list.add(new CorrelationCoefficient("Кореляційне відношення", ratioCoefficient,
                getCorrelationRatioTest(ratioCoefficient, N),
                AnalysisFunctions.quantileFishersDistribution(1 - A, k - 1, N - k)));
        return list;
    }

    public static List<StatisticalCharacteristic> getStaticalCharacteristics(List<Double> sample){
        List<StatisticalCharacteristic> list = new ArrayList<>();

        list.add(getStatisticalCharacteristic(
                "Середнє арифметичне", getMean(sample), getMeanDeviation(sample), sample.size())
        );
        list.add(new StatisticalCharacteristic(
                "Медіана", getMedian(sample), null,
                getMedianIntervalBorder(sample,true), getMedianIntervalBorder(sample, false))
        );
        list.add(getStatisticalCharacteristic(
                "Середньоквадратичне відхилення",
                getStandardDeviation(sample, false), getStandardDeviationDeviation(sample), sample.size())
        );
        list.add(getStatisticalCharacteristic(
                "Коефіцієнт асиметрії",
                getSkewnessCoefficient(sample), getSkewnessCoefficientDeviation(sample), sample.size())
        );
        list.add(getStatisticalCharacteristic(
                "Коефіцієнт ексцесу",
                getKurtosisCoefficient(sample), getKurtosisCoefficientDeviation(sample), sample.size())
        );
        list.add(new StatisticalCharacteristic(
                "Мінімум", getMin(sample), null, null, null)
        );
        list.add(new StatisticalCharacteristic(
                "Максимум", getMax(sample), null, null, null)
        );
        return list;
    }

    private static StatisticalCharacteristic getStatisticalCharacteristic(
            String name, Double value, Double deviation, int size){
        double p = AnalysisFunctions.quantileStudentsDistribution(P, size - 1) * deviation;
        return new StatisticalCharacteristic(
                name, value, deviation,
                value - p, value + p
        );
    }

}
