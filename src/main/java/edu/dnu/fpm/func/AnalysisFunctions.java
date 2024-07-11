package edu.dnu.fpm.func;

import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;

public class AnalysisFunctions {
    public static double kernelGaussianFunction(double u){
        return 1 / Math.sqrt(2 * Math.PI) * Math.exp(-u * u / 2);
    }

    public static double kolmogorovDistributionFunction(double z){
        double sum = 0;
        for (int k = 1; k < 6; k++) {
            double e = Math.exp(-2 * k * k * z * z);
            double p = k % 2 == 0 ? e : -e;
            sum += p;
        }
        return 1 + 2 * sum;
    }

    public static double normalDistributionEstimationFunction(double u){
        double r = 0.2316419;
        double t = 1 / (1 + r * u);
        double b1 = 0.31938153;
        double b2 = -0.356563782;
        double b3 = 1.781477937;
        double b4 = -1.821255978;
        double b5 = 1.330274429;
        double g = kernelGaussianFunction(u);
        double t2 = t * t;
        double t3 = t2 * t;
        double t4 = t3 * t;
        double t5 = t4 * t;

        return 1 - g * (b1 * t + b2 * t2 + b3 * t3 + b4 * t4 + b5 * t5);
    }

    public static double quantileFishersDistributionCalc(double p, double v1, double v2){
        double u = quantileNormalDistribution(p);
        double s = 1 / v1 + 1 / v2;
        double d = 1 / v1 - 1 / v2;

        double z1 = u * Math.sqrt(s / 2) - d * (u * u + 2) / 6;
        double z2 = Math.sqrt(s / 2) * (s * (u * u + 3 * u) / 24 + d * d * (u * u * u + 11 * u) / (72 * s));
        double z3 = d * s * (Math.pow(u, 4) + 9 * u * u + 8) / 120;
        double z4 = d * d * d * (3 * Math.pow(u, 4) + 7 * u * u - 16) / (3240 * s);
        double z5 = s * s * (Math.pow(u, 5) + 20 * u * u * u + 15 * u) / 1920;
        double z6 = Math.pow(d, 4) * (Math.pow(u, 5) + 44 * u * u * u + 183 * u) / 2880;
        double z7 = Math.pow(d, 4) * (9 * Math.pow(u, 5) - 284 * u * u * u - 1513 * u) / (155_520 * s * s);
        double z = z1 + z2 - z3 + z4 + Math.sqrt(s / 2) * (z5 + z6 + z7);

        return Math.exp(2 * z);
    }

    public static double quantileFishersDistribution(double p, double v1, double v2){
        if(v1 <= 0 || v2 <= 0){
            return -1;
        }
        return new FDistribution(v1, v2).inverseCumulativeProbability(p);
    }

    public static double CDF_F(double x, double v1, double v2){
        return new FDistribution(v1, v2).cumulativeProbability(x);
    }

    public static double CDF_S(double x, double v){
        return new TDistribution(v).cumulativeProbability(x);
    }

    public static double CDF_N(double x){
        return new NormalDistribution().cumulativeProbability(x);
    }


    public static double quantileStudentsDistributionCalc(double p, double v){
        double u = quantileNormalDistribution(p);
        double u2 = u * u;
        double u3 = u * u2;
        double u5 = u3 * u2;
        double u7 = u5 * u2;
        double u9 = u7 * u2;

        double g1 = (u3 + u) / 4;
        double g2 = (5 * u5 + 16 * u3 + 3 * u) / 96;
        double g3 = (3 * u7 + 19 * u5 + 17 * u3 - 15 * u) / 384;
        double g4 = (79 * u9 + 779 * u7 + 1482 * u5 - 1920 * u3 - 945 * u) / 92160;

        return u + g1 / v + g2 / (v * v) + g3 / (v * v * v) + g4 / (v * v * v * v);
    }

    public static double quantileStudentsDistribution(double p, double v){
        if(v <= 0){
            return -1;
        }
        return new TDistribution(v).inverseCumulativeProbability(p);
    }

    public static double quantileNormalDistributionCalc(double p){
        return (p > 0.5) ? quantileNormalDistributionFunction(1 - p) :
                -quantileNormalDistributionFunction(p);
    }

    public static double quantileNormalDistribution(double p){
        return new NormalDistribution().inverseCumulativeProbability(p);
    }

    public static double quantileNormalDistributionFunction(double a){
        double t = Math.sqrt(-2 * Math.log(a));
        double c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
        double d1 = 1.432788, d2 = 0.1892659, d3 = 0.001308;
        return t - (c0 + c1 * t + c2 * t * t) /
                (1 + d1 * t + d2 * t * t + d3 * t * t * t);
    }
}

