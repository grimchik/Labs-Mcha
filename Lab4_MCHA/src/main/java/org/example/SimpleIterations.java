package org.example;

import java.lang.Math;

public class SimpleIterations {

    private static final double EPSILON = 1e-6;
    private static final int MAX_ITERATIONS = 1000;

    public static double function1(double x, double y) {
        return Math.atan(x * y + 0.1);
    }

    public static double function2(double x) {
        return Math.sqrt((1 - 0.7 * x * x) / 2);
    }

    public static double[] simpleIterationMethod(double x0, double y0) {
        double x = x0;
        double y = y0;
        int iterations = 0;

        while (iterations < MAX_ITERATIONS) {
            double newX = function1(x, y);
            double newY = function2(x);

            if (Math.abs(newX - x) < EPSILON && Math.abs(newY - y) < EPSILON) {
                break;
            }

            x = newX;
            y = newY;
            iterations++;
        }

        return new double[]{x, y};
    }
    public static double f1(double x, double y) {
        return Math.tan(x * y + 0.1) - x;
    }

    public static double f2(double x, double y) {
        return 0.7 * x * x + 2 * y * y - 1;
    }

    public static double df1dx(double x, double y) {
        return y / (Math.cos(x * y + 0.1) * Math.cos(x * y + 0.1)) - 1;
    }

    public static double df1dy(double x, double y) {
        return x / (Math.cos(x * y + 0.1) * Math.cos(x * y + 0.1));
    }

    public static double df2dx(double x) {
        return 1.4 * x;
    }

    public static double df2dy(double y) {
        return 4 * y;
    }

    public static double[] newtonMethod(double x0, double y0) {
        double x = x0;
        double y = y0;

        for (int i = 0; i < MAX_ITERATIONS; i++) {
            double jacobianDet = df1dx(x, y) * df2dy(y) - df1dy(x, y) * df2dx(x);
            if (Math.abs(jacobianDet) < EPSILON) {
                throw new RuntimeException("Jacobian determinant is too small.");
            }

            double deltaX = (-f1(x, y) * df2dy(y) + f2(x, y) * df1dy(x, y)) / jacobianDet;
            double deltaY = (f1(x, y) * df2dx(x) - f2(x, y) * df1dx(x, y)) / jacobianDet;

            x += deltaX;
            y += deltaY;

            if (Math.abs(deltaX) < EPSILON && Math.abs(deltaY) < EPSILON) {
                break;
            }
        }

        return new double[]{x, y};
    }
    public static void main(String[] args) {
        double x0 = 1.0;
        double y0 = 1.0;

        double[] result1 = simpleIterationMethod(x0, y0);
        double[] result2 = newtonMethod(x0, y0);
        System.out.println("Solution by Simple Iterations:");
        System.out.println("x = " + result2[0]);
        System.out.println("y = " + result2[1]);// начальное приближение


        System.out.println("Solution by Newton's method:");
        System.out.println("x = " + result2[0]);
        System.out.println("y = " + result2[1]);
    }
}
