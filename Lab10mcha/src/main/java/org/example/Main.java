import org.apache.commons.math3.analysis.DifferentiableUnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.AdamsBashforthIntegrator;

public class Main {
    private static final double a = 0.9;
    private static final double m = 2.0;

    public static void main(String[] args) {
        double a = 0.9;
        double m = 2.0;
        double[] y0 = { 0.0 }; // начальные условия
        double tStart = 0.0;
        double tEnd = 1.0;

        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 1;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                yDot[0] = a * (1 - y[0] * y[0]) / ((1 + m) * t * t + y[0] * y[0] + 1);
            }
        };

        AdamsBashforthIntegrator integrator = new AdamsBashforthIntegrator(2, tStart, tEnd, 1.0e-10, 1.0e-6);
        integrator.integrate(ode, tStart, y0, tEnd, y0);

        System.out.println("Результат с использованием библиотеки Apache commons math: " + y0[0]);
        double[] y = solveAdams(0, 1, 0.1);
        printSolution(y);
    }

    private static double[] solveAdams(double start, double end, double step) {
        int n = (int) ((end - start) / step) + 1;
        double[] x = new double[n];
        double[] y = new double[n];

        // Используем метод Рунге-Кутта для получения первых четырех точек
        RungeKuttaMethod rungeKutta = new RungeKuttaMethod();
        x[0] = start;
        y[0] = 0.0;

        for (int i = 1; i <= 3; i++) {
            x[i] = x[i - 1] + step;
            y[i] = rungeKutta.step(x[i - 1], y[i - 1], step);
        }

        // Применяем метод Адамса
        for (int i = 4; i < n; i++) {
            x[i] = x[i - 1] + step;
            double predictor = y[i - 1] + step / 24.0 * (55 * function(x[i - 1], y[i - 1]) -
                    59 * function(x[i - 2], y[i - 2]) + 37 * function(x[i - 3], y[i - 3]) -
                    9 * function(x[i - 4], y[i - 4]));

            y[i] = y[i - 1] + step / 24.0 * (9 * function(x[i], predictor) +
                    19 * function(x[i - 1], y[i - 1]) - 5 * function(x[i - 2], y[i - 2]) +
                    function(x[i - 3], y[i - 3]));
        }

        return y;
    }

    private static double function(double x, double y) {
        return a * (1 - y * y) / ((1 + m) * x * x + y * y + 1);
    }

    private static void printSolution(double[] y) {
        System.out.println("Результат методом Адамса:");
        for (int i = 0; i < y.length; i++) {
            if (i == 10)
            System.out.println("x = " + i * 0.1 + ", y = " + y[i]);
        }
    }

    private static class RungeKuttaMethod {
        public double step(double x, double y, double h) {
            double k1 = h * function(x, y);
            double k2 = h * function(x + h / 2, y + k1 / 2);
            double k3 = h * function(x + h / 2, y + k2 / 2);
            double k4 = h * function(x + h, y + k3);

            return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        }
    }
}
