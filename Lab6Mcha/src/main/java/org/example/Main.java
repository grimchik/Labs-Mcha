import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.analysis.solvers.PolynomialSolver;

public class Main {
    public static void main(String[] args) {
        // Заданные значения
        int k = 3;
        double m = 1.5;

        //double[] x = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
        //double[] p = {0, 0.41, 0.79, 1.13, 1.46, 1.76, 2.04, 2.3, 2.55, 2.79, 3.01};
        double[] x = {0,2,3};
        double[] p = {8,4,9};
        double[] y = new double[p.length];
        for (int i = 0; i < p.length; i++) {
            y[i] = p[i] + Math.pow(-1, k) * m;
        }

        double point = 0.47;

        // Интерполяция в форме Лагранжа
        PolynomialFunction lagrangePoly = lagrangeInterpolation(x, y);
        System.out.println("Interpolation (Lagrange): " + lagrangePoly);

        // Интерполяция в форме Ньютона
        PolynomialFunction newtonPoly = newtonInterpolation(x, y);
        System.out.println("Interpolation (Newton): " + newtonPoly);

        // Вычисление значения многочленов в точке
        double lagrangeResult = lagrangePoly.value(point);
        double newtonResult = newtonInterpolation2(x, y, point);

        System.out.println("Interpolation result (Lagrange) at x = " + point + ": " + lagrangeResult);
        System.out.println("Interpolation result (Newton) at x = " + point + ": " + newtonResult);
    }

    private static PolynomialFunction lagrangeInterpolation(double[] x, double[] y) {
        int n = x.length;
        PolynomialFunction lagrangePoly = new PolynomialFunction(new double[]{0.0});

        for (int i = 0; i < n; i++) {
            PolynomialFunction term = new PolynomialFunction(new double[]{y[i]});
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    double denominator = x[i] - x[j];
                    term = term.multiply(new PolynomialFunction(new double[]{-x[j] / denominator, 1.0 / denominator}));
                }
            }
            lagrangePoly = lagrangePoly.add(term);
        }

        return lagrangePoly;
    }
    private static double newtonInterpolation2(double[] x, double[] y, double xi) {
        double result = y[0];
        double[] dividedDifference = computeDividedDifference(x, y);

        for (int i = 1; i < x.length; i++) {
            double term = dividedDifference[i];
            for (int j = 0; j < i; j++) {
                term *= (xi - x[j]);
            }
            result += term;
        }

        return result;
    }
    private static double[] computeDividedDifference(double[] x, double[] y) {
        int n = x.length;
        double[] dividedDifference = new double[n];
        System.arraycopy(y, 0, dividedDifference, 0, n);

        for (int i = 1; i < n; i++) {
            for (int j = n - 1; j >= i; j--) {
                dividedDifference[j] = (dividedDifference[j] - dividedDifference[j - 1]) / (x[j] - x[j - i]);
            }
        }

        return dividedDifference;
    }
    private static PolynomialFunction newtonInterpolation(double[] x, double[] y) {
        int n = x.length;
        PolynomialFunction newtonPoly = new PolynomialFunction(new double[]{y[0]});

        for (int i = 1; i < n; i++) {
            newtonPoly = newtonPoly.multiply(new PolynomialFunction(new double[]{-x[i - 1], 1}))
                    .add(new PolynomialFunction(new double[]{y[i]}));
        }

        return newtonPoly;
    }
}
