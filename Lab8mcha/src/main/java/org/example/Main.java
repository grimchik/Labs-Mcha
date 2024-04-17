import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math3.random.RandomDataGenerator;

import java.util.function.DoubleUnaryOperator;

public class Main {

    private static final double epsilon = 0.000001;

    public static double firstDerivative(UnivariateFunction f, double x, UnivariateFunction f2, UnivariateFunction f3) {
        double M2 = Math.abs(f2.value(x));
        double M3;
        double h;
        if (M2 == 0) {
            M3 = Math.abs(f3.value(x));
            h = Math.pow((6 * epsilon / M3), 0.5);
        } else {
            double h2 = 2 * epsilon / M2;
            M3 = Math.abs(f3.value(x));
            double h3 = Math.pow((6 * epsilon / M3), 0.5);
            h = Math.min(h2, h3);
        }
        return (f.value(x + h) - f.value(x - h)) / (2 * h);
    }

    public static double secondDerivative(UnivariateFunction f, double x, UnivariateFunction f4) {
        double M4 = Math.abs(f4.value(x));
        double h2 = Math.abs((12 * epsilon / M4));
        double h = Math.abs(Math.pow(h2, 0.5));
        return (f.value(x + h) - 2 * f.value(x) + f.value(x - h)) / h2;
    }

    public static double secondDerivative2(UnivariateFunction f, double x, UnivariateFunction f2, UnivariateFunction f4) {
        double M4 = Math.abs(f4.value(x));
        double h2 = Math.abs((12 * epsilon / M4));
        double h = Math.abs(Math.pow(h2, 0.5));
        double a = (f.value(x + h) - 2 * f.value(x) + f.value(x - h)) / h2;
        while (a - f2.value(x) > epsilon) {
            h -= Math.pow(epsilon, 3);
            a = (f.value(x + h) - 2 * f.value(x) + f.value(x - h)) / h2;
        }
        return a;
    }

    private static final double EPSILON = 0.000001;

    public static double intMiddleRectangle(UnivariateFunction f, double l, double r, int n) {
        double h = (r - l) / n;
        double x = l + h / 2;
        double s = 0.0;
        while (x < r) {
            s += f.value(x) * h;
            x += h;
        }
        return s;
    }

    public static double integralTrapezoid(UnivariateFunction f, double l, double r, int n) {
        double h = (r - l) / n;
        double x = l + h / 2;
        double s = 0.0;
        while (x < r) {
            s += ((f.value(x - h / 2) + f.value(x + h / 2)) / 2) * h;
            x += h;
        }
        return s;
    }

    public static double integralSimpson(UnivariateFunction f, double l, double r, int n) {
        SimpsonIntegrator integrator = new SimpsonIntegrator();
        return integrator.integrate(n, f, l, r);
    }

    public static double integralRandSegments(UnivariateFunction f, double l, double r) {
        double leftCoeff = 1 / 3.0;
        double rightCoeff = 1 / 2.0;
        double hPrev = r - l;
        double ansPrev = integralTrapezoid(f, l, r, 1);
        RandomDataGenerator randomDataGenerator = new RandomDataGenerator();
        while (true) {
            double hNew = hPrev * (leftCoeff + (rightCoeff - leftCoeff) * randomDataGenerator.nextUniform(0, 1));
            int n = (int) Math.floor((r - l) / hNew);
            double m = l + hNew * n;
            double ansNew = integralTrapezoid(f, l, m, n) + integralTrapezoid(f, m, r, 1);
            if (Math.abs(ansNew - ansPrev) < EPSILON) {
                System.out.println("\nN : " + n);
                return ansNew;
            }
            ansPrev = ansNew;
            hPrev = hNew;
        }
    }

    public static double integralMiddleRectangleViaEstimation(UnivariateFunction f, double l, double r, double epsilon, double m2deLR) {
        if (m2deLR > 0.0) {
            double M2 = m2deLR;
            double h = Math.sqrt(24 * epsilon / (r - l) / M2);
            int n = (int) Math.ceil((r - l) / h);
            return intMiddleRectangle(f, l, r, n);
        } else {
            return integralRandSegments(f, l, r);
        }
    }

    public static double integralTrapezoidViaEstimation(UnivariateFunction f, double l, double r, double epsilon, double m2deLR) {
        if (m2deLR > 0.0) {
            double M2 = m2deLR;
            double h = Math.sqrt(12 * epsilon / (r - l) / M2);
            int n = (int) Math.ceil((r - l) / h);
            return integralTrapezoid(f, l, r, n);
        } else {
            return integralRandSegments(f, l, r);
        }
    }

    public static double integralSimpsonViaEstimation(DoubleUnaryOperator f, double l, double r, double epsilon, double m4deLR) {
        if (m4deLR > 0.0) {
            double M4 = m4deLR;
            double h = Math.pow((180 * epsilon / (r - l) / M4), 0.25);
            int n = (int) Math.ceil((r - l) / h);
            return integralSimpson(f, l, r, n);
        } else {
            int n = 1; // Добавлено объявление переменной n
            return integralRandSegments(x -> integralSimpson(f, l, x, n), l, r);
        }
    }


    private static double integralSimpson(DoubleUnaryOperator f, double l, double r, int n) {
        double h = (r - l) / n;
        double x = l + h / 2;
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            s += f.applyAsDouble(x) * h;
            x += h;
        }
        return s;
    }

    private static double integralRandSegments(DoubleUnaryOperator integrator, double l, double r) {
        double leftCoeff = 1 / 3.0;
        double rightCoeff = 1 / 2.0;
        double hPrev = r - l;
        double ansPrev = integrator.applyAsDouble(l);
        RandomDataGenerator randomDataGenerator = new RandomDataGenerator();
        while (true) {
            double hNew = hPrev * (leftCoeff + (rightCoeff - leftCoeff) * randomDataGenerator.nextUniform(0, 1));
            int n = (int) Math.floor((r - l) / hNew);
            double m = l + hNew * n;
            double ansNew = integrator.applyAsDouble(l) + integrator.applyAsDouble(m);
            if (Math.abs(ansNew - ansPrev) < EPSILON) {
                System.out.println("\nN : " + n);
                return ansNew;
            }
            ansPrev = ansNew;
            hPrev = hNew;
        }
    }
}
