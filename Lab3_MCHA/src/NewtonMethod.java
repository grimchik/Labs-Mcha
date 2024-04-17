import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

public class NewtonMethod {

    // Максимальное количество итераций
    private static final int MAX_ITERATIONS = 1000;
    // Точность
    private static final double EPSILON = 1e-6;

    // Функция f(x)
    public static double f(double x) {
        //return x*x-0.9*x-0.36;
        //return x * x * x - 19.7997 * x * x + 28.9378 * x + 562.833;
        return x*x - 4*x;
        //return x;
    }
    // Производная функции f'(x)
    public static double df(double x) {
        //return 2 * x  - 0.9 ;
        return 2*x-4;
        //return 1;
    }

    public static double newtonMethod(double x0) {
        double x = x0;
        int iterations = 0;

        while (iterations < MAX_ITERATIONS) {
            double fValue = f(x);
            double dfValue = df(x);

            if (Math.abs(dfValue) < EPSILON) {
                throw new RuntimeException("Derivative is too small, solution might not converge.");
            }

            x = x - fValue / dfValue;

            if (Math.abs(fValue) < EPSILON) {
                break;
            }

            iterations++;
        }

        return x;
    }

    public static List<Function<Double, Double>> generateSturmSequence() {
        List<Function<Double, Double>> sequence = new ArrayList<>();
        sequence.add(NewtonMethod::f);
        sequence.add(NewtonMethod::df);

        for (int i = 1; i < sequence.size() && sequence.get(i - 1).apply(1.0) != 0; i++) {
            Function<Double, Double> fPrev = sequence.get(i - 1);
            Function<Double, Double> fCur = sequence.get(i);
            sequence.add(x -> -fPrev.apply(x) % fCur.apply(x));
        }

        return sequence;
    }

    public static int signChanges(double x, List<Function<Double, Double>> sequence) {
        int changes = 0;
        double prevSign = Math.signum(sequence.get(0).apply(x));

        for (int i = 1; i < sequence.size(); i++) {
            double currentSign = Math.signum(sequence.get(i).apply(x));
            if (prevSign != currentSign) changes++;
            prevSign = currentSign;
        }

        return changes;
    }
    public static double bisectMethod(double a, double b) {
       //if (f(a) * f(b) > 0) {
       //     throw new IllegalArgumentException("f(a) and f(b) must have opposite signs.");
       //}

        double mid = a;
        while ((b - a) / 2.0 > EPSILON) {
            mid = (a + b) / 2.0;

            if (Math.abs(f(mid)) < EPSILON) break;

            if (f(a) * f(mid) < 0) {
                b = mid;
            } else {
                a = mid;
            }
        }

        return mid;
    }
    public static double chordMethod(double x0, double x1) {
        while (Math.abs(x1 - x0) > EPSILON) {
            double temp = x1;
            x1 = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));
            x0 = temp;
        }
        return x1;
    }
    public static int rootsCount(double a, double b) {
        List<Function<Double, Double>> sequence = generateSturmSequence();
        return signChanges(a, sequence) - signChanges(b, sequence);
    }

    public static void main(String[] args) {
        double a1 = -2;
        double b1 = 6;
        int count = rootsCount(a1, b1);
        System.out.println("Number of roots in the interval [" + a1 + ", " + b1 + "]: " + count);
        double x0 = -10;
        double root = newtonMethod(x0);
        System.out.println("Root found by Newton's method: " + root);

        double a = -10;
        double b = 1.5;
        double rootBisect = bisectMethod(a, b);
        System.out.println("Root found by Bisect method: " + rootBisect);
        double x1 = 0.0001;
        double rootChord = chordMethod(x0, x1);
        System.out.println("Root found by Chord method: " + rootChord);
    }
}
