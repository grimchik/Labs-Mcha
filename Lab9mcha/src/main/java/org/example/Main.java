import java.util.function.Function;

public class Main {
    public static void main(String[] args) {
        double m = 2.0;
        double a = 0.9;
        double h = 0.1; // начальный шаг интегрирования
        double tolerance = 0.0001; // требуемая точность

        // Решение с методом Эйлера
        double resultEuler = solveEuler(m, a, h, tolerance);
        System.out.println("Метод Эйлера: y(1) = " + resultEuler);

        // Решение с методом Рунге-Кутты
        double resultRungeKutta = solveRungeKutta(m, a, h, tolerance);
        System.out.println("Метод Рунге-Кутты: y(1) = " + resultRungeKutta);
    }

    private static double solveEuler(double m, double a, double h, double tolerance) {
        double x = 0.0;
        double y = 0.0;

        while (x < 1.0 - h/2) {
            double yNext = y + h * equation(m, a, x, y);
            double yHalf = y + 0.5 * h * equation(m, a, x, y);
            double yNextHalf = yHalf + 0.5 * h * equation(m, a, x + 0.5 * h, yHalf);

            double error = Math.abs(yNextHalf - yNext) / 15.0;

            if (error < tolerance) {
                y = yNext;
                x += h;
            }

            // Адаптивное изменение шага
            h = Math.min(h * 2.0, Math.max(0.1, 0.9 * h * Math.sqrt(tolerance / error)));
        }

        return y;
    }

    private static double solveRungeKutta(double m, double a, double h, double tolerance) {
        double x = 0.0;
        double y = 0.0;

        while (x < 1.0 - h/2) {
            double k1 = h * equation(m, a, x, y);
            double k2 = h * equation(m, a, x + 0.5 * h, y + 0.5 * k1);
            double k3 = h * equation(m, a, x + 0.5 * h, y + 0.5 * k2);
            double k4 = h * equation(m, a, x + h, y + k3);

            double yNext = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
            double yHalf = y + 0.5 * h * equation(m, a, x, y);
            double yNextHalf = yHalf + 0.5 * h * equation(m, a, x + 0.5 * h, yHalf);

            double error = Math.abs(yNextHalf - yNext) / 15.0;

            if (error < tolerance) {
                y = yNext;
                x += h;
            }

            // Адаптивное изменение шага
            h = Math.min(h * 2.0, Math.max(0.1, 0.9 * h * Math.sqrt(tolerance / error)));
        }

        return y;
    }

    private static double equation(double m, double a, double x, double y) {
        return (a * (1 - y * y)) / ((1 + m) * x * x + y * y + 1);
    }
}
