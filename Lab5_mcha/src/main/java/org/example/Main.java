import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import java.text.DecimalFormat;

public class Main {
    public static void main(String[] args) {
        // Ваши матрицы C и B
        double[][] matrixCres = {
                {0.2, 0, 0.2, 0, 0},
                {0, 0.2, 0, 0.2, 0},
                {0.2, 0, 0.2, 0, 0.2},
                {0, 0.2, 0, 0.2, 0},
                {0, 0, 0.2, 0, 0.2}
        };
        double[][] matrixBres = {
                {2.33, 0.81, 0.67, 0.92, -0.53},
                {0.81, 2.33, 0.81, 0.67, 0.92},
                {0.67, 0.81, 2.33, 0.81, 0.92},
                {0.92, 0.67, 0.81, 2.33, -0.53},
                {-0.53, 0.92, 0.92, -0.53, 2.33}
        };

        int scalar = 3;

        // Умножаем матрицу C на число
        int rowsC = matrixCres.length;
        int columnsC = matrixCres[0].length;

        double[][] result = new double[rowsC][columnsC];

        for (int i = 0; i < rowsC; i++) {
            for (int j = 0; j < columnsC; j++) {
                result[i][j] = matrixCres[i][j] * scalar;
            }
        }

        // Складываем результат с матрицей B
        int rowsB = matrixBres.length;
        int columnsB = matrixBres[0].length;

        for (int i = 0; i < rowsB; i++) {
            for (int j = 0; j < columnsB; j++) {
                result[i][j] = result[i][j] + matrixBres[i][j];
            }
        }
        double[][] resultres2 =
                {
                        {5,1,2},
                        {1,4,1},
                        {2,1,3}
                };
        RealMatrix matrixA = new Array2DRowRealMatrix(resultres2);

        // Вычисление собственных значений и векторов
        EigenDecomposition eigenDecomposition = new EigenDecomposition(matrixA);

        // Получение собственных значений
        double[] eigenvalues = eigenDecomposition.getRealEigenvalues();

        // Получение собственных векторов
        RealVector[] eigenvectors = new RealVector[matrixA.getColumnDimension()];
        for (int i = 0; i < eigenvectors.length; i++) {
            eigenvectors[i] = eigenDecomposition.getEigenvector(i);
        }

        // Вывод результатов
        DecimalFormat format = new DecimalFormat("0.0000");

        System.out.println("\nСобственные значения (Apache Commons Math):");
        for (double value : eigenvalues) {
            System.out.println(format.format(value));
        }

        System.out.println("\nСобственные векторы (Apache Commons Math):");
        for (RealVector vector : eigenvectors) {
            System.out.println(formatVector(vector.toArray()));
        }
        // Ваши матрицы C и B из первой части
        double[][] matrixC = {
                {0.2, 0, 0.2, 0, 0},
                {0, 0.2, 0, 0.2, 0},
                {0.2, 0, 0.2, 0, 0.2},
                {0, 0.2, 0, 0.2, 0},
                {0, 0, 0.2, 0, 0.2}
        };
        double[][] matrixB = {
                {2.33, 0.81, 0.67, 0.92, -0.53},
                {0.81, 2.33, 0.81, 0.67, 0.92},
                {0.67, 0.81, 2.33, 0.81, 0.92},
                {0.92, 0.67, 0.81, 2.33, -0.53},
                {-0.53, 0.92, 0.92, -0.53, 2.33}
        };

        int scalarres = 3;

        // Умножаем матрицу C на число
        int rows1res = matrixC.length;
        int columns1res = matrixC[0].length;

        double[][] resultres = new double[rows1res][columns1res];

        for (int i = 0; i < rows1res; i++) {
            for (int j = 0; j < columns1res; j++) {
                resultres[i][j] = matrixC[i][j] * scalarres;
            }
        }

        // Складываем результат с матрицей B
        int rowsres = matrixB.length;
        int columnsres = matrixB[0].length;

        for (int i = 0; i < rowsres; i++) {
            for (int j = 0; j < columnsres; j++) {
                resultres[i][j] = resultres[i][j] + matrixB[i][j];
            }
        }

        double epsilon = 0.0001; // Критерий остановки (может потребоваться настройка)

        JacobiResult matrixres = jacobi(resultres2, epsilon);

        System.out.println("\nСобственные значения :");
        for (double eigenvalue : matrixres.eigenvalues) {
            System.out.println(format.format(eigenvalue));
        }

        System.out.println("\nСобственные векторы :");
        for (double[] eigenvector : matrixres.eigenvectors) {
            System.out.println(formatVector(eigenvector));
        }
        System.out.println("\nКоличество итераций:");
        System.out.println(matrixres.iterations);

    }
    private static String formatVector(double[] vector) {
        DecimalFormat format = new DecimalFormat("0.0000");
        StringBuilder result = new StringBuilder();
        for (double value : vector) {
            result.append(format.format(value)).append(" ");
        }
        return result.toString();
    }
    static class JacobiResult {
        int iterations;
        double[] eigenvalues;
        double[][] eigenvectors;

        JacobiResult(double[] eigenvalues, double[][] eigenvectors, int iterations) {
            this.iterations = iterations;
            this.eigenvalues = eigenvalues;
            this.eigenvectors = eigenvectors;
        }
    }

    static JacobiResult jacobi(double[][] matrix, double epsilon) {
        int n = matrix.length;

        double[][] eigenVectors = new double[n][n];
        for (int i = 0; i < n; i++) {
            eigenVectors[i][i] = 1.0;
        }

        int maxIterations = 1000; // Максимальное количество итераций (может потребоваться настройка)
        int iterations = 0;

        while (true) {
            int p = 0, q = 0;
            double maxOffDiagonal = 0.0;

            // Находим максимальное по модулю значение внедиагональных элементов
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    double offDiagonal = Math.abs(matrix[i][j]);
                    if (offDiagonal > maxOffDiagonal) {
                        maxOffDiagonal = offDiagonal;
                        p = i;
                        q = j;
                    }
                }
            }

            // Проверяем критерий остановки
            if (maxOffDiagonal < epsilon || iterations >= maxIterations) {
                break;
            }

            // Вычисляем угол поворота (theta)
            double theta = 0.5 * Math.atan2(2 * matrix[p][q], matrix[p][p] - matrix[q][q]);

            // Создаем матрицу вращения
            double[][] rotationMatrix = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i == j) {
                        rotationMatrix[i][j] = 1.0;
                    } else {
                        rotationMatrix[i][j] = 0.0;
                    }
                }
            }
            rotationMatrix[p][p] = Math.cos(theta);
            rotationMatrix[q][q] = Math.cos(theta);
            rotationMatrix[p][q] = -Math.sin(theta);
            rotationMatrix[q][p] = Math.sin(theta);

            // Выполняем вращение матрицы
            matrix = multiply(transpose(rotationMatrix), multiply(matrix, rotationMatrix));
            eigenVectors = multiply(eigenVectors, rotationMatrix);

            iterations++;
        }

        double[] eigenvalues = new double[n];
        for (int i = 0; i < n; i++) {
            eigenvalues[i] = matrix[i][i];
        }

        return new JacobiResult(eigenvalues, eigenVectors, iterations);
    }

    static double[][] multiply(double[][] a, double[][] b) {
        int rowsA = a.length;
        int colsA = a[0].length;
        int colsB = b[0].length;

        double[][] result = new double[rowsA][colsB];

        for (int i = 0; i < rowsA; i++) {
            for (int j = 0; j < colsB; j++) {
                for (int k = 0; k < colsA; k++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        return result;
    }

    static double[][] transpose(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;

        double[][] result = new double[cols][rows];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[j][i] = matrix[i][j];
            }
        }

        return result;
    }
}
