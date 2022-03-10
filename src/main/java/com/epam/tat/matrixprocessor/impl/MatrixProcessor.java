package com.epam.tat.matrixprocessor.impl;

import com.epam.tat.matrixprocessor.IMatrixProcessor;
import com.epam.tat.matrixprocessor.exception.MatrixProcessorException;
import java.util.logging.Logger;
import java.math.BigDecimal;
import java.math.RoundingMode;

/**
 * Function Description:
 * Complete the functions below. All methods must work with matrices of the double type.
 * <p>
 * Constraints:
 * 0 < m < 10
 * 0 < n < 10
 * where m - number of rows in matrix
 * where n - number of columns in matrix
 * <p>
 * In case of incorrect input values or inability to perform a calculation, the method should throw an appropriate
 * exception.
 */
public class MatrixProcessor implements IMatrixProcessor {
    private static final String ACT_1 = "Incorrect matrix value";
    private static final Logger logger = Logger.getGlobal();
    public static void main(String[] args) {
        int m = ((int) ((Math.random() + 0.1) * 9));
        int n = ((int) ((Math.random() + 0.1) * 9));
        double[][] matrix = new double[m][n];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] = ((Math.random() - Math.random()) * 100);
            }
        }
        logger.info("Start matrix");
        MatrixProcessor.printMatrix(matrix);
        MatrixProcessor processor = new MatrixProcessor();
        processor.transpose(matrix);
        processor.turnClockwise(matrix);
        double[][] firstMatrix = new double[m][n];
        for (int i = 0; i < firstMatrix.length; i++) {
            for (int j = 0; j < firstMatrix[i].length; j++) {
                firstMatrix[i][j] = ((Math.random()) - (Math.random()) * 100);
            }
        }
        double[][] secondMatrix = new double[m][n];
        for (int i = 0; i < secondMatrix.length; i++) {
            for (int j = 0; j < secondMatrix[i].length; j++) {
                secondMatrix[i][j] = ((Math.random() - Math.random()) * 100);
            }
        }
        processor.multiplyMatrices(firstMatrix, secondMatrix);
        processor.getInverseMatrix(matrix);
        processor.getMatrixDeterminant(matrix);
    }

    //public MatrixProcessor() {}

    public static void printMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.printf("%.5f", matrix[i][j]);
                System.out.print("\t");
            }
            System.out.println();
        }
    }

    /**
     * Matrix transpose is an operation on a matrix where its rows become columns with the same numbers.
     * Ex.:
     * |1 2|			|1 3 5|
     * |3 4|   ====>	|2 4 6|
     * |5 6|
     *
     * @param matrix - matrix for transposition
     * @return the transposed matrix
     */
    @Override // срабатывает только для симметричных матриц, но цикл проходит, а матрица не транспонируется
    public double[][] transpose(double[][] matrix) {
        if (matrix != null && matrix.length > 0 && matrix[0].length > 0) {
            double[][] temp = new double[matrix[0].length][matrix.length];
            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix[i].length; j++) {
                    temp[j][i] = matrix[i][j];
                }
            }
            matrix = temp;
            logger.info("Transposed matrix");
            printMatrix(matrix);
            return matrix;
        } else {
            throw new MatrixProcessorException(ACT_1);
        }
    }

    /**
     * The method flips the matrix clockwise.
     * Ex.:
     * * |1 2|			|5 3 1|
     * * |3 4|   ====>	|6 4 2|
     * * |5 6|
     *
     * @param matrix - rotation matrix
     * @return rotated matrix
     */
    @Override
    public double[][] turnClockwise(double[][] matrix) {
        if (matrix != null && matrix.length > 0 && matrix[0].length > 0) {
            double[][] temp = new double[matrix[0].length][matrix.length];
            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix[i].length; j++) {
                    temp[j][matrix.length - i - 1] = matrix[i][j];
                }
            }
            matrix = temp;
            logger.info("Turned clockwise matrix");
            printMatrix(matrix);
            return matrix;
        } else {
            throw new MatrixProcessorException(ACT_1);
        }
    }

    /**
     * This method multiplies matrix firstMatrix by matrix secondMatrix
     * <p>
     * See {https://en.wikipedia.org/wiki/Matrix_multiplication}
     *
     * @param firstMatrix  - first matrix to multiply
     * @param secondMatrix - second matrix to multiply
     * @return result matrix
     */
    @Override
    public double[][] multiplyMatrices(double[][] firstMatrix, double[][] secondMatrix) {
        if (firstMatrix != null && secondMatrix != null && firstMatrix.length > 0 && secondMatrix.length > 0 && firstMatrix[0].length > 0 && secondMatrix[0].length > 0 && firstMatrix[0].length == secondMatrix.length  ) {// возможно ошибка здесь и надо поменять местами стоки и столбцы
            double[][] matrix = new double[firstMatrix[0].length][secondMatrix.length];// здесь тоже поменял
            for (int i = 0; i < firstMatrix.length; i++) {
                for (int j = 0; j < secondMatrix[0].length; j++) {
                    for (int k = 0; k < secondMatrix.length; k++) {
                        matrix[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
                        BigDecimal value = BigDecimal.valueOf(matrix[i][j]);
                        value = value.setScale(3, RoundingMode.HALF_UP);
                        matrix[i][j] = value.doubleValue();

                    }
                }
            }
            logger.info("First");
            printMatrix(firstMatrix);
            logger.info("Second");
            printMatrix(secondMatrix);
            logger.info("Multiplied Matrix");
            printMatrix(matrix);
            return matrix;
        } else {
            throw new MatrixProcessorException(ACT_1);
        }
    }

    /**
     * This method returns the inverse of the matrix
     * <p>
     * See {https://en.wikipedia.org/wiki/Invertible_matrix}
     *
     * @param matrix - input matrix
     * @return inverse matrix for input matrix
     */
    @Override
    public double[][] getInverseMatrix(double[][] matrix) {
        if (matrix != null && matrix.length > 0 && matrix[0].length > 0 && matrix.length == matrix[0].length  ) {
            double temp;
            double[][] E = new double[matrix.length][matrix.length];
            for (int i = 0; i < matrix.length; i++)
                for (int j = 0; j < matrix.length; j++) {
                    E[i][j] = 0;
                    if (i == j)
                        E[i][j] = 1;
                }
            for (int k = 0; k < matrix.length; k++) {
                temp = matrix[k][k];
                for (int j = 0; j < matrix.length; j++) {
                    matrix[k][j] /= temp;
                    E[k][j] /= temp;
                }
                for (int i = k + 1; i < matrix.length; i++) {
                    temp = matrix[i][k];
                    for (int j = 0; j < matrix.length; j++) {
                        matrix[i][j] -= matrix[k][j] * temp;
                        E[i][j] -= E[k][j] * temp;
                    }
                }
            }
            for (int k = matrix.length - 1; k > 0; k--) {
                for (int i = k - 1; i >= 0; i--) {
                    temp = matrix[i][k];
                    for (int j = 0; j < matrix.length; j++) {
                        matrix[i][j] -= matrix[k][j] * temp;
                        E[i][j] -= E[k][j] * temp;
                    }
                }
            }
            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix.length; j++) {
                    matrix[i][j] = E[i][j];
                    BigDecimal value = BigDecimal.valueOf(matrix[i][j]);
                    value = value.setScale(3, RoundingMode.HALF_UP);
                    matrix[i][j] = value.doubleValue();
                }
            }
            logger.info("Inverted Matrix");
            printMatrix(matrix);
            return matrix;
        } else {
            throw new MatrixProcessorException(ACT_1);
        }
    }

    /**
     * This method returns the determinant of the matrix
     * <p>
     * See {https://en.wikipedia.org/wiki/Determinant}
     *
     * @param matrix - input matrix
     * @return determinant of input matrix
     */
    @Override
    public double getMatrixDeterminant(double[][] matrix) {
        if (matrix != null && matrix.length > 0 && matrix[0].length > 0 && matrix.length == matrix[0].length  ) {
            double calcResult = 0.0;
            if (matrix.length == 2) {
                calcResult = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
            } else {
                int koeff = 1;
                for (int i = 0; i < matrix.length; i++) {
                    if (i % 2 == 1) {  //я решил не возводить в степень, а просто поставить условие - это быстрее. Т.к. я раскладываю всегда по первой (читай - "нулевой") строке, то фактически я проверяю на четность значение i+0.
                        koeff = -1;
                    } else {
                        koeff = 1;
                    }
                    //собственно разложение:
                    calcResult += koeff * matrix[0][i] * this.getMatrixDeterminant(this.GetMinor(matrix, 0, i));
                }
            }
            //возвращаем ответ
            return calcResult;
        } else {
            throw new MatrixProcessorException(ACT_1);
        }
    }

    //функция, к-я возвращает нужный нам минор. На входе - определитель, из к-го надо достать минор и номера строк-столбцов, к-е надо вычеркнуть.
    private double[][] GetMinor(double[][] matrix, int row, int column) {
        int minorLength = matrix.length - 1;
        double[][] minor = new double[minorLength][minorLength];
        int dI = 0;//эти переменные для того, чтобы "пропускать" ненужные нам строку и столбец
        int dJ = 0;
        for (int i = 0; i <= minorLength; i++) {
            dJ = 0;
            for (int j = 0; j <= minorLength; j++) {
                if (i == row) {
                    dI = 1;
                } else {
                    if (j == column) {
                        dJ = 1;
                    } else {
                        minor[i - dI][j - dJ] = matrix[i][j];
                    }
                }
            }
        }
        return minor;
    }

}
