package com.kazakov;

import org.la4j.LinearAlgebra;
import org.la4j.inversion.MatrixInverter;
import org.la4j.matrix.Matrix;
import org.la4j.matrix.dense.Basic2DMatrix;
import org.la4j.vector.Vector;
import org.la4j.vector.dense.BasicVector;

public class DoubleSimplex {
    public static Vector c = new BasicVector(new double[]{ -1, 2, -3, -1, 2, 3, 2 });
    public static Vector b = new BasicVector(new double[]{ 4, 5, 6 });
    public static Vector y = new BasicVector(new double[]{ 1, 1, 1 });
    public static Vector Jb = new BasicVector(new double[]{ 2, 4, 7 });
    public static Matrix A = new Basic2DMatrix(new double[][]{
            { 1, 2, 2, 0, -2, 8, 4 },
            { 1, 0, -1, 0, 4, 5, 1 },
            { 1, 0, 0, -1, 3, 6, 0 }
    });

    public static Matrix Ab = new Basic2DMatrix(new double[A.rows()][Jb.length()]); // Строки х Столбцы
    public static Matrix B;
    public static Vector PsevdoPlan;
    public static double E = 0.000001;

    public static void main(String[] args) {
        Init();

        Step1();
        if(Step2())
            System.out.println("план найден");
    }

    public static void Init() {
        for (int i = 0; i < Jb.length(); i++) {
            Ab.setColumn(i, A.getColumn((int) Jb.get(i)));
        }

        MatrixInverter inverter = Ab.withInverter(LinearAlgebra.GAUSS_JORDAN);
        B = inverter.inverse(LinearAlgebra.DENSE_FACTORY);
    }

    public static void Step1() {
        PsevdoPlan = B.multiply(b);
    }

    public static boolean Step2() {
        for (int i = 0; i < PsevdoPlan.length(); i++) {
            if(PsevdoPlan.get(i) < E)
                return false;
        }
        return true;
    }

    public static void Step3() {

    }

    public static void Step4() {

    }

    public static void Step5() {

    }
}
