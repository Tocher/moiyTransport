package com.kazakov;

import org.la4j.LinearAlgebra;
import org.la4j.inversion.MatrixInverter;
import org.la4j.matrix.Matrix;
import org.la4j.matrix.dense.Basic1DMatrix;
import org.la4j.matrix.dense.Basic2DMatrix;
import org.la4j.vector.Vector;
import org.la4j.vector.dense.BasicVector;

public class KvadrZadacha {
    public static Vector c = new BasicVector(new double[]{ -8, -6, -4, -6 });
    public static Vector b = new BasicVector(new double[]{ 2, 3 });
    public static Matrix A = new Basic2DMatrix(new double[][]{
            { 1, 0, 2, 1 },
            { 0, 1, -1, 2 }
    });
    public static Matrix D = new Basic2DMatrix(new double[][]{
            { 2, 1, 1, 0 },
            { 1, 1, 0, 0 },
            { 1, 0, 1, 0 },
            { 0, 0, 0, 0 }
    });

    public static Vector x = new BasicVector(new double[]{ 2, 3, 0, 0 });
    public static Vector Jon = new BasicVector(new double[]{ 1, 2 });
    public static Vector Jstar = new BasicVector(new double[]{ 1, 2 });
    public static Vector Jn = new BasicVector(new double[x.length() - Jstar.length()]);

    public static Vector Delta = new BasicVector(new double[x.length()]);
    public static Vector Teta = new BasicVector(new double[x.length()]);
    public static Vector Lstar = new BasicVector(new double[Jstar.length()]);
    public static double DeltaSymbol;
    public static double Teta0;
    public static double E = 0.000001;
    public static double J0;
    public static double Jz;
    public static boolean toStep4 = false;
    public static boolean debbug = false;

    public static void main(String[] args) {
        JnInit();

        while(true) {
            if(!toStep4) {
                Step1();
                if (Step2()) {
                    System.out.println("решение найдено");
                    break;
                }
                Step3();
                toStep4 = false;
            }
            Step4();
            if (!Step5()) {
                System.out.println("задача не имеет решения в силу неограниченности снизу целевой функции на множестве планов");
                break;
            }
            Step6();
            switch(Step7()) {
                case 4:
                    toStep4 = true;
                    break;
                case 0:
                    debbug = true;
                    break;
            }
            if(debbug)
                break;
        }

    }

    public static void Step1() {
        Vector cx = c.add(D.multiply(x));
        // Вектор потенциалов
        Vector cOn = new BasicVector(new double[Jon.length()]);
        Matrix AOn = new Basic2DMatrix(new double[Jon.length()][A.rows()]); // первым колонки или строки?
        for (int i = 0; i < Jon.length(); i++) {
            cOn.set(i, cx.get((int) Jon.get(i) - 1)); // отсчет с 0
            AOn.setColumn(i, A.getColumn((int) Jon.get(i) - 1));
        }
        MatrixInverter inverter = AOn.withInverter(LinearAlgebra.GAUSS_JORDAN);
        AOn = inverter.inverse(LinearAlgebra.DENSE_FACTORY);

        Vector ux = cOn.multiply(-1).multiply(AOn);

        for (int i = 0; i < Jn.length(); i++) {
            int k = (int) Jn.get(i) - 1;

            Vector tmp = A.getColumn(k);
            tmp = tmp.multiply(ux.toColumnMatrix());
            Delta.set(k, tmp.get(0) + cx.get(k));
        }
    }

    public static boolean Step2() {
        for (int i = 0; i < Delta.length(); i++) {
            if(Delta.get(i) < E)
                return false;
        }
        return true;
    }

    public static void Step3() {
        for (int i = 0; i < Jn.length(); i++) {
            if(Delta.get((int) Jn.get(i) - 1) < -E)
            {
                J0 = (int) Jn.get(i); // первое нахождение
                break;
            }
        }
        System.out.println("Jn" + Jn);
        System.out.println("Delta" + Delta);
        System.out.println("J0" + J0);
    }

    public static void Step4() {
        Vector DstarJ0 = new BasicVector(getDJ0());
        Vector AJ0 = A.getColumn((int) J0 - 1);
        Matrix Astar = new Basic2DMatrix(new double[A.rows()][Jstar.length()]);
        System.out.println("Jstar" + Jstar);
        for (int i = 0; i < Jstar.length(); i++) {
            Astar.setColumn(i, A.getColumn((int) Jstar.get(i) - 1));
        }
        Matrix Dstar = new Basic2DMatrix(new double[Jstar.length()][Jstar.length()]);
        for (int i = 0; i < Jstar.length(); i++) {
            for (int j = 0; j < Jstar.length(); j++) {
                Dstar.set(i,j, D.get(i, j));
            }
        }
        System.out.println("Dstar" + Dstar);
        System.out.println("Astar" + Astar);
        System.out.println("DstarJ0" + DstarJ0);
        System.out.println("AJ0" + AJ0);

        MatrixInverter inverter = Astar.withInverter(LinearAlgebra.GAUSS_JORDAN);
        Matrix Bstar = inverter.inverse(LinearAlgebra.DENSE_FACTORY);

        Lstar = AJ0.multiply(-1).multiply(Bstar);

        Vector tmp = Dstar.multiply(Lstar);
        tmp = tmp.add(DstarJ0).multiply(-1);

        Bstar = Bstar.transpose();
        Vector y = tmp.multiply(Bstar);

        int j0 = (int) J0 - 1;
        DeltaSymbol = DstarJ0.multiply(Lstar.toColumnMatrix()).get(0) + AJ0.multiply(y.toColumnMatrix()).get(0) + D.get(j0, j0);
    }

    public static boolean Step5() {
        int j0 = (int) J0 - 1;
        if( DeltaSymbol < E && DeltaSymbol > -E)
            Teta.set(j0, Double.POSITIVE_INFINITY);
        else {
            System.out.println();
            Teta.set(j0, Math.abs(Delta.get(j0) / DeltaSymbol));
        }
        Teta0 = Teta.get(j0);
        Jz = j0 + 1;

        for (int i = 0; i < Jstar.length(); i++) {
            if (Lstar.get(i) < -E) {
                Teta.set(i, (-1 * x.get((int) Jstar.get(i) - 1) / Lstar.get(i))); // МОГУТ БЫТЬ КОСЯКИ
            } else {
                Teta.set(i, Double.POSITIVE_INFINITY);
            }

            if(Teta0 > Teta.get(i)) {
                Teta0 = Teta.get(i);
                Jz = i + 1;
            }
        }

        if(Teta0 == Double.POSITIVE_INFINITY)
            return false;
        return true;
    }

    public static void Step6() {
        for (int i = 0; i < Jstar.length(); i++) {
            int j = (int) Jstar.get(i) - 1;
            x.set(j, (x.get(j) + Lstar.get(j) * Teta0));
        }
        for (int i = 0; i < Jn.length(); i++) {
            int j = (int) Jn.get(i) - 1;
            x.set(j, 0);
        }
        x.set((int) J0 - 1, (x.get((int) J0 - 1) + Teta0));
    }

    public static int Step7() {
        if(Jz == J0) { // Случай а
            System.out.println("Случай а");
            Vector tmp = new BasicVector(new double[Jstar.length()+1]);
            boolean J0inserted = false;
            int k = 0;
            for (int i = 0; i < Jstar.length(); i++) {
                if(Jstar.get(i) > J0 && !J0inserted) {
                    tmp.set(k, J0);
                    J0inserted = true;
                    k++;
                }
                tmp.set(k, Jstar.get(i));
                k++;
            }
            if(!J0inserted)
                tmp.set(k, J0);
            Jstar = tmp;
            tmp = new BasicVector(new double[Jn.length()-1]);
            k = 0;
            for (int i = 0; i < Jn.length(); i++) {
                if(Jn.get(i) != J0) {
                    tmp.set(k, Jn.get(i));
                    k++;
                }
            }
            return 1;
        }
        else if (isCaseB()) {
            System.out.println("Случай б");
            Delta.set((int) J0 - 1 ,Delta.get((int) J0 - 1) + Teta0 * DeltaSymbol);
            Vector tmp = new BasicVector(new double[Jstar.length()-1]);
            int k = 0;
            for (int i = 0; i < Jstar.length(); i++) {
                if(Jstar.get(i) != J0) {
                    tmp.set(k, Jstar.get(i));
                    k++;
                }
            }
            Jstar = tmp;
            return 4;
        }
        else if(isCaseC()) {
            System.out.println("Случай в");
            return 0;
        }
        else {
            System.out.println("Случай г");
            return 0;
        }
    }

    public static boolean isJSinJON(double JS) {
        for (int j = 0; j < Jon.length(); j++) {
            if(JS == Jon.get(j))
                return true;
        }
        return false;
    }

    public static boolean isCaseB() {
        for (int i = 0; i < Jstar.length(); i++) {
            if(isJSinJON(Jstar.get(i)))
                continue;
            if(Jstar.get(i) == Jz)
                return true;
        }
        return false;
    }

    public static boolean isCaseC() {
        Matrix Aon = new Basic1DMatrix(new double[Jon.length()][Jon.length()]);
        for (int i = 0; i < Jon.length(); i++) {
            Aon.setColumn(i, A.getColumn((int) Jon.get(i) - 1));
        }
        MatrixInverter inverter = Aon.withInverter(LinearAlgebra.GAUSS_JORDAN);
        Aon = inverter.inverse(LinearAlgebra.DENSE_FACTORY);

        for (int i = 0; i < Jon.length(); i++) {
            if (Jon.get(i) == Jz) {
                for (int j = 0; j < Jstar.length(); j++) {
                    if (isJSinJON(Jstar.get(j)))
                        continue;
                    Vector Ajplus = A.getColumn((int) Jstar.get(j) - 1); // Aj+
                    Vector es = getIdentityColumn(i, Jon.length());
                    Vector tmp = es.multiply(Aon).toColumnMatrix().multiply(Ajplus);
                    if (tmp.get(0) != 0)
                        return true;
                }
            }
        }
        return false;
    }

    public static double[] getDJ0() {
        double[] arr = new double[Jstar.length()];
        for(int i = 0; i < Jstar.length(); i++)
            arr[i] = D.get(i, (int) J0 - 1);
        return arr;
    }

    public static void JnInit() {
        int Ji = 0;
        int Jni = 0;
        for (int i = 0; i < x.length(); i++) {
            if( Jstar.length() > Ji) {
                if (Jstar.get(Ji) == i + 1)
                    Ji++;
            }
            else {
                Jn.set(Jni, i + 1);
                Jni++;
            }
        }
    }

    public static Vector getIdentityColumn(int col, int size) {
        Matrix e = new Basic1DMatrix(new double[size][size]);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if(i == j)
                    e.set(i, j, 1);
            }
        }
        return e.getColumn(col);
    }
}
