package com.kazakov;

import org.la4j.matrix.Matrix;
import org.la4j.matrix.dense.Basic2DMatrix;
import org.la4j.vector.Vector;
import org.la4j.vector.dense.BasicVector;

import java.util.ArrayList;

public class Main {
    public static Vector A = new BasicVector(new double[]{ 100, 120, 150, 50 }); // Запасы
    public static Vector B = new BasicVector(new double[]{ 70, 50, 60, 30, 50, 70, 90 }); // Потребности
    public static Matrix C = new Basic2DMatrix(new double[][]{
            { 15, 23, 28, 13, 12, 16, 21 },
            { 24, 17, 24, 11, 11, 12, 18 },
            { 8, 19, 16, 18, 9, 17, 15 },
            { 12, 28, 18, 16, 21, 20, 25 }
    });
/*
    public static Vector A = new BasicVector(new double[]{ 40, 35, 10, 25 }); // Запасы
    public static Vector B = new BasicVector(new double[]{ 5, 15, 10, 10, 20, 30, 20 }); // Потребности
    public static Matrix C = new Basic2DMatrix(new double[][]{
            { 10, 15, 3, 10, 8, 1, -15 },
            { 7, 8, 4, 1, 6, -6, 7 },
            { 4, -3, 6, 9, 12, 6, 13 },
            { -10, 5, 8, -4, 2, 5, 11 }
    });
*/
/*
    public static Vector A = new BasicVector(new double[]{ 130, 55, 80, 65, 135 }); // Запасы
    public static Vector B = new BasicVector(new double[]{ 130, 75, 65, 60, 75, 60 }); // Потребности
    public static Matrix C = new Basic2DMatrix(new double[][]{
        { 6, 6, 8, 5, 4, 3 },
        { 2, 4, 3, 9, 8, 5 },
        { 3, 5, 7, 9, 6, 11 },
        { 3, 5, 4, 4, 2, 1 },
        { 2, 5, 6, 3, 2, 8 },
    });
*/

/*
    public static Vector A = new BasicVector(new double[]{ 50, 40, 60, 31 }); // Запасы | Предложение
    public static Vector B = new BasicVector(new double[]{ 30, 50, 20, 40, 30, 11 }); // Потребности | Спрос
    public static Matrix C = new Basic2DMatrix(new double[][]{
            { 2, 1, 3, 3, 2, 5 },
            { 3, 2, 2, 4, 3, 4 },
            { 3, 5, 4, 2, 4, 1 },
            { 4, 2, 2, 1, 2, 2 }
    });
*/
    public static int M = A.length();
    public static int N = B.length();
    public static Matrix D;
    public static int numCount;
    public static int maxI;
    public static int maxJ;
    public static int max_iterations;
    public static Matrix BasisCells = new Basic2DMatrix(new double[A.length()][B.length()]);
    public static ArrayList<TranspNode> path = new ArrayList<TranspNode>();
    public static Vector U = new BasicVector(new double[A.length()]);
    public static Vector V = new BasicVector(new double[B.length()]);

    public static void metodSeveroZapad() {
        D = new Basic2DMatrix(new double[A.length()][B.length()]);
        Vector tempA = A;
        Vector tempB = B;
        double min = 0;
        numCount = 0;
        for(int i = 0; i < A.length(); i++) {
            for(int j = 0; j < B.length(); j++) {
                if(tempA.get(i) == 0) {
                    break;
                }
                if(tempB.get(j) == 0)
                    continue;

                if(tempA.get(i) <= tempB.get(j)) {
                    min = tempA.get(i);
                    tempA.set(i, 0);
                    tempB.set(j, tempB.get(j)-min);
                }
                else {
                    min = tempB.get(j);
                    tempB.set(j, 0);
                    tempA.set(i, tempA.get(i)-min);
                }
                D.set(i, j, min);
                //BasisCells.set(i, j, 1);
                numCount++;
            }
        }
    }

    public static boolean getPotencial() {
        for (int i = 0; i < U.length(); i++) {
            U.set(i, Double.NEGATIVE_INFINITY);
        }
        for (int i = 0; i < V.length(); i++) {
            V.set(i, Double.NEGATIVE_INFINITY);
        }
        U.set(0, 0);

        max_iterations = A.length() * B.length();
        getPotencialHorizontal(0);

        // Проверка вычислений
        for (int i = 0; i < U.length(); i++) {
            if( U.get(i) == Double.NEGATIVE_INFINITY ){
                System.out.println("Не удалось вычислить потенциал u["+i+"]");
                return false;
            }
        }
        for (int i = 0; i < V.length(); i++) {
            if( V.get(i) == Double.NEGATIVE_INFINITY ){
                System.out.println("Не удалось вычислить потенциал v["+i+"]");
                return false;
            }
        }
        return true;
    }

    public static void getPotencialVertical(int j) {
        if( V.get(j) == Double.NEGATIVE_INFINITY) {
            System.out.println("Ошибка получения потенциала v["+j+"]");
        }
        for (int i = 0; i < A.length(); i++) {
            if(BasisCells.get(i,j) == 0)
                continue;
            if(U.get(i) != Double.NEGATIVE_INFINITY)
                continue;
            else {
                U.set(i, C.get(i,j) - V.get(j));
                getPotencialHorizontal(i);
            }
        }
    }

    public static void getPotencialHorizontal(int i) {
        max_iterations--;
        if( max_iterations == 0) {
            System.out.println("Зацикливание при вычислении потенциалов");
            return;
        }
        if( U.get(i) == Double.NEGATIVE_INFINITY) {
            System.out.println("Ошибка получения потенциала u["+i+"]");
        }
        for (int j = 0; j < B.length(); j++) {
            if(BasisCells.get(i,j) == 0)
                continue;
            if(V.get(j) != Double.NEGATIVE_INFINITY)
                continue;
            else {
                V.set(j, C.get(i,j) - U.get(i));
                getPotencialVertical(j);
            }
        }
    }

    public static boolean getCycle(int i0, int j0) {
        max_iterations = A.length() * B.length();
        path.clear();

        if( getHorizontalCycle(i0, j0) )
            return true;
        return false;
    }

    public static boolean getHorizontalCycle(int i0, int j0) {
        max_iterations--;
        if(max_iterations == 0) {
            System.out.println("Слишком много итераций поиска цикла.");
            return false;
        }
        for (int j = 0; j < B.length(); j++) {
            if(j == j0)
                continue;
            if(BasisCells.get(i0,j) == 0) {
                continue;
            }
            if(getVerticalCycle(i0, j)) {
                path.add(new TranspNode(i0, j));
                return true;
            }
        }
        return false;
    }

    public static boolean getVerticalCycle(int i0, int j0) {
        for (int i = 0; i < A.length(); i++) {
            if(j0 == maxJ && i == maxI) {
                path.add(new TranspNode(i, j0));
                return true;
            }
            if(i == i0) {
                continue;
            }
            if(BasisCells.get(i,j0) == 0) {
                continue;
            }
            if(getHorizontalCycle(i, j0)) {
                path.add(new TranspNode(i, j0));
                return true;
            }
        }
        return false;
    }

    public static double getPrice() {
        double Sum = 0;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                Sum += (D.get(i,j)*C.get(i,j));
            }
        }
        return Sum;
    }

    public static boolean refreshPlan() {
        double vol;
        double teta = Double.MIN_VALUE;
        boolean sign = true;
        for (int q = 0; q < path.size(); q++) {
            if(!sign) {
                vol = D.get(path.get(q).i, path.get(q).j);
                if( teta == Double.MIN_VALUE)
                    teta = vol;
                else {
                    if( vol < teta )
                        teta = vol;
                }
            }
            sign = !sign;
        }
        if( teta == Double.MIN_VALUE)
            System.out.println("ошибка вычислений тета");
        if( teta == 0)
            return false;
        sign = true;
        for (int q = 0; q < path.size(); q++) {
            int i = path.get(q).i;
            int j = path.get(q).j;
            if(!sign) {
                D.set(i, j, D.get(i, j) - teta);
            }
            else {
                D.set(i, j, D.get(i, j) + teta);
            }
            sign = !sign;
        }
        return true;
    }

    public static boolean isOptimum() {
        boolean result = true;
        double delta;
        double minDelta = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if(BasisCells.get(i,j) == 1) {
                    delta = 0;
                }
                else
                    delta = C.get(i,j) - U.get(i) - V.get(j);

                if(delta < 0)
                    result = false;
                if(minDelta == Double.NEGATIVE_INFINITY) {
                    minDelta = delta;
                    maxI = i;
                    maxJ = j;
                }
                else {
                    if( delta < minDelta) {
                        minDelta = delta;
                        maxI = i;
                        maxJ = j;
                    }
                }
            }
        }
        return result;
    }

    public static void printCycle() {
        System.out.print("Cycle");
        for (int i = 0; i < path.size(); i++) {
            System.out.print("[" + path.get(i).i + "," + path.get(i).j + "]");
        }
        System.out.println();
    }

    public static void addEmptyBasisCell() {
        while(true) {
            int i = (int) Math.round(Math.random() * (M-1));
            int j = (int) Math.round(Math.random() * (N-1));
            if(BasisCells.get(i,j) == 1) {
                continue;
            }
            if(D.get(i,j) != 0) {
                System.out.println("Ненулевые отгрузки для не базисной ячейки");
                break;
            }
            BasisCells.set(i, j, 1);
            break;
        }
    }

    public static void main(String[] args) {
        metodSeveroZapad(); // Первоначальный план
        System.out.println("Первоначальный план методом северо-западного угла");
        System.out.println(D);
        System.out.println("Стоимость перевозки: " + getPrice());

        while(true) {
            int BasisCellsCount = 0;
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    if (D.get(i, j) > 0) {
                        BasisCells.set(i, j, 1);
                        BasisCellsCount++;
                    } else if (D.get(i, j) < 0) {
                        System.out.println("Отгрузки не должны быть отрицательными");
                    }
                    else
                        BasisCells.set(i, j, 0);
                }
            }

            while (BasisCellsCount < (M + N - 1)) {
                System.out.println("Решение вырождено");
                addEmptyBasisCell();
                BasisCellsCount++;
            }

            if (!getPotencial())
                continue;
            //System.out.println("U: " + U);
            //System.out.println("V: " + V);
            if (isOptimum()) {
                System.out.println(D);
                System.out.println("решение оптимально.");
                break;
            }
            System.out.println("решение неоптимально.");

            System.out.println("Delta: ["+ maxI + "," + maxJ + "]");
            if(!getCycle(maxI, maxJ)) {
                System.out.println("Не удалось найти цикл");
                break;
            }
            //printCycle();
            refreshPlan();
            //System.out.println(D);
            //System.out.println(BasisCells);

            System.out.println("Стоимость перевозки: " + getPrice());
            System.out.println();
        }
    }
}
