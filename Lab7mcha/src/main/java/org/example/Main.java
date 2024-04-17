package org.example;

import java.util.Scanner;

public class Main {
    public static void main(String[] args) {
        Scanner scanner =new Scanner(System.in);
        int n = scanner.nextInt();
        int m = scanner.nextInt();
        int [][] matrix = new int[n][m];
        for (int i = 0 ; i < n ; i++)
        {
            for (int j =0 ; j <m;j++)
            {
                matrix[i][j]=scanner.nextInt();
            }
        }
        for(int i = 0;i <n ;i++)
        {
            for (int j = m--; j>=0;j--)
            {
                System.out.print(matrix[i][j]+" ");
            }
            System.out.println();
        }
    }
}
