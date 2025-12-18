/*
 * heat2d-i-naive2.c --- 2次元熱方程式を陰解法で素朴に解く (version 3)
 *                       境界にある格子点における値は未知数に含めない方法
 * http://nalab.mind.meiji.ac.jp/~mk/program/fdm/heat2d-i-naive.c
 * http://nalab.mind.meiji.ac.jp/~mk/program/fdm/smallmatrix.h
 * http://nalab.mind.meiji.ac.jp/~mk/program/linear/lu.c
 * http://nalab.mind.meiji.ac.jp/~mk/program/linear/lu.h
 * To compile
 *     cc -O3 -c lu.c
 *     cglsc heat2d-i-naive2.c lu.o
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef OLD
#include <matrix.h>
#else
#include "smallmatrix.h"
#endif
#ifndef G_DOUBLE
#define G_DOUBLE
#endif
#include <glsc.h>
#include "lu.h"

#define psi(i,j) (((j)-1)*mm+(i)-1)

int main(void)
{
    double a, b, c, d;
    int N_x, N_y, mm, NN, i, j, p, q, L, n, nMax;
    matrix Uk, A;
    double *B, *vector_U, cond;
    int *iwork, skip;
    double h_x, h_y, lambda_x, lambda_y, lambda, lambda_limit, tau, Tmax, dt;
    double f(double, double), alpha(double, double);
    double theta, beta_0, beta_1, beta_2, beta_00, beta_10, beta_20;
    double x, y, t;

    /* 問題を考える区間 [a,b]×[c,d] */
    a = 0.0; b = 1.0; c = 0.0; d = 1.0;

    /* 区間の分割数 */
    printf("Nx, Ny: "); scanf("%d %d", &N_x, &N_y);

    mm = N_x - 1;
    NN = (N_x - 1) * (N_y - 1);
    /* 空間の刻み幅 */
    h_x = (b - a) / N_x;
    h_y = (d - c) / N_y;

    /* 行列、ベクトルを記憶する変数のメモリー割り当て */
    if ((Uk = new_matrix(N_x + 1, N_y + 1)) == NULL) {
        fprintf(stderr, "数列 U^k を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((A = new_matrix(NN, NN)) == NULL) {
        fprintf(stderr, "係数行列 A を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((B = malloc(sizeof(double) * NN)) == NULL) {
        fprintf(stderr, "B を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((vector_U = malloc(sizeof(double) * NN)) == NULL) {
        fprintf(stderr, "vector_U を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((iwork = malloc(sizeof(int) * NN)) == NULL) {
        fprintf(stderr, "iwork を記憶する領域の確保に失敗\n");
        exit(1);
    }

    /* θ法の重みの決定 */
    printf("θ (0≦θ≦1): "); scanf("%lf", &theta);

    if (theta == 1.0) {
        printf("τ: "); scanf("%lf", &tau);
    } else {
        printf("τ(≦%g≡最大値ノルムに関する安定性条件を満たすτの上限): ",
               0.5 / (1 - theta) / (1 / (h_x * h_x) + 1 / (h_y * h_y)));
        scanf("%lf", &tau);
    }

    lambda_x = tau / (h_x * h_x);
    lambda_y = tau / (h_y * h_y);
    lambda = lambda_x + lambda_y;

    /* 最大値ノルムに関する安定性を満たすλの上限 */
    lambda_limit = 1.0 / (2.0 * (1.0 - theta));

    if (lambda > lambda_limit)
        printf("注意: λ=%g>1/2(1-θ) となっています。\n", lambda);
    else
        printf("λ=%g\n", lambda);

    /* 初期値の設定 */
    for (i = 0; i <= N_x; i++) {
      x = a + i * h_x;
      for (j = 0; j <= N_y; j++)
        Uk[i][j] = f(x, c + j * h_y);
    }

    /* 連立1次方程式に現れる係数 */
    beta_0 = 1.0 + 2.0 * theta * lambda;
    beta_1 = -theta * lambda_x;
    beta_2 = -theta * lambda_y;

    beta_00 = 1.0 - 2.0 * (1.0 - theta) * lambda;
    beta_10 = (1.0 - theta) * lambda_x;
    beta_20 = (1.0 - theta) * lambda_y;

    /* 係数行列の作成 */
    /* まず 0 クリア */
    for (p = 0; p < NN; p++) {
        for (q = 0; q < NN; q++)
            A[p][q] = 0.0;
    }
    for (i = 1; i < N_x; i++)
        for (j = 1; j < N_y; j++) {
            L = psi(i, j);
            if (j != 1) A[L][L - mm] = beta_2;
            if (i != 1) A[L][L -  1] = beta_1;
            A[L][L] = beta_0;
            if (i != N_x - 1) A[L][L +  1] = beta_1;
            if (j != N_y - 1) A[L][L + mm] = beta_2;
        }

    /* 連立1次方程式の係数行列を表示する */
    if (NN < 20) {
      printf("連立1次方程式の行列\n");
      for (p = 0; p < NN; p++) {
        for (q = 0; q < NN; q++)
          printf(" %4.1f", A[p][q]);
        printf("\n");
      }
    }

    printf("備考: 1+2θλ=%4.1f, -θλx=%5.1f, -θλy=%5.1f\n",
           beta_0, beta_1, beta_2);

    printf("Tmax: "); scanf("%lf", &Tmax);
    printf("Δt: ");  scanf("%lf", &dt);
    skip = rint(dt / tau);
    if (skip == 0) {
        skip = 1;
    }
    dt = skip * tau;

    nMax = rint(Tmax / tau);

    /* グラフィックス・ライブラリィ GLSC の呼び出し */
    g_init("Meta", 250.0, 160.0);
    g_device(G_BOTH);
    g_def_scale(0,
                0.0, 1.0, 0.0, 1.0,
                30.0, 70.0, 100.0, 72.0);
    g_def_scale(4,
                -1.0, 1.0, -1.0, 1.0,
                30.0, 30.0, 100.0, 100.0);
    g_def_line(0, G_BLACK, 0, G_LINE_SOLID);
    g_sel_scale(0);

    g_cls();
#ifdef OLD
    g_hidden2(1.0, 1.0, 0.4, -1.0, 1.0, 5.0, 25.0, 20.0, 20.0, 20.0,
              150.0, 100.0, Uk, N_x + 1, N_y + 1,
              1, G_SIDE_NONE, 2, 1);
#else
    g_hidden(1.0, 1.0, 0.4, -1.0, 1.0, 5.0, 25.0, 20.0, 20.0, 20.0,
             150.0, 100.0, &Uk[0][0], N_x + 1, N_y + 1,
             1, G_SIDE_NONE, 2, 1);
#endif
    /* 係数行列 LU 分解 */
    decomp(NN, A, &cond, iwork, B);
    if (cond + 1 == cond) {
        /* 条件数が大きければ、計算をあきらめる */
        printf("MATRIX IS SINGULAR TO WORKING PRECISION\n");
        return 0;
    }

    /* 時間に関するループ */
    for (n = 1; n <= nMax; n++) {

        /* まず、素朴な連立1次方程式の右辺を用意する */
        /* 第 n ステップにおける値の寄与 */
        for (i = 1; i < N_x; i++)
            for (j = 1; j < N_y; j++) {
                L = psi(i,j);
                B[L] = beta_00 * Uk[i][j]
                    + beta_10 * (Uk[i + 1][j] + Uk[i - 1][j])
                    + beta_20 * (Uk[i][j + 1] + Uk[i][j - 1]);
            }
        /* 第 n+1 ステップの境界上の格子点における値 (既知) の寄与 */
        for (j = 1; j < N_y; j++) {
            y = c + j * h_y;
            /* (1,j) のとき */
            L = psi(1, j);
            B[L] -= beta_1 * alpha(a, y);
            /* (N_x-1,j) のとき */
            L = psi(N_x - 1, j);
            B[L] -= beta_1 * alpha(b, y);
        }
        for (i = 1; i < N_x; i++) {
            x = a + i * h_x;
            /* (i,1) のとき */
            L = psi(i, 1);
            B[L] -= beta_2 * alpha(x, c);
            /* (i,N_y-1) のとき */
            L = psi(i, N_y - 1);
            B[L] -= beta_2 * alpha(x, d);
        }

        /* A vector_U = B を解く */
        solve(NN, A, B, iwork);
        /* */
        for (i = 1; i < N_x; i++)
            for (j = 1; j < N_y; j++)
                Uk[i][j] = B[psi(i,j)];

        /* 下の辺、上の辺にある格子点 (角の点も含める) */
        for (i = 0; i <= N_x; i++) {
            x = a + i * h_x;
            Uk[i][0]   = alpha(x, c);
            Uk[i][N_y] = alpha(x, d);
        }
        /* 左の辺、右の辺にある格子点 (角の点は含めない) */
        for (j = 1; j < N_y; j++) {
            y = c + j * h_y;
            Uk[0][j]   = alpha(a, y);
            Uk[N_x][j] = alpha(b, y);
        }

        /* データを数値で表示 */
        if (n % skip == 0) {
#ifdef  PRINT
            for (i = 0; i <= N_x; i++) {
                for (j = 0; j <= N_y; j++)
                    printf(" %5.1f", Uk[i][j]);
                printf("\n");
            }
#endif
            /* 鳥瞰図を描く */
            g_cls();
#ifdef OLD
            g_hidden2(1.0, 1.0, 0.4, -1.0, 1.0, 5.0, 25.0, 20.0, 20.0, 20.0,
                      150.0, 100.0, Uk, N_x + 1, N_y + 1,
                      1, G_SIDE_NONE, 2, 1);
#else
            g_hidden(1.0, 1.0, 0.4, -1.0, 1.0, 5.0, 25.0, 20.0, 20.0, 20.0,
                     150.0, 100.0, &Uk[0][0], N_x + 1, N_y + 1,
                     1, G_SIDE_NONE, 2, 1);
#endif
        }
    }
    /* マウスでクリックされるのを待つ */
    g_sleep(-1.0);
    /* ウィンドウを消す */
    g_term();

    return 0;
}

/* 初期値 */
double f(double x, double y)
{
    /* ピラミッド型の関数 */
    if (y > 0.5)
        y = 1 - y;
    if (x > 0.5)
        x = 1 - x;
    if (y < x)
        return 5 * y;
    else
        return 5 * x;
}

/* 境界値 */
double alpha(double x, double y)
{
    /* 同次 Dirichlet 境界条件 */
    return 0.0;
}
