/*
 * heat2n-i-naive2.c --- 2次元熱方程式 (Neumann 境界条件) を陰解法で素朴に解く
 *
 * http://nalab.mind.meiji.ac.jp/~mk/program/fdm/heat2n-i-naive2.c
 * http://nalab.mind.meiji.ac.jp/~mk/program/fdm/smallmatrix.h
 * http://nalab.mind.meiji.ac.jp/~mk/program/linear/lu.c
 * http://nalab.mind.meiji.ac.jp/~mk/program/linear/lu.h
 *
 * To compile
 *     cglsc heat2n-i-naive.c lu.c
 *
 * このプログラムについては 1998 年度卒研の学生だった深石君に感謝します。
 * このプログラムは効率的ではないので、N_x, N_y はあまり大きくしないこと。
 * 50程度にしておくことを勧めます。
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

// 格子点に番号をつける。
// いわゆる column major mode である。
// m = N_x + 1
// 行列の最初のm行, 最後のm行は j=0,N_y (長方形領域の下辺、上辺) に対応する。
#define phi(i,j) (j)*m+(i)

// Neumann境界値
double b_t(double, double, double);
double b_b(double, double, double);
double b_l(double, double, double);
double b_r(double, double, double);

int main(void)
{
    double a, b, c, d;
    int N_x, N_y, m, N, i, j, p, q, L, n, nMax;
    matrix Uk, A;
    double *B, *vector_U, cond;
    int *iwork, skip;
    double h_x, h_y, lambda_x, lambda_y, lambda, lambda_limit, tau, Tmax, dt;
    double f(double, double), Phi(double, double, double);
    double theta, gamma, alpha, beta, gamma_p, alpha_p, beta_p;
    double x, y, tnp1, tn;

    /* 問題を考える区間 [a,b]×[c,d] */
    a = 0.0; b = 1.0; c = 0.0; d = 1.0;

    /* 区間の分割数 */
    printf("Nx, Ny: "); scanf("%d %d", &N_x, &N_y);

    m = N_x + 1;
    N = (N_x + 1) * (N_y + 1);
    /* 空間の刻み幅 */
    h_x = (b - a) / N_x;
    h_y = (d - c) / N_y;

    /* 行列、ベクトルを記憶する変数のメモリー割り当て */
    if ((Uk = new_matrix(N_x + 1, N_y + 1)) == NULL) {
        fprintf(stderr, "数列 U^k を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((A = new_matrix(N, N)) == NULL) {
        fprintf(stderr, "係数行列 A を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((B = malloc(sizeof(double) * N)) == NULL) {
        fprintf(stderr, "B を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((vector_U = malloc(sizeof(double) * N)) == NULL) {
        fprintf(stderr, "vector_U を記憶する領域の確保に失敗\n");
        exit(1);
    }
    if ((iwork = malloc(sizeof(int) * N)) == NULL) {
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

    /* 連立1次方程式に現れる係数 
       γU_{ij}^{n+1}
         +α(U_{i+1,j}^{n+1}+U_{i-1,j}^{n+1})
         +β(U_{i,j+1}^{n+1}+U_{i,j-1}^{n+1})
      =γ'U_{ij}^n
         +α'(U_{i+1,j}^n+U_{i-1,j}^n)
         +β'(U_{i,j+1}^n+U_{i,j-1}^n)
     と書いたときのα,β,γ,α',β',γ' */
    gamma = 1.0 + 2.0 * theta * lambda;
    alpha = - theta * lambda_x;
    beta  = - theta * lambda_y;

    gamma_p = 1.0 - 2.0 * (1.0 - theta) * lambda;
    alpha_p = (1.0 - theta) * lambda_x;
    beta_p  = (1.0 - theta) * lambda_y;

    /* 係数行列の作成 */
    /* まず 0 クリア */
    for (p = 0; p < N; p++) {
        for (q = 0; q < N; q++)
            A[p][q] = 0.0;
    }
    for (i = 0; i <= N_x; i++)
        for (j = 0; j <= N_y; j++) {
            L = phi(i, j);
            if (j != 0)   A[L][L - m] = beta;
            if (i != 0)   A[L][L - 1] = alpha;
                          A[L][L    ] = gamma;
            if (i != N_x) A[L][L + 1] = alpha;
            if (j != N_y) A[L][L + m] = beta;
        }

    /* 左辺または右辺 (角点を除く) */
    for (j = 1; j < N_y; j++) {
        /* 左辺にある (0,j) */
        L = phi(0,j);
#ifdef ORIGINAL
        // 右の格子点の係数を2倍
        A[L][L + 1] *= 2.0;
#else
        // その格子点と、上と下の格子点の係数を2で割る。
        // 左の格子点(仮想格子点)の係数は行列にない。
        A[L][L - m] /= 2.0; A[L][L] /= 2.0; A[L][L + m] /= 2.0;
#endif
        /* 右辺にある (N_x,j) */
        L = phi(N_x,j);
#ifdef ORIGINAL
        // 左の格子点の係数を2倍
        A[L][L - 1] *= 2.0;
#else
        // その格子点と、上と下の格子点の係数を2で割る。
        // 右の格子点(仮想格子点)の係数は行列にない。
        A[L][L - m] /= 2.0; A[L][L] /= 2.0; A[L][L + m] /= 2.0;
#endif
    }

    /* 下辺または上辺 (角点を除く) */
    for (i = 1; i < N_x; i++) {
        /* 下辺にある (i,0) */
        L = phi(i,0);
#ifdef ORIGINAL
        // 上の格子点の係数を2倍
        A[L][L + m] *= 2.0;
#else
        // その格子点と、左と右の格子点の係数を2で割る。
        // 下の格子点(仮想格子点)の係数は行列にない。
        A[L][L - 1] /= 2.0; A[L][L] /= 2.0; A[L][L + 1] /= 2.0;
#endif
        /* 上辺にある (i,N_y) */
        L = phi(i,N_y);
#ifdef ORIGINAL
        // 下の格子点の係数を2倍
        A[L][L - m] *= 2.0;
#else
        // その格子点と、左と右の格子点の係数を2で割る。
        // 上の格子点(仮想格子点)の係数は行列にない。
        A[L][L - 1] /= 2.0; A[L][L] /= 2.0; A[L][L + 1] /= 2.0;
#endif
    }
    /* 角の点 */
#ifdef ORIGINAL
    /* 左下の角点 */
    L = phi(0,0);
    A[L][L+1] *= 2.0; A[L][L+m] *= 2.0;
    /* 左上の角点 */
    L = phi(0,N_y);
    A[L][L+1] *= 2.0; A[L][L-m] *= 2.0;
    /* 右下の角点 */
    L = phi(N_x,0);
    A[L][L-1] *= 2.0; A[L][L+m] *= 2.0;
    /* 右上の角点 */
    L = phi(N_x,N_y);
    A[L][L-1] *= 2.0; A[L][L-m] *= 2.0;
#else
    /* 左下の角点 -- 右と上の格子点の係数を2で、その格子点の係数を4で割る */
    L = phi(0,0);
    A[L][L] /= 4.0; A[L][L+1] /= 2.0; A[L][L+m] /= 2.0;
    /* 左上の角点 -- 右と下の格子点の係数を2で、その格子点の係数を4で割る */
    L = phi(0,N_y);
    A[L][L] /= 4.0; A[L][L+1] /= 2.0; A[L][L-m] /= 2.0;
    /* 右下の角点 -- 左と上の格子点の係数を2で、その格子点の係数を4で割る */
    L = phi(N_x,0);
    A[L][L] /= 4.0; A[L][L-1] /= 2.0; A[L][L+m] /= 2.0;
    /* 右上の角点 -- 左と下の格子点の係数を2で、その格子点の係数を4で割る */
    L = phi(N_x,N_y);
    A[L][L] /= 4.0; A[L][L-1] /= 2.0; A[L][L-m] /= 2.0;
#endif

    /* 連立1次方程式の係数行列を表示する */
    if (N < 20) {
      printf("素朴に作った連立1次方程式の行列\n");
      for (p = 0; p < N; p++) {
        for (q = 0; q < N; q++)
          printf("%7.3g", A[p][q]);
        printf("\n");
      }
    }

#ifndef ORIGINAL
    /* 対称性のチェック */
    for (p = 0; p < N; p++)
      for (q = 0; q <= p; q++)
        if (A[p][q] != A[q][p])
          printf("A[%d][%d]=%g, A[%d][%d]=%g\n", p, q, A[p][q], q, p, A[q][p]);
#endif

    printf("備考: 1+2θλ=%5.2f, -θλx=%5.2f, -θλy=%5.2f\n",
           gamma, alpha, beta);

    printf("Tmax: "); scanf("%lf", &Tmax);
    printf("Δt: ");  scanf("%lf", &dt);
    skip = rint(dt / tau);
    if (skip == 0) skip = 1;

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
    decomp(N, A, &cond, iwork, B);
    if (cond + 1 == cond) {
        /* 条件数が大きければ、計算をあきらめる */
        printf("MATRIX IS SINGULAR TO WORKING PRECISION\n");
        return 0;
    }

    /* 時間に関するループ */
    for (n = 0; n < nMax; n++) {
        tn = n * tau; tnp1 = (n + 1) * tau; 

        /* まず、素朴な連立1次方程式の右辺を用意する */
        /* 内部の格子点 */
        for (i = 1; i < N_x; i++)
            for (j = 1; j < N_y; j++) {
                L = phi(i,j);
                B[L] = gamma_p * Uk[i][j]
                    + alpha_p * (Uk[i + 1][j] + Uk[i - 1][j])
                    + beta_p * (Uk[i][j + 1] + Uk[i][j - 1]);
            }
        /* 以下、境界にある格子点での方程式を立てる。
         * 仮想格子点での値は境界条件を中心差分近似した方程式を
         * 用いて消去する
         *   ……右辺にも移項する量があるので、右辺について後で処理する */
        /* 下の辺、上の辺にある格子点 (角の点は含めない) */
        for (i = 1; i < N_x; i++) {
            x = a + i * h_x;
            /* (i, 0) */
            L = phi(i,0); y = c;
            B[L] = gamma_p * Uk[i][0]
                 + alpha_p * (Uk[i + 1][0] + Uk[i - 1][0])
                 + 2 * beta_p * Uk[i][1]
                 + 2 * h_y  * (beta_p * b_b(x, y, tn) - beta * b_b(x, y, tnp1));
#ifndef ORIGINAL
            B[L] /= 2;
#endif
            /* (i, N_y) */
            L = phi(i,N_y); y = d;
            B[L] = gamma_p * Uk[i][N_y]
                 + alpha_p * (Uk[i + 1][N_y] + Uk[i - 1][N_y])
                 + 2 * beta_p * Uk[i][N_y - 1]
                 + 2 * h_y * (beta_p * b_t(x, y, tn) - beta * b_t(x, y, tnp1));
#ifndef ORIGINAL
            B[L] /= 2;
#endif
        }
        /* 左の辺、右の辺にある格子点 (角の点は含めない) */
        for (j = 1; j < N_y; j++) {
            y = c + j * h_y;
            /* (0, j) */
            L = phi(0,j); x = a;
            B[L] = gamma_p * Uk[0][j]
                 + 2 * alpha_p * Uk[1][j]
                 + beta_p * (Uk[0][j + 1] + Uk[0][j - 1])
                 + 2 * h_x * (alpha_p * b_l(x, y, tn) - alpha * b_l(x, y, tnp1));
#ifndef ORIGINAL
            B[L] /= 2;
#endif
            /* (N_x, j) */
            L = phi(N_x,j); x = b;
            B[L] = gamma_p * Uk[N_x][j]
                 + 2 * alpha_p * Uk[N_x - 1][j]
                 + beta_p * (Uk[N_x][j + 1] + Uk[N_x][j - 1])
                 + 2 * h_x * (alpha_p * b_r(x, y, tn) - alpha * b_r(x, y, tnp1));
#ifndef ORIGINAL
            B[L] /= 2;
#endif
        }
        /* 左下 */
        L = phi(0,0); x = a; y = c;
        B[L] = gamma_p * Uk[0][0]
             + 2 * alpha_p * Uk[1][0] + 2 * beta_p * Uk[0][1]
             + 2 * h_x * (alpha_p * b_l(x,y,tn) - alpha * b_l(x,y,tnp1))
             + 2 * h_y * (beta_p * b_b(x,y,tn) - beta * b_b(x,y,tnp1));
#ifndef ORIGINAL
            B[L] /= 4;
#endif
        /* 左上 */
        L = phi(0,N_y); x = a; y = d;
        B[L] = gamma_p * Uk[0][N_y]
             + 2 * alpha_p * Uk[1][N_y] + 2 * beta_p * Uk[0][N_y - 1]
             + 2 * h_x * (alpha_p * b_l(x,y,tn) - alpha * b_l(x,y,tnp1))
             + 2 * h_y * (beta_p * b_t(x,y,tn) - beta * b_t(x,y,tnp1));
#ifndef ORIGINAL
            B[L] /= 4;
#endif
        /* 右下 */
        L = phi(N_x,0); x = b; y = c;
        B[L] = gamma_p * Uk[N_x][0]
             + 2 * alpha_p * Uk[N_x - 1][0] + 2 * beta_p * Uk[N_x][1]
             + 2 * h_x * (alpha_p * b_r(x,y,tn) - alpha * b_r(x,y,tnp1))
             + 2 * h_y * (beta_p * b_b(x,y,tn) - beta * b_b(x,y,tnp1));
#ifndef ORIGINAL
            B[L] /= 4;
#endif
        /* 右上 */
        L = phi(N_x,N_y); x = b; y = d;
        B[L] = gamma_p * Uk[N_x][N_y]
             + 2 * alpha_p * Uk[N_x - 1][N_y] + 2 * beta_p * Uk[N_x][N_y - 1]
             + 2 * h_x * (alpha_p * b_r(x,y,tn) - alpha * b_r(x,y,tnp1))
             + 2 * h_y * (beta_p * b_t(x,y,tn) - beta * b_t(x,y,tnp1));
#ifndef ORIGINAL
            B[L] /= 4;
#endif

        /* A vector_U = B を解く */
        solve(N, A, B, iwork);
        /* */
        for (i = 0; i <= N_x; i++)
            for (j = 0; j <= N_y; j++)
                Uk[i][j] = B[phi(i,j)];

        /* データを数値で表示 */
        if (n % skip == 0) {
#ifdef  PRINT
            for (i = 0; i <= N_x; i++) {
                for (j = 0; j <= N_y; j++)
                    printf(" %6.2f", Uk[i][j]);
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

/* Neumann 境界値 */
/* 上の辺での Neumann 境界値 ∂u/∂y (boundary value on the top side */
double b_t(double x, double y, double t)
{
  return 0.0; /* 同次 Neumann 境界条件 */
}

/* 下の辺での Neumann 境界値 -∂u/∂y (boundary value on the bottom side */
double b_b(double x, double y, double t)
{
  return 0.0; /* 同次 Neumann 境界条件 */
}

/* 左の辺での Neumann 境界値 ∂u/∂y (boundary value on the left side */
double b_l(double x, double y, double t)
{
  return 0.0; /* 同次 Neumann 境界条件 */
}

/* 右の辺での Neumann 境界値 -∂u/∂y (boundary value on the right side */
double b_r(double x, double y, double t)
{
  return 0.0; /* 同次 Neumann 境界条件 */
}
