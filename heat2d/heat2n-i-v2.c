/*
 * heat2n-i-v2.c --- 2次元熱方程式 (Neumann 境界条件) を陰解法で解く
 *                   Neumann 境界値の扱いを変更
 *
 * http://nalab.mind.meiji.ac.jp/~mk/program/fdm/heat2n-i-v2.c
 * http://nalab.mind.meiji.ac.jp/~mk/program/fdm/smallmatrix.h
 * http://nalab.mind.meiji.ac.jp/~mk/program/linear/symbaldlu.c
 * http://nalab.mind.meiji.ac.jp/~mk/program/linear/symbaldlu.h
 *
 * scp -p heat2n-i-v2.c nalab.mind.meiji.ac.jp:Sites/misc/20240821
 *
 * To compile
 *     cglsc heat2n-i-v2.c symbandlu.c
 * or
 *     ccmg heat2n-i-v2.c symbandlu.c
 *
 *   このプログラムについては 1998 年度卒研の学生だった深石君に感謝します。
 *
 * version 2 での変更点
 *  Neumann 境界値を返す関数をこれまで
 *    double Phi(double x, double y, double t)
 *  としていたのを
 *    double b_t(double x, double y, double t) // -∂u/∂y(x,d,t)
 *    double b_b(double x, double y, double t) // -∂u/∂y(x,c,t)
 *    double b_l(double x, double y, double t) // -∂u/∂x(a,y,t)
 *    double b_r(double x, double y, double t) //  ∂u/∂y(b,y,t)
 *  という4つの関数にした。
 *  引数の数を減らせるが、x,y,t のままにすることを選んだ(1つダミーになるが)。
 *    これまで、例えば左下の角点では
 *
 *      L = phi(0,0);
 *      B[L] = gamma_p * Uk[0][0]
 *           + 2 * alpha_p * Uk[1][0] + 2 * beta_p * Uk[0][1]
 *           + 2 * h_x * (alpha * Phi(a,c,t) - alpha_p * Phi(a,c,t_p))
 *           + 2 * h_y * (beta * Phi(a,c,t) - beta_p * Phi(a,c,t_p));
 *           B[L] /= 4;
 *
 *   としていたのを
 *
 *      L = phi(0,0);
 *      B[L] = gamma_p * Uk[0][0]
 *           + 2 * alpha_p * Uk[1][0] + 2 * beta_p * Uk[0][1]
 *           + h_x * (alpha * b_l(a,c,t) - alpha_p * b_l(a,c,t_p))
 *           + h_y * (beta * b_b(a,c,t) - beta_p * b_b(a,c,t_p));
 *          B[L] /= 4;
 *
 *   と改めた。これまでは(温めすぎ|冷やしすぎ)で、こうするのが正しいと
 *   信じている。
 *
 *     下の辺での "修正"
            // (i, 0)
            L = phi(i,0);
            B[L] = gamma_p * Uk[i][0]
                 + alpha_p * (Uk[i + 1][0] + Uk[i - 1][0])
                 + 2 * beta_p * Uk[i][1]
                 + 2 * h_y  * (beta * b_b(x, c, t) - beta_p * b_b(x, c, t_p));
            B[L] /= 2;
 *     左の辺での "修正"
            // (0, j)
            L = phi(0,j);
            B[L] = gamma_p * Uk[0][j]
                 + 2 * alpha_p * Uk[1][j]
                 + beta_p * (Uk[0][j + 1] + Uk[0][j - 1])
                 + 2 * h_x * (alpha * b_l(a, y, t) - alpha_p * b_l(a, y, t_p));
            B[L] /= 2;

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
#include "symbandlu.h"

// 格子点に番号をつける。
// いわゆる column major mode である。
// m = N_x + 1
// 行列の最初のm行, 最後のm行は j=0,N_y (長方形領域の下辺、上辺) に対応する。
#define phi(i,j) (j)*m+(i)

double sqr(double x) { return x * x; }
double max(double a, double b) { if (a > b) return a; else return b; }
void print_matrix(matrix, int, int);
  
int p_id = -1; //
double pi;

int main(int argc, char **argv)
{
    double a, b, c, d;
    int N_x, N_y, m, N, i, j, p, q, L, n, nMax;
    matrix Uk, A;
    double *B, *vector_U, cond;
    int *iwork, skip;
    double h_x, h_y, lambda_x, lambda_y, lambda, lambda_limit, tau, Tmax, dt;
    double f(double, double);
    double b_t(double, double, double), b_b(double, double, double);
    double b_l(double, double, double), b_r(double, double, double);
    double theta, gamma, alpha, beta, gamma_p, alpha_p, beta_p;
    double x, y, tnp1, tn;
    // 誤差
    double exact_sol(double, double, double);
    void print_problem(int);

    /* 円周率 */
    pi = 4.0 * atan(1.0);

    if (argc == 2) {
       p_id = atoi(argv[1]);
       printf("p_id=%d\n", p_id);
       print_problem(p_id);
    }

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
    if ((A = new_matrix(N, m + 1)) == NULL) {
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
    printf("gamma=%f, alpha=%f, beta=%f\n", gamma, alpha, beta);

    gamma_p = 1.0 - 2.0 * (1.0 - theta) * lambda;
    alpha_p = (1.0 - theta) * lambda_x;
    beta_p  = (1.0 - theta) * lambda_y;
    printf("gamma'=%f, alpha'=%f, beta'=%f\n", gamma_p, alpha_p, beta_p);

    /* 係数行列の作成 */
    /* まず 0 クリア */
    for (p = 0; p < N; p++) {
        for (q = 0; q <= m; q++)
            A[p][q] = 0.0;
    }
    for (i = 0; i <= N_x; i++)
        for (j = 0; j <= N_y; j++) {
            L = phi(i, j);
            /*
            if (j != 0)   A[L][L - m] = beta;
            if (i != 0)   A[L][L - 1] = alpha;
            */
                          A[L][0] = gamma; /* A[L][L    ] = gamma; */
            if (i != N_x) A[L][1] = alpha; /* A[L][L + 1] = alpha; */
            if (j != N_y) A[L][m] = beta;  /* A[L][L + m] = beta;  */
        }

    if (N < 20) {
      printf("連立1次方程式の行列 (対称化前)\n");
      print_matrix(A, N, m);
    }

    /* 左辺または右辺にある格子点 */
    for (j = 1; j < N_y; j++) {
        /* (0,j) */
        L = phi(0,j);
        A[L][0] /= 2.0; A[L][m] /= 2.0;
        /* A[L][L - m] /= 2.0; A[L][L] /= 2.0; A[L][L + m] /= 2.0; */
        /* (N_x,j) */
        L = phi(N_x,j);
        A[L][0] /= 2.0; A[L][m] /= 2.0;
        /* A[L][L - m] /= 2.0; A[L][L] /= 2.0; A[L][L + m] /= 2.0; */
    }

    /* 下辺または上辺にある格子点 */
    for (i = 1; i < N_x; i++) {
        /* (i,0) */
        L = phi(i,0);
        A[L][0] /= 2.0; A[L][1] /= 2.0;
        /* A[L][L - 1] /= 2.0; A[L][L] /= 2.0; A[L][L + 1] /= 2.0; */
        /* (i,N_y) */
        L = phi(i,N_y);
        A[L][0] /= 2.0; A[L][1] /= 2.0;
        /* A[L][L - 1] /= 2.0; A[L][L] /= 2.0; A[L][L + 1] /= 2.0; */
    }
    /* 角の点 (対角成分は4で割る、それ以外は2で割る) */
    /* 左下 */
    L = phi(0,0);
    A[L][0] /= 4.0; A[L][1] /= 2.0; A[L][m] /= 2.0;
    /* A[L][L] /= 4.0; A[L][L+1] /= 2.0; A[L][L+m] /= 2.0; */
    /* 左上 */
    L = phi(0,N_y);
    A[L][0] /= 4.0; A[L][1] /= 2.0;
    /* A[L][L] /= 4.0; A[L][L+1] /= 2.0; A[L][L-m] /= 2.0; */
    /* 右下 */
    L = phi(N_x,0);
    A[L][0] /= 4.0; A[L][m] /= 2.0;
    /* A[L][L] /= 4.0; A[L][L-1] /= 2.0; A[L][L+m] /= 2.0; */
    /* 右上 */
    L = phi(N_x,N_y);
    A[L][0] /= 4.0;
    /* A[L][L] /= 4.0; A[L][L-1] /= 2.0; A[L][L-m] /= 2.0; */

    /* 連立1次方程式の係数行列を表示する */
    if (N < 20) {
      printf("連立1次方程式の行列\n");
      print_matrix(A, N, m);
    }

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
    printf("一時停止します。再開にはウィンドウをクリックして下さい。\n");
    g_sleep(-1.0);

    /* 係数行列 LU 分解 */
    symbandlu(A, N, m + 1);
    if (cond + 1 == cond) {
        /* 条件数が大きければ、計算をあきらめる */
        printf("MATRIX IS SINGULAR TO WORKING PRECISION\n");
        return 0;
    }

    /* 時間に関するループ */
    for (n = 0; n < nMax; n++) {
        tnp1 = (n + 1) * tau; tn = n * tau;

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
                 + 2 * h_y * (beta_p * b_b(x, y, tn) - beta * b_b(x, y, tnp1));
            B[L] /= 2;
            /* (i, N_y) */
            L = phi(i,N_y); y = d;
            B[L] = gamma_p * Uk[i][N_y]
                 + alpha_p * (Uk[i + 1][N_y] + Uk[i - 1][N_y])
                 + 2 * beta_p * Uk[i][N_y - 1]
                 + 2 * h_y * (beta_p * b_t(x, y, tn) - beta * b_t(x, y, tnp1));
            B[L] /= 2;
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
            B[L] /= 2;
            /* (N_x, j) */
            L = phi(N_x,j); x = b;
            B[L] = gamma_p * Uk[N_x][j]
                 + 2 * alpha_p * Uk[N_x - 1][j]
                 + beta_p * (Uk[N_x][j + 1] + Uk[N_x][j - 1])
                 + 2 * h_x * (alpha_p * b_r(x, y, tn) - alpha * b_r(x, y, tnp1));
            B[L] /= 2;
        }
        /* 左下(left,bottom)の角点 */
        L = phi(0,0); x = a; y = c;
        B[L] = gamma_p * Uk[0][0]
             + 2 * alpha_p * Uk[1][0] + 2 * beta_p * Uk[0][1]
             + 2 * h_x * (alpha_p * b_l(x,y,tn) - alpha * b_l(x,y,tnp1))
             + 2 * h_y * (beta_p * b_b(x,y,tn) - beta * b_b(x,y,tnp1));
        B[L] /= 4;
        /* 左上(left,top)の角点 */
        L = phi(0,N_y); x = a; y = d;
        B[L] = gamma_p * Uk[0][N_y]
             + 2 * alpha_p * Uk[1][N_y] + 2 * beta_p * Uk[0][N_y - 1]
             + 2 * h_x * (alpha_p * b_l(x,y,tn) - alpha * b_l(x,y,tnp1))
             + 2 * h_y * (beta_p * b_t(x,y,tn) - beta * b_t(x,y,tnp1));
        B[L] /= 4;
        /* 右下(right,bottom)の角点 */
        L = phi(N_x,0); x = b; y = c;
        B[L] = gamma_p * Uk[N_x][0]
             + 2 * alpha_p * Uk[N_x - 1][0] + 2 * beta_p * Uk[N_x][1]
             + 2 * h_x * (alpha_p * b_r(x,y,tn) - alpha * b_r(x,y,tnp1))
             + 2 * h_y * (beta * b_b(x,y,tn) - beta_p * b_b(x,y,tnp1));
        B[L] /= 4;
        /* 右上(right,top)の角点 */
        L = phi(N_x,N_y); x = b; y = d;
        B[L] = gamma_p * Uk[N_x][N_y]
             + 2 * alpha_p * Uk[N_x - 1][N_y] + 2 * beta_p * Uk[N_x][N_y - 1]
             + 2 * h_x * (alpha_p * b_r(x,y,tn) - alpha * b_r(x,y,tnp1))
             + 2 * h_y * (beta_p * b_t(x,y,tn) - beta * b_t(x,y,tnp1));
        B[L] /= 4;

        /* A vector_U = B を解く */
        symbandsolve(A, B, N, m + 1);
        /* */
        for (i = 0; i <= N_x; i++)
            for (j = 0; j <= N_y; j++)
                Uk[i][j] = B[phi(i,j)];

        /* データを数値で表示 */
        if (n % skip == 0) {
          double errorp, error, sum;
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
            error = 0.0; sum = 0.0;
            for (i = 0; i <= N_x; i++) {
              x = a + i * h_x;
              for (j = 0; j <= N_y; j++) {
                y = c + j * h_y;
                // error += sqr(exact_sol(x, y, t) - Uk[i][j]);
                errorp = sqr(exact_sol(x, y, tnp1) - Uk[i][j]);
                if (errorp > error)
                  error = errorp;
                sum += Uk[i][j];
              }
            }
            printf("error=%e, sum=%e\n", error / ((N_x+1)*(N_y+1)),
                                         h_x * h_y *sum);
        }
    }
    /* マウスでクリックされるのを待つ */
    g_sleep(-1.0);
    /* ウィンドウを消す */
    g_term();

    return 0;
}

void error(int p_id)
{
  fprintf(stderr, "unknown p_id=%d\n", p_id);
  exit(-1);
}

void print_problem(int id)
{
  if (id == 0)
    printf("初期値: cos(πx)cos(πy), Neumann境界値: 0\n");
  else if (id == 1)
    printf("初期値: (x-1/2)^2-(y-1/2)^2+cos(πx)cos(πy), Neumann境界値: 1\n");
  else if (id == 2)
    printf("初期値: (x-1/2)(y-1/2)+cos(πx)cos(πy), Neumann境界値: 1\n");
}

/* 初期値 */
double f(double x, double y)
{
  if (p_id == 0)
    return cos(pi * x) * cos(pi * y);
  else if (p_id == 1)
    return sqr(x - 0.5) - sqr(y - 0.5) + cos(pi * x) * cos(pi * y);
  else if (p_id == 2)
    return (x - 0.5) * (y - 0.5) + cos(pi * x) * cos(pi * y);
  else
    error(p_id);
  return 0.0; // dummy
}

/* 上の辺での Neumann 境界値 ∂u/∂y (boundary value on the top side */
double b_t(double x, double y, double t)
{
  if (p_id == 0)
    return 0.0; /* 同次 Neumann 境界条件 */
  else if (p_id == 1)
    return - 1.0;
  else if (p_id == 2)
    return x - 0.5;
  else
    error(p_id);
  return 0.0; // dummy
}

/* 下の辺での Neumann 境界値 -∂u/∂y (boundary value on the bottom side */
double b_b(double x, double y, double t)
{
  if (p_id == 0)
    return 0.0; /* 同次 Neumann 境界条件 */
  else if (p_id == 1)
    return - 1.0;
  else if (p_id == 2)
    return 0.5 - x;
  else
    error(p_id);
  return 0.0; // dummy
}

/* 左の辺での Neumann 境界値 ∂u/∂y (boundary value on the left side */
double b_l(double x, double y, double t)
{
  if (p_id == 0)
    return 0.0; /* 同次 Neumann 境界条件 */
  else if (p_id == 1)
    return 1.0;
  else if (p_id == 2)
    return 0.5 - y;
  else
    error(p_id);
  return 0.0; // dummy
}

/* 右の辺での Neumann 境界値 -∂u/∂y (boundary value on the right side */
double b_r(double x, double y, double t)
{
  if (p_id == 0)
    return 0.0; /* 同次 Neumann 境界条件 */
  else if (p_id == 1)
    return 1.0;
  else if (p_id == 2)
    return y - 0.5;
  else
    error(p_id);
  return 0.0; // dummy
}

double exact_sol(double x, double y, double t)
{
  if (p_id == 0)
    return exp(- 2 * pi * pi * t) * cos(pi * x) * cos(pi * y);
  else if (p_id == 1)
    return sqr(x - 0.5) - sqr(y - 0.5)
      + exp(- 2.0 * pi * pi * t) * cos(pi * x) * cos(pi * y);
  else if (p_id == 2)
    return (x - 0.5) * (y - 0.5)
      + exp(- 2.0 * pi * pi * t) * cos(pi * x) * cos(pi * y);
  else
    error(p_id);
  return 0.0; // dummy
}

void print_matrix(matrix A, int N, int m)
{
  int i, j;
  double t;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (i - j < - m || i - j > m)
        t = 0;
      else if (i <= j)
        t = A[i][j-i];
      else
        t = A[j][i-j];
      printf("%8.4g ", t);
    }
    printf("\n");
  }
}
