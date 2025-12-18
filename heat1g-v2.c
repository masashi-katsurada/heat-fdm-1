/*
 * heat1g-v2.c 1次元熱方程式, 境界条件は Robin (含む Neumann) と Dirichlet
 *
 * Robin境界条件に取り組む機会を与えてくれた A.Y. に感謝
 * 2005年3月10日 作成
 * 2006年1月2日  注釈を書き足す
 * 2024年8月21日 URLなど注釈を書き直す
 *
 *   u_t(x,t)=u_{xx}(x,t) ((x,t)∈(0,1)×(0,∞))
 *   -u_x(0,t)+γu(0,t)=T0  (t∈(0,∞)) または u(0,t)=α (t∈(0,∞))
 *    u_x(1,t)+γu(1,t)=T1  (t∈(0,∞)) または u(1,t)=β (t∈(0,∞))
 *   u(x,0)=u0(x)         (x∈[0,1])
 *
 *  差分方程式の説明は次の文書に書き足した (2024年8月21日時点で第1章5節)。
 *  『熱方程式に対する差分法 I --- 区間における熱方程式 ---』
 *   https://m-katsurada.sakura.ne.jp/labo/text/heat-fdm-1.pdf
 *
 *  このプログラム自体は
 *   curl -O https://m-katsurada.sakura.ne.jp/program/fdm/heat1g-v2.c
 *  とすれば入手できる
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double exact(double x, double t);
double u0(double x);
void trilu(int n, double *al, double *ad, double *au);
void trisol(int n, double *al, double *ad, double *au, double *b);

/* left, right は != 0 だと Robin 境界条件
 *
 *  ∂u
 *  ---- + γ u=T
 *  ∂n
 *
 * を課す (γ=0だとNeumann境界条件)。== 0 だとDirichlet境界条件
 *
 *  u=α (x=0), u=β (x=1)
 *
 * を課す。
 */

int left = 0, right = 0;
double pi;

int main(void)
{
  int i,Nx,n,nMax;
  double *u, *uu;
  double *ad,*al,*au;
  double x,h,lambda,theta,tau,gamma,t,error,maxerror;
  double T0,T1,alpha,beta;
  double Tmax = 0.1;

  pi = 4 * atan(1.0);
  gamma = 1;

  if (left)
    T0 = 0;
  else
    alpha = 0;
  if (right)
    T1 = 0;
  else
    beta = 0;

  printf("Nx="); scanf("%d", &Nx);
  h = 1.0 / Nx;
  lambda = 0.5;
  theta = 0.5;
  tau = lambda * h * h;

  u = malloc(sizeof(double) * (Nx+1));
  uu = malloc(sizeof(double) * (Nx+1));
  ad = malloc(sizeof(double) * (Nx+1));
  au = malloc(sizeof(double) * (Nx+1));
  al = malloc(sizeof(double) * (Nx+1));

  /* 係数行列を作る */
  for (i = 1; i < Nx; i++) {
    ad[i] = 1 + 2 * theta * lambda;
    au[i] = al[i] = - theta * lambda;
  }

  if (left) {
    ad[0] = 1 + 2 * theta * lambda * (1 + gamma * h);
    au[0] = - 2 * theta * lambda;
  }
  else {
    ad[0] = 1; au[0] = 0;
  }
  if (right) {
    ad[Nx] = 1 + 2 * theta * lambda * (1 + gamma * h);
    al[Nx] = - 2 * theta * lambda;
  }
  else {
    ad[Nx] = 1; al[Nx] = 0;
  }
  /* LU 分解 */
  trilu(Nx + 1, al, ad, au);

  /* 初期条件 */
  for (i = 0; i <= Nx; i++)
    u[i] = u0(i * h);

  /* 繰り返し */
  nMax = Tmax / tau + 0.5;
  for (n = 1; n <= nMax; n++) {
    for (i = 1; i < Nx; i++)
      uu[i] = (1-2*(1-theta)*lambda)*u[i]+(1-theta)*lambda*(u[i+1]+u[i-1]);
    if (left)
      uu[0] = (1 - 2 * (1 - theta) * lambda * (1 + gamma * h)) * u[0]
               + 2 * (1 - theta) * lambda * u[1]
               + 2 * h * lambda * T0;
    else
      uu[0] = U0;
    if (right)
      uu[Nx]  = (1 - 2 * (1 - theta) * lambda * (1 + gamma * h)) * u[Nx]
                 + 2 * (1 - theta) * lambda * u[Nx - 1]
                 + 2 * h * lambda * T1;
    else
      uu[Nx] = U1;

    trisol(Nx + 1, al, ad, au, uu);
    
    for (i = 0; i <= Nx; i++)
      u[i] = uu[i];
    t = n * tau;
#ifdef VERBOSE
    printf("%f\n", t);
    for (i = 0; i <= Nx; i++)
      printf("%5.2f ", u[i]);
    printf("\n");
    for (i = 0; i <= Nx; i++)
      printf("%5.2f ", exact(i * h, t));
    printf("\n");
#endif
    maxerror = 0;
    for (i = 0; i <= Nx; i++) {
      error = fabs(exact(i * h, t) - u[i]);
      if (error > maxerror)
        maxerror = error;
    }
    printf("error=%e\n", maxerror);
  }
}

#define lambda1 2.028757838110434
#define lambda2 4.913180439434883

double exact(double x, double t)
{
  if (left) {
    /* まだ厳密解を知りません */
  }
  else {
    if (right)
      return sin(lambda1 * x) * exp(- lambda1 * lambda1 * t)
        + sin(lambda2 * x) * exp(- lambda2 * lambda2 * t);
    else
      return sin(pi * x) * exp(- pi * pi * t)
        + sin(2 * pi * x) * exp(- 4 * pi * pi * t);
  }
}

double u0(double x)
{
  return exact(x, 0.0);
}

/* 3項方程式 (係数行列が三重対角である連立1次方程式のこと) Ax=b を解く
 *
 *   入力
 *     n: 未知数の個数
 *     al,ad,au: 連立1次方程式の係数行列
 *       (al: 対角線の下側 i.e. 下三角部分 (lower part)
 *        ad: 対角線       i.e. 対角部分 (diagonal part)
 *        au: 対角線の上側 i.e. 上三角部分 (upper part)
 *        つまり
 *
 *         ad[0] au[0]   0   .................. 0
 *         al[1] ad[1] au[1]   0  ............. 0
 *           0   al[2] ad[2] au[2]  0 ......... 0
 *                         ....................
 *                         al[n-2] ad[n-2] au[n-2]
 *                           0     al[n-1] ad[n-1]

 *        al[i] = A_{i,i-1}, ad[i] = A_{i,i}, au[i] = A_{i,i+1},
 *        al[0], au[n-1] は意味がない)
 *
 *     b: 連立1次方程式の右辺の既知ベクトル
 *        (添字は 0 から。i.e. b[0],b[1],...,b[n-1] にデータが入っている。)
 *   出力
 *     al,ad,au: 入力した係数行列を LU 分解したもの
 *     b: 連立1次方程式の解
 *   能書き
 *     一度 call すると係数行列を LU 分解したものが返されるので、
 *     以後は同じ係数行列に関する連立1次方程式を解くために、
 *     関数 trisol() が使える。
 *   注意
 *     Gauss の消去法を用いているが、ピボットの選択等はしていな
 *     いので、ピボットの選択をしていないので、係数行列が正定値である
 *     などの適切な条件がない場合は結果が保証できない。
 */

/* 三重対角行列の LU 分解 (pivoting なし) */
void trilu(int n, double *al, double *ad, double *au)
{
  int i, nm1 = n - 1;
  /* 前進消去 (forward elimination) */
  for (i = 0; i < nm1; i++) {
    al[i + 1] /= ad[i];
    ad[i + 1] -= au[i] * al[i + 1];
  }
}

/* LU 分解済みの三重対角行列を係数に持つ3項方程式を解く */
void trisol(int n, double *al, double *ad, double *au, double *b)
{
  int i, nm1 = n - 1;
  /* 前進消去 (forward elimination) */
  for (i = 0; i < nm1; i++) b[i + 1] -= b[i] * al[i + 1];
  /* 後退代入 (backward substitution) */
  b[nm1] /= ad[nm1];
  for (i = n - 2; i >= 0; i--) b[i] = (b[i] - au[i] * b[i + 1]) / ad[i];
}
