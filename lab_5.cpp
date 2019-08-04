#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <cstring>

const double a = 0;//7;
const double b = 1.2;//20;

double f(double x){
    return 10*x*x*cosh(x)*sin(13*x);//0.053*x*log(fabs(cos(x/2))/3)-1;
}

const double eps_i = 1e-8;
const int N = 30;
const double d_eps = 1e-2;
const int numb = 100;

#define FILE_NAME "table.csv"

double _vector[N];

double phi(double x, int k)
{
    double t0, t1, buff;
    double x2;
    int i;
    if (k == 0) {
        return 1.0;
    }
    if (k == 1) {
        return x;
    }
    t0 = 1.0; t1 = x;
    x2 = x + x;
    for (i = 1; i < k; i++) {
        buff = x2*t1 - t0;
        t0 = t1;
        t1 = buff;
    }
    return t1;
}

double get_Pm(double x, int m, double _vector[N])
{
    double var = 0.0;
    int i;
    for (i = 0; i < m; i++) {
        var += _vector[i]*phi(x, i);
    }
    return var;
}

double Sigma(double x, int m, int)
{
    double var = f(x) - get_Pm(x, m, _vector);
    return var*var;
}

double phi_product(double x, int i, int j)
{
    double var;
    if (i == j) {
        var = phi(x, i);
        return var*var;
    }
    if (j == N) {
        return f(x)*phi(x, i);
    }
    return phi(x, i)*phi(x, j);
}

double Simpson(double a, double b, double h, int n, double (*func)(double,int,int), int i, int j)
{
    double sigma1 = 0.0, sigma2 = 0.0;
    int k;
    for (k = 1; k < n; k += 2) {
        sigma1 += func(a + k*h, i, j);
    }
    for (k = 2; k < n - 1; k += 2) {
        sigma2 += func(a + k*h, i, j);
    }
    return h/3*(func(a, i, j) + func(b, i, j) + 4*sigma1 + 2*sigma2);
}

double integral(double a, double b, double (*func)(double,int,int), int i, int j, double eps)
{
    long long int n = (b - a)/sqrt(sqrt(eps));
    double In = 0.0, I2n = 0.0, h = (b - a)/n;
    n = (b - a)/h;
    if(n % 2 == 1) {
        n++;
    }
    h = (b - a)/n;
    In = Simpson(a, b, h, n, func, i, j);
    n*=2; h /= 2;
    I2n = Simpson(a, b, h, n, func, i, j);
    eps *= 15;
    while (fabs((In - I2n)/I2n) > eps) {
        In = I2n;
        n*=2; h /= 2;
        I2n = Simpson(a, b, h, n, func, i, j);
    }
    return I2n;
}

void get_matrix(double matrix[][N + 1], double a, double b, double eps)
{
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N + 1; j++) {
            matrix[i][j] = integral(a, b, phi_product, i, j, eps);
        }
    }
}

void GaussianElimination(const double _matrix[][N + 1], double _vector[N], const int n)
{
    double **matrix;
    double mx1;
    int i, j, k;
    matrix = new double* [n];
    for (i = 0; i < n; i++) {
        matrix[i] = new double[n + 1];
        memcpy(matrix[i], _matrix[i], n*sizeof(double));
        matrix[i][n] = _matrix[i][N];
    }
    for(i = 0; i<n; i++){
        for(j = n; j>=i; j--){
            matrix[i][j]/=matrix[i][i];
        }
        for(k = i+1; k<n; k++){
            mx1 = matrix[k][i];
            for(j = 0; j<=n; j++){
                matrix[k][j]-=(mx1*matrix[i][j]);
            }
        }
    }
    for(i = n-1; i>=0; i--){
        _vector[i]=matrix[i][n];
        for(j = i+1; j<n; j++){
            _vector[i]-=(_vector[j]*matrix[i][j]);
        }
        delete matrix[i];
    }delete matrix;
}

double Delta(double a, double b, double eps, int m)
{
    return sqrt(integral(a, b, Sigma, m, -1, eps)/(b - a));
}

int main()
{
    double matrix[N][N + 1];
    int i;
    double delta;
    double x;
    double h;
    std::ofstream file(FILE_NAME);
    printf("Сalculation of polynomial\n");
    get_matrix(matrix, a, b, eps_i);
    printf("Searching polynomial degree where delta P(m) < %.0e\n", d_eps);
    for (i = 1; i < N; i++) {
        GaussianElimination(matrix, _vector, i);
        delta = Delta(a,b,eps_i, i + 1);
        if (delta < d_eps) {
            printf("m = %d\n", i + 1);
            printf("\nRoots of a generalized polynomial:\n");
            for (int j = 0; j < i + 1; j++) {
                printf("a%d = %e\n", j, _vector[j]);
            }
            printf("\nWriting XY coordinates to file \"%s\" for delta P(%d) = %.2e\n", FILE_NAME, i + 1, delta);
            printf("XY coordinates:\n");
            h = (double)(b-a)/numb;
            for (x = a; x <= b; x += h)
            {
                printf("%.20f, %.20f\n", x, get_Pm(x, N, _vector));
                file << std::setprecision(20) << x << ";" << get_Pm(x, N, _vector) << std::endl;
            }
            file.close();
            break;
        }
    }
    if (i == N) {
        printf("Something went wrong. Try to reduce d_eps or increase N.\n");
    }
    return 0;
}
