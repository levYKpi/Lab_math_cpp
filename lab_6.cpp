#include <stdio.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

std::ofstream _file;
const double a = 0, b = 0.5;

double yy(double x){
    return 4*exp(x)+2*exp(3*x);
}

double funcQ(double x, double y, double t){
    return 4*t-3*y;
}

double funcT(double x, double y, double t){
    return t;
}

void runge_kutt(double x, double y, double t,\
                double& _x, double& _y, double& _t,\
                double h){
    double k1,k2,k3,k4,q1,q2,q3,q4;
    k1 = funcT(x,y,t);
    q1 = funcQ(x,y,t);
    k2 = funcT(x+h/2,y+h*k1/2,t+h*q1/2);
    q2 = funcQ(x+h/2,y+h*k1/2,t+h*q1/2);
    k3 = funcT(x+h/2,y+h*k2/2,t+h*q2/2);
    q3 = funcQ(x+h/2,y+h*k2/2,t+h*q2/2);
    k4 = funcT(x+h,y+h*k3,t+h*q3);
    q4 = funcQ(x+h,y+h*k3,t+h*q3);
    _x = x + h;
    _y = y + h*(k1+2*k2+2*k3+k4)/6;
    _t = t + h*(q1+2*q2+2*q3+q4)/6;
}

double RK4(int N, double x, double y, double t, double eps, void (*output)(double, double)){
    int i;
    double h, deps;
    while(1){
        h = (b-a)/N;
        output(x,y);
        for(i = 0; i<N; i++){
            runge_kutt(x,y,t,x,y,t,h);
            output(x,y);
        }
        if((deps = fabs(y - yy(b)))<eps){
            break;
        }else{
            x = a; y = 6; t = 10;
            N++;
        }
    }
    printf("method RG4\neps: %.20f\nN: %d\n",deps,N);
    return deps;
}

double *double_count(double eps, double x0, double y0, double t0, double h, int N){
    double eta, y, y2, x, xi, xi_1, t, T[N]; int _count = 0, n, i, k;
    double *Y;
    Y = new double[N];
    Y[0] = y0;
    xi = a + h;
    xi_1 = a;
    T[0] = t0;
    for(i = 1; i<N; i++){
        n = 1 + (b-a)/sqrt(sqrt(eps));
        eta = (xi - xi_1)/n;
        x = xi_1;
        y = Y[i-1];
        t = T[i-1];
        for(k = 0; k<n; k++){
            runge_kutt(x,y,t,x,y,t,eta);
        }
        n*=2;
        eta=(xi-xi_1)/n;
        y2 = Y[i-1];
        t = T[i-1];
        for(k = 0; k<n; k++){
            runge_kutt(x,y2,t,x,y2,t,eta);
        }
        while(fabs(y2 - y) > 15*eps){
            y = y2;
            n*=2;
            eta = (xi - xi_1)/n;
            x = xi_1;
            y2 = Y[i-1];
            t = T[i-1];
            for(k = 0; k<n; k++){
                runge_kutt(x,y2,t,x,y2,t,eta);
            }
        }
        Y[i] = y2;
        T[i] = t;
        xi_1 = xi;
        xi += h;
        _count+=n;
    }
    printf("method double_count\neps: %.20f\nN: %d\n",fabs(y2 - y)/15.0,_count);
    return Y;
}

void print_table(double x, double y){
    printf("%.10f;%.20f\n",x,y);
    if(_file.is_open()){
        _file << std::setprecision(5) << x << ";" << std::setprecision(20) << y << std::endl;
    }
}


int main(int argc, char **argv)
{
    double eps;
    double *Y;
    _file.open("table.csv");
    eps = RK4(10,a,6,10,0.01,print_table);
    Y = double_count(eps,a,6,10,(b-a),2);
    printf("\n");
    _file.close();
    delete Y;
    return 0;
}
