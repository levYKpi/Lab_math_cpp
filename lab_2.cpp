#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <iomanip>

using namespace std;
//x^2 - sin^2(x) + cos(x) - 1.5 = 0
//x = +-sqrt(sin^2(x) - cos(x) + 1.5)

double f(double x) {
	return x * x - sin(x)*sin(x) + cos(x) - 1.5;
}

double f1(double x) {
	return 2 * x - 2 * sin(x)*cos(x) - sin(x);
}

double phy(double x) {
	if(f1(x) >= 0)
		return sqrt(sin(x)*sin(x) - cos(x) + 1.5);
	else
		return -sqrt(sin(x)*sin(x) - cos(x) + 1.5);
}

double dphy(double x) {
	if (f1(x) >= 0)
		return (2 * sin(x)*cos(x) + sin(x)) / (2 * sqrt(sin(x)*sin(x) - cos(x) + 1.5));
	else
		return -(2 * sin(x)*cos(x) + sin(x)) / (2 * sqrt(sin(x)*sin(x) - cos(x) + 1.5));
}

double f2(double x) {
	return 2 - 2 * cos(2 * x) - cos(x);
}

double find(double a, double b, double eps, bool p) {
	double q;
	double x0 = a;
	double xk = phy(x0);
	int i = 1;
	if (dphy(a) > dphy(b)) q = dphy(a); else q = dphy(b);
	if (q >= 1) { cout << "error iter find() q \n"; return xk; }
	while (abs(xk - x0) > (1 - q) / q * eps) {
		i++;
		x0 = xk;
		xk = phy(x0);
	}
	cout << setprecision(14) << setw(14);
	if (!p)
		cout << eps << "  " << xk << "  " << abs(xk - x0) << endl;
	else
		cout << setw(7) << eps << "  " << i;
	return xk;
}

double Nfind(double a, double b, double eps, bool p) {
	double ml;
	double xk;
	int i = 1;
	if (abs(f1(a)) < abs(f1(b))) ml = abs(f1(a)); else ml = abs(f1(b));
	if (f2(a)*f(a) > 0) xk = a; else xk = b;
	if (ml <= 0) { cout << "error Nfind() ml \n"; return xk; }
	while (abs(f(xk) / ml) > eps) {
		i++;
		xk -= f(xk) / f1(xk);
	}
	cout << setprecision(14) << setw(14);
	if (!p)
		cout << eps << "  " << xk << "  " << abs(f(xk) / ml) << endl;
	else
		cout << setw(7) << "  " << i << endl;
	return xk;
}

int main() {
	double eps = 0.01;
	cout << "First x\nIter.\n          eps[i] X                accuracy" << endl;
	for (int i = 0; i < 5; i++, eps *= 1e-3) {
		find(1.4, 1.6, eps, false);
	}cout << endl;
	eps = 0.01;
	cout << "Newton\n         eps[i] X                accuracy" << endl;
	for (int i = 0; i < 5; i++, eps *= 1e-3) {
		Nfind(1.4, 1.6, eps, false);
	}cout << endl;
	eps = 0.01;
	cout << " eps[i]  I(i)   N(i)" << endl;
	for (int i = 0; i < 5; i++, eps *= 1e-3) {
		find(1.4, 1.6, eps, true);
		Nfind(1.4, 1.6, eps, true);
	}cout << endl;
	eps = 0.01;
	cout << "Second x\nIter.\n          eps[i] X                accuracy" << endl;
	for (int i = 0; i < 5; i++, eps *= 1e-3) {
		find(-1.6, -1.4, eps, false);
	}cout << endl;
	eps = 0.01;
	cout << "Newton\n         eps[i] X                accuracy" << endl;
	for (int i = 0; i < 5; i++, eps *= 1e-3) {
		Nfind(-1.6, -1.4, eps, false);
	}cout << endl;
	eps = 0.01;
	cout << " eps[i]  I(i)   N(i)" << endl;
	for (int i = 0; i < 5; i++, eps *= 1e-3) {
		find(-1.6, -1.4, eps, true);
		Nfind(-1.6, -1.4, eps, true);
	}cout << endl;
	system("pause");
	return 0;
}
