#include <iostream>
#include <iomanip>
#include <math.h>
#define M_PI 3.14159265358979323846

using namespace std;

double errors[4] = { 0 };

double funk(double x) {
	return 1 + 0.01*exp(x)*(sin(x) + cos(x));
}

double Funk(double x) {
	return 0.01*exp(x)*sin(x) + x;
}
/*
double funk2(double x) {
	return 0.02*exp(x)*cos(x);
}*/

double Integral(double B, double A) {
	return Funk(B) - Funk(A);
}
/*
int nIter(double b, double a, double h) {
	return (int)((b-a)/h) + 1;
}

double hIter(double b, double a, double e) {
	return sqrt(12 * e / (b - a) / funk2(9 * M_PI / 4));
}*/

double tIntegral(double a, double h, int n) {
	double Int = 0.0;
	for (int i = 1; i<n; i++) {
		Int += h * funk(a + h * i);
	}
	Int += (funk(a) + funk(a + n * h))*h / 2;
	return Int;
}

void table1(double a, double b) {
	double h, int2, integral = Integral(b, a);
	int n, i = 0;
	cout << setw(14) << "eps" << setw(14) << "step" << setw(14) << \
		"Int" << setw(14) << "Int2" << setw(20) << "error\n";
	for (double e = 1e-2; e >= 1e-8; e *= 1e-2) {
		h = sqrt(e);//hIter(b, a, e);
		n = (int)((b - a) / h);//nIter(b, a, h);
		h = (b - a) / n;
		int2 = tIntegral(a, h, n);
		errors[i] = abs(integral - int2);
		cout << setw(14) << e << setw(14) << n << setw(18) << \
			setprecision(12) << integral << setw(18) << \
			int2 << setw(20) << errors[i] << endl;
		i++;
	}
}

void RefinedCalculation(double a, double b, double eps) {
	double integral = Integral(b, a);
	double Int1 = 0.0;
	double Int2 = 0.0;
	double h = sqrt(eps);
	int n = (int)((b - a) / h);
	h = (b - a) / n;
	Int1 = tIntegral(a, h, n);
	n *= 2;
	h /= 2;
	Int2 = tIntegral(a, h, n);
	while (abs(Int1 - Int2) > 3 * eps) {
		Int1 = Int2;
		n *= 2; h /= 2;
		Int2 = tIntegral(a, h, n);
	}
	cout << setw(20) << eps << setw(14) << n << setw(20) << abs(integral - Int2) << endl;
}

void table2(double a, double b) {
	cout << setw(14) << "eps" << setw(14) << "step" << setw(14) << "error\n";
	for (int i = 0; i<4; i++) {
		RefinedCalculation(a, b, errors[i]);
	}
}

int main() {
	table1(2, 11);
	cout << endl;
	table2(2, 11);
	system("pause");
	return 0;
}
