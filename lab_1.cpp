#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

const long double a = -7, b = 11, X =(a+b)/2, h = (b - a) / 10;

long double ret_R(long double u) {
	return 2*u/3;
}

long double u_next(long double x, long double u, int k) {
	return x*x*u/(2*k*(2*k-1));
}

int main() {
	long double Sum, SumU = 0, u_n = 1, u = 1, R, Y = cosh(X);
	int n = 0, k = 1;

	cout << setw(5) << "eps" << setw(5) << "n" << setw(20) <<\
			"Absolute error" << setw(20) << "Remainder term" << endl;

	for (long double eps = 1e-2; eps >= 1e-14; eps /= 1000) {
		do {
			SumU += u;
			u = u_next(X, u, k);
			k++;
		} while (u > eps);
		//SumU += u;
		R = ret_R(u);
		Sum = SumU + u;

		if (eps == 1e-8) n = k;

		cout << setw(5) << eps << " " << setw(5) << k-1 << " " <<\
				setw(20) << setprecision(10) << abs(Sum - Y) <<\
				" " << setw(20) << R << endl;
	}
	cout << endl;

	cout << setw(5) << "Xi" << setw(20) << "Absolute error" <<\
			setw(20) << "Remainder term" << endl;

	long double x_i;
	for (int i = 0; i <= 10; i++) {
		x_i = a + h*i;
		SumU = 0; u = 1;
		for (int j = 1; j < n; j++) {
			SumU += u;
			u = u_next(x_i, u, j);
		}
		//SumU += u;
		R = ret_R(u);
		Sum = SumU + u;
		cout << setw(5) << x_i << setw(20) <<\
				abs(Sum - cosh(x_i)) << setw(20) << R << endl;
	}

	system("pause");
	return 0;
}