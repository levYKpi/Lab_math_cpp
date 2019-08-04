#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

void printArr(vector<vector <double> > *arr, vector<double> *b, int n) {
	for (int i = 0; i<n; i++) {
		for (int j = 0; j<n; j++) {
			cout << setprecision(3);
			cout << setw(6);
			cout << (*arr)[i][j];
			if (j == n - 1) cout << "  =  " << (*b)[i] << endl;
			else cout << "  +  ";
		}
	}
}

void JGArr(vector<vector <double> > *arr, vector<double> *b, int n) {
	for (int k = 0; k<n; k++) {
		if ((*arr)[k][k] == 0) break;
		double akk = (*arr)[k][k];
		for (int j = 0; j<n; j++) {
			(*arr)[k][j] /= akk;
		}(*b)[k] /= akk;
		for (int i = 0; i<n; i++) {
			if (i == k) continue;
			double aik = (*arr)[i][k];
			for (int j = 0; j<n; j++) {
				(*arr)[i][j] -= (*arr)[k][j] * aik;
			}(*b)[i] -= (*b)[k] * aik;
		}
	}
}

double mNormB(vector<double> *a, vector<double> *b, int n) {
	//vector<double> X(n, 0);
	double max = 0, x;
	for (int i = 0; i<n; i++) {
		x = abs((*a)[i] - (*b)[i]);
		/*X[i] = (*b)[i] - (*a)[i];
		if (X[i]<0) X[i] *= -1;
		if (X[i]>max) max = X[i];*/
		if (max < x) max = x;
	}
	return max;
}

double mNormA(vector<vector<double> > *arr, int n) {
	double sumb = 0;
	double maxa = 0;
	for (int i = 0; i<n; i++) {
		sumb = 0;
		for (int j = 0; j<n; j++) {
			sumb += abs((*arr)[i][j]);
		}if (sumb > maxa) maxa = sumb;
	}
	return maxa;
}

void Saidel(vector<vector <double> > *arr, vector<double> *b, int n, double eps) {
	vector<double> beta(n, 0), X(n, 0), X0(n, 0);
	vector<vector <double> > alpha(n, vector<double>(n, 0));
	double q = 0, norm = 0;
	int k = 1;
	for (int i = 0; i<n; i++) {
		beta[i] = (*b)[i] / (*arr)[i][i];
		for (int j = 0; j<n; j++) {
			if (i == j) continue;
			alpha[i][j] = -(*arr)[i][j] / (*arr)[i][i];
		}
	}
	cout << endl;
	X0 = beta;
	q = mNormA(&alpha, n);
	for (int i = 0; i<n; i++) {
		double Sum1 = 0, Sum2 = 0;
		for (int j = 0; j<i; j++) {
			if (i == j) continue;
			Sum1 += alpha[i][j] * X0[j];
		}for (int j = i; j<n; j++) {
			Sum2 += alpha[i][j] * (*b)[j];
		}
		X[i] = beta[i] + Sum1 + Sum2;
	}X0 = X;
	for (k = 2; k<100; k++) {
		for (int i = 0; i<n; i++) {
			double Sum1 = 0, Sum2 = 0;
			for (int j = 0; j<i; j++) {
				if (i == j) continue;
				Sum1 += alpha[i][j] * X[j];
			}for (int j = i; j<n; j++) {
				Sum2 += alpha[i][j] * X0[j];
			}
			X[i] = beta[i] + Sum1 + Sum2;
		}
		norm = mNormB(&X, &X0, n);
		if (norm <= eps * (1 - q) / q) break;
		X0 = X;
	}
	cout << "\nK: " << k << endl;
	cout << "Seidel X out:\n";
	for (int i = 0; i<n; i++) {
		cout << " X" << i << " :" << X[i];
	}cout << endl;
	if (q >= 1) cout << "not convergent\n";
}

int main() {
	double ARR[4][4] = { 2,1,2,9,16,48,15,16,1,4,12,18,10,7,36,18 };
	double B[] = { 13,46,42,90 };
	vector<vector <double> > arr, arr_r;
	vector<double> b, b_r;
	for (int i = 0; i<4; i++) {
		arr.push_back(vector<double>(ARR[i], ARR[i] + 4));
	}
	b = vector<double>(B, B + 4);
	printArr(&arr, &b, 4); cout << endl;

	arr_r = arr; b_r = b;
	JGArr(&arr_r, &b_r, 4);
	printArr(&arr_r, &b_r, 4); cout << endl;

	for (int j = 0; j<4; j++) {
		arr_r[0][j] = (arr[0][j] * 2 - arr[2][j])*4.5 + arr[3][j];
	}b_r[0] = (b[0] * 2 - b[2])*4.5 + b[3];
	arr_r[1] = arr[1]; b_r[1] = b[1];
	arr_r[2] = arr[3]; b_r[2] = b[3];
	arr_r[3] = arr[2]; b_r[3] = b[2];
	cout << "\nmodified matrix\n";
	printArr(&arr_r, &b_r, 4);

	Saidel(&arr_r, &b_r, 4, 1e-3);
	system("pause");
	return 0;
}
