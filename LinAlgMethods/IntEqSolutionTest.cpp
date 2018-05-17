#include "LinAlgLib.cpp"
#include <iostream>
#include <cmath>

using namespace linalg;
using namespace std;

#define N 3

int main()
{
	double* b = new double[N]; b[0] = 1; b[1] = 5; b[2] = 9;
	double* x = new double[N]; x[0] = 1; x[1] = 1; x[2] = 1;

	double** A = new double*[N];
	for (int i = 0; i < N; ++i)
		A[i] = new double[N];

	A[0][0] = 3; A[0][1] = 4; A[0][2] = 0;
	A[1][0] = 4; A[1][1] = -3; A[1][2] = 0;
	A[2][0] = 0; A[2][1] = 0; A[2][2] = 5;

	conjugate_gradient_method(A, b, x, 0.001, N);/*
	cout << "(" << x[0];
	for (int i = 0; i < N; ++i)
		cout << ", " << x[i];
	cout << ")" << endl;
*/
	return 0;
}