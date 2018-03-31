#include "LinAlgLib.cpp"
#include <iostream>
#include <cmath>

using namespace linalg;
using namespace std;

#define N 100


int main()
{
	/*
	!!! N = 3 !!!
	double** A = create_zero_matrix(N);
	A[0][0] = 3;	A[0][1] = 2;	A[0][2] = 1;
	A[1][0] = 2;	A[1][1] = 2;	A[1][2] = 1;
	A[2][0] = 1;	A[2][1] = 1;	A[2][2] = 1;

	double** A_inv = create_zero_matrix(N);
	A_inv[0][0] = 1;	A_inv[0][1] = -1;	A_inv[0][2] = 0;
	A_inv[1][0] = -1;	A_inv[1][1] = 2;	A_inv[1][2] = -1;
	A_inv[2][0] = 0;	A_inv[2][1] = -1;	A_inv[2][2] = 2;

	double** A_star = gauss_invert(A, N);
	print_matrix(A_star, N);
	*/

	double** A = create_zero_matrix(N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			A[i][j] = 2 + cos(i + j + 3);
		}
	}

	double** A_star = gauss_invert(A, N);
	double** E_star = mult_m(A, A_star, N);
	double** E = create_unit_matrix(N);

/*
	double error = l_inf_mnorm(subtract_m(E_star, E, N), N);
	double relative_error = error / l_inf_mnorm(A, N);

	print_matrix(A, N);
	cout << "------------------------" << endl;
	print_matrix(A_star, N);
	cout << "------------------------" << endl;
	print_matrix(E_star, N);
	cout << "------------------------" << endl;
	cout << "error = " << error << endl;
	cout << "relative error = " << relative_error << endl;
	*/

	return 0;
}
