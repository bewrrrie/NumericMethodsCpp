#include "LinAlgLib.cpp"
#include "Mygen.cpp"// <<---- I don't know what in this file but it works nice.
//                            Just use wrapper functions to generate matrices.
#include <iostream>
#include <cmath>

using namespace linalg;
using namespace std;

#define N 100


//!! Matrix generation function wrappers !!//

void gen_symmetric(double** A, double** A_inv, size_t n, double alpha, double beta)
{
	mygen (A, A_inv, n, alpha, beta, 1, 2, 0, 1); // symmetric
}

void gen_simple_struct(double** A, double** A_inv, size_t n, double alpha, double beta)
{
	mygen (A, A_inv, n, alpha, beta, 1, 2, 1, 1); // simple structure
}

void gen_Jordan(double** A, double** A_inv, size_t n, double alpha, double beta)
{
	mygen (A, A_inv, n, alpha, beta, 0, 0, 2, 1); // Jordan form
}




int main()
{
	double alpha, beta;
	double** A = new double*[N];
	double** A_copy = new double*[N];
	double** A_inv = new double*[N];
	double** error = new double*[N];
	double** residual = new double*[N];
	double** E = create_unit_matrix(N);
	double** A_inv_star;

	for (int i = 0; i < N; ++i)
	{
		A[i] = new double[N];
		A_inv[i] = new double[N];
		error[i] = new double[N];
		residual[i] = new double[N];
	}

	cout << "--------------------------------" << endl;
	cout << "Enter ALPHA value:";
	cin >> alpha;
	cout << "Enter BETA value:";
	cin >> beta;
	cout << endl;

	cout << "ALPHA = " << alpha << ", BETA = " << beta << endl;
	gen_Jordan(A, A_inv, N, alpha, beta);

	A_copy = copy_matrix(A, N);
	A_inv_star = gauss_invert(A, N);
	error = subtract_m(A_inv_star, A_inv, N);
	residual = subtract_m( mult_m(A_copy, A_inv_star, N), E, N );

	cout << endl << "abs.error = " << l_inf_mnorm(error, N) << endl;
	cout << "abs.residual = " << l_inf_mnorm(residual, N) << endl;
	cout << "rel.error = " << l_inf_mnorm(error, N) / l_inf_mnorm(A_inv, N) << endl;
	cout << "--------------------------------" << endl;

	return 0;
}