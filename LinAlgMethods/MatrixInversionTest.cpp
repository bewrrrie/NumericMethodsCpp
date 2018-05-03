#include "LinAlgLib.cpp"
#include "Mygen.cpp"// <<---- I don't know what in this file but it works nice.
//                            Just use wrapper functions to generate matrices.
#include <iostream>
#include <cmath>

using namespace linalg;
using namespace std;

#define N 100
#define ALPHA 1.e-3
#define BETA 1.


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
	double alpha = 1;
	double beta = 100;

	//printing data

	return 0;
}