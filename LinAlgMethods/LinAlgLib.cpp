#include <cmath>
#include <iostream>

using namespace std;

namespace linalg
{
	double** copy_matrix(double** A, size_t n)
	{
		double** A_copy = new double*[n];

		for(int i = 0; i < n; ++i)
		{
			A_copy[i] = new double[n];

			for (int j = 0; j < n; ++j)
				A_copy[i][j] = A[i][j];
		}

		return A_copy;
	}


	double** create_zero_matrix(size_t n)
	{
		double** A = new double*[n];

		for(int i = 0; i < n; ++i)
			A[i] = new double[n];

		return A;
	}


	double** create_unit_matrix(size_t n)
	{
		double** E = create_zero_matrix(n);

		for (int i = 0; i < n; ++i)
		{
			E[i][i] = 1;
		}

		return E;
	}


	double l_inf_mnorm(double** A, size_t n)
	{
		double norm = 0;

		for (int i = 0; i < n; ++i)
		{
			double current = 0;

			for (int j = 0; j < n; ++j)
				current += fabs(A[i][j]);

			if (current > norm)
				norm = current;
		}

		return norm;
	}


	double** subtract_m(double** A, double** B, size_t n)
	{
		double** C = create_zero_matrix(n);

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				C[i][j] = A[i][j] - B[i][j];

		return C;
	}


	double** add_m(double** A, double** B, size_t n)
	{
		double** C = create_zero_matrix(n);

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				C[i][j] = A[i][j] + B[i][j];

		return C;
	}


	double** mult_m(double** A, double** B, size_t n)
	{
		double** C = create_zero_matrix(n);

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				for (int k = 0; k < n; ++k)
					C[i][j] += A[i][k] * B[k][j];

		return C;
	}


	void print_matrix(double** A, size_t n)
	{
		for (int i = 0; i < n; ++i)
		{
			cout << A[i][0];
			for (int j = 1; j < n; ++j)
				cout << " " << A[i][j];
			cout << endl;
		}
	}


	void transpose_columns(double** A, size_t n, size_t i, size_t j)
	{
		for (int k = 0; k < n; ++k)
		{
			double tmp = A[k][i];
			A[k][i] = A[k][j];
			A[k][j] = tmp;
		}
	}

	void transpose_lines(double** A, size_t n, size_t i, size_t j)
	{
		double* pointer = A[i];
		A[i] = A[j];
		A[j] = pointer;
	}


	double** invert_lower_triangular(double** A, size_t n)
	{
		double** A_inv = create_zero_matrix(n);

		for (int i = 0; i < n; ++i)
		{
			A_inv[i][i] = 1 / A[i][i];
		}

		for (int i = 1; i < n; ++i)
		{
			for (int j = i - 1; j > -1; --j)
			{
				double summ = 0;

				for (int k = 0; k < i; ++k)
				{
					summ -= A_inv[k][j] * A[i][k];
				}

				A_inv[i][j] = summ / A[i][i] ;
			}
		}

		return A_inv;
	}

	double** invert_upper_triangular(double** A, size_t n)
	{
		double** A_inv = create_zero_matrix(n);

		for (int i = 0; i < n; ++i)
		{
			A_inv[i][i] = 1 / A[i][i];
		}

		for (int i = n - 2; i > -1; --i)
		{
			for (int j = i + 1; j < n; ++j)
			{
				double summ = 0;

				for (int k = i + 1; k < n; ++k)
				{
					summ -= A_inv[k][j] * A[i][k];
				}

				A_inv[i][j] = summ / A[i][i] ;
			}
		}

		return A_inv;
	}


	double** gauss_invert(double** A, size_t n)
	{
		double** E = create_unit_matrix(n);
		int* transposes = new int[n];

		for (int i = 0; i < n; ++i)
		{
			double k = A[i][i];
			int tmp_j = i;

			for (int j = 0; j < n; ++j)
				if (fabs(A[i][j]) > fabs(k))
				{
					k = A[i][j];
					tmp_j = j;
				}

			if (i != tmp_j)
			{
				transpose_columns(A, n, i, tmp_j);
			}

			transposes[i] = tmp_j;

			for (int j = 0; j < n; ++j)
			{
				A[i][j] /= k;
				E[i][j] /= k;
			}


			for (int l = i + 1; l < n; ++l)
			{
				k = A[l][i];

				for (int j = 0; j < n; ++j)
				{
					A[l][j] -= k * A[i][j];
					E[l][j] -= k * E[i][j];
				}
			}
		}

		double** A_triangular_inv = invert_upper_triangular(A, n);

		for (int i = n - 1; i > -1; --i)
			transpose_lines(A_triangular_inv, n, i, transposes[i]);

		return mult_m(A_triangular_inv, E, n);
	}
}