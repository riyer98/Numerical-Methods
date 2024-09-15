//SOLVING POISSON'S EQUATION BY USING LU DECOMPOSITION OF THE TRIDIAGONAL MATRIX.

//This is an extra code to demonstrate that solving a very large matrix by LU decomposition is not feasible.
//I ran this code to solve a 1000x1000 tridiagonal matrix by LU decomposition.
#include<iostream>
#include<cmath>
#include<fstream>
#include "time.h"
using namespace std;

void decomp(int N, double** A, double** L, double** U, double* b, double* x)
{
	int i, j, k; double* y, sum;
	y = new double[N];

	for (i = 0; i < N; i++) L[i][i] = 1;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < i; j++)
		{
			U[i][j] = 0;
			L[j][i] = 0;
		}
	}

	for (i = 0; i < N; i++)
	{
		for (j = i; j < N; j++)
		{
			sum = 0;
			for (k = 0; k < i; k++) sum += L[i][k] * U[k][j];
			U[i][j] = A[i][j] - sum;
			sum = 0;
			for (k = 0; k < i; k++)  sum += L[j][k] * U[k][i];
			L[j][i] = (A[j][i] - sum) / U[i][i];
		}
	}

	for (i = 0; i < N; i++)
	{
		sum = 0;
		for (j = 0; j < i; j++)  sum += L[i][j] * y[j];
		y[i] = b[i] - sum;
	}

	for (i = N - 1;i >= 0;i--)
	{
		sum = 0;
		for (j = i + 1;j < N; j++) sum += U[i][j] * x[j];
		x[i] = (y[i] - sum) / U[i][i];
	}
}

int main(int argc, char** argv)
{
	double** A, ** L, ** U, * b, * x, * analytic, tri[3], * ux, h, err, errmax; int i, j, N;
	clock_t start, stop;
	cout << "enter no. of steps: ";
	cin >> N;

	cout << "enter the 3 tridiagonal matrix coefficients: ";
	for (i = 0; i < 3; i++) cin >> tri[i];

	start = clock();

	A = new double* [N];
	for (i = 0; i < N; i++) A[i] = new double[N];

	L = new double* [N];
	for (i = 0; i < N; i++) L[i] = new double[N];

	U = new double* [N];
	for (i = 0; i < N; i++) U[i] = new double[N];

	b = new double[N];
	ux = new double[N];
	x = new double[N];
	analytic = new double[N];

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++) A[i][j] = 0;
	}

	A[0][0] = tri[1]; A[0][1] = tri[2];
	A[N - 1][N - 2] = tri[0]; A[N - 1][N - 1] = tri[1];

	for (i = 1; i < N - 1; i++)
	{
		for (j = -1; j < 2; j++) A[i][i + j] = tri[j + 1];
	}

	h = 1. / (N + 1); 

	for (i = 0; i < N; i++)
	{
		x[i] = (i + 1) * h;
		b[i] = h * h * (3 * x[i] + x[i] * x[i]) * exp(x[i]);
		analytic[i] = x[i] * (1 - x[i]) * exp(x[i]);
	}
	
	decomp(N, A, L, U, b, ux);
	
	errmax = log10(abs((ux[0] - analytic[0]) / analytic[0]));
	for (i = 1; i < N; i++);
	{
		err = log10(abs((ux[i] - analytic[i]) / analytic[i]));
		if (err > errmax) errmax = err;
	}
	stop = clock();

	if (argc > 1)
	{
		ofstream data;

		data.open(argv[1], ios::trunc);

		data << "#x value" << "	" << "u(x) numerical" << "	" <<"u(x) analytical" << endl;
		data << 0 << "		" << 0 << "		" << 0 << endl;
		for (i = 0; i < N; i++) data << x[i] << "	" << ux[i] << "		" << analytic[i] << endl;
		data << 1 << "		" << 0 << "		" << 0 << endl;

		data << endl << "Maximum error(log) = " << errmax;
		data.close();

		cout << endl << "Values of x and u(x) (both numerical and analytical) have been stored in " << argv[1] << endl;
	}
	cout << "Maximum error(log) = " << errmax << endl;
	cout << endl << "time taken = " << ((float)(stop - start)) / CLOCKS_PER_SEC << " s" << endl;
}
/*TIME COMPARISON WITH TRIDIAGONAL SOLVER
This code took ~2.5 to 3 s to give result for a 1000x1000 matrix, as compared to the tridiagonal solver, which did it in ~1 ms.
We can see that for large matrices of order >10^5 or so, this method is not feasible.*/