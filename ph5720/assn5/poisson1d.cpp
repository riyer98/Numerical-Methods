//Solving Poisson's equation by solving tridiagonal matrix
#include<iostream>
#include<cmath>
#include<fstream>
#include "time.h"
using namespace std;

//Function to solve for matrix A, an Nx3 matrix (see main fn for why it is Nx3)
void trisolve(double** A, double* b, int N)
{
	int i, j;
	//I use row operations to convert A (augmented with b) into a completely diagonal matrix, making off-diagonal elements 0.
	//Note: Each row only has 2 off-diagonal elements, and dividing by the diagonal element gives you the answer, we only require
	//approximately 3n~O(n) operations for an nxn tridiagonal matrix. In comparison, the standard Gauss elimination requires O(n^3) operations, much longer.
	for (i = N - 1; i > 1; i--)
	{
		b[i - 1] -= b[i] * A[i - 1][2] / A[i][1];
		for (j = 0; j < 2; j++) A[i - 1][j + 1] -= A[i][j] * A[i - 1][2] / A[i][1];

	}
	b[0] -= b[1] * A[0][1] / A[1][1];
	for (j = 0; j < 2; j++) A[0][j] -= A[1][j] * A[0][1] / A[1][1];

	b[1] -= b[0] * A[1][0] / A[0][0];
	A[1][0] = 0;
	for (i = 1; i < N - 1; i++)
	{
		b[i + 1] -= A[i + 1][0] * b[i] / A[i][1];
		A[i + 1][0] = 0;
	}
	
	b[0] /= A[0][0];
	for (i = 1; i < N; i++) b[i] /= A[i][1];
}

//I am using command line arguments for the main function to input the file name("file.dat") in which I store the data.
//e.g. To save data in a file named "poisson.dat", in terminal we need to write "./a.out poisson.dat" (see line 108)
int main(int argc, char** argv)
{
	double** A, * b, * x, * analytic, err, errmax, tri[3], h; int i, j, N;
	clock_t start, stop;
	
	cout << "enter no. of steps: ";
	cin >> N;

	cout << "enter the 3 tridiagonal matrix coefficients: ";
	for (i = 0; i < 3; i++) cin >> tri[i];

	start = clock();

	//Instead of defining A as an NxN matrix I have defined it as an Nx3 matrix (last and first row have only 2 elements) 
	//because at most 3 elements in each row are non-zero. By avoiding storing zeroes, I have greatly reduced the memory allocation, 
	//thereby increasing speed of the program.
	A = new double* [N];
	A[0] = new double[2];
	for (i = 1; i < N - 1; i++) A[i] = new double[3];
	A[N - 1] = new double[2];

	b = new double[N];
	x = new double[N];
	analytic = new double[N];


	A[0][0] = tri[1]; A[0][1] = tri[2];
	A[N - 1][0] = tri[0]; A[N - 1][1] = tri[1];

	for (i = 1; i < N - 1; i++)
	{
		for (j = 0; j < 3; j++) A[i][j] = tri[j];
	}
	
	h = 1. / (N + 1);

	for (i = 0; i < N; i++)
	{
		x[i] = (i + 1) * h;
		b[i] = h * h * (3 * x[i] + x[i] * x[i]) * exp(x[i]);
		analytic[i] = x[i] * (1 - x[i]) * exp(x[i]);
	}
	
	trisolve(A, b, N);

	//Finding the maximum error here(log) defined as log10 (|vi - ui|/|ui|), where ui is the analytic value.
	errmax = log10(abs((b[0] - analytic[0]) / analytic[0]));
	for (i = 1; i < N; i++)
	{
		err = log10(abs((b[i] - analytic[i]) / analytic[i]));
		if (err > errmax) errmax = err;
	}

	stop = clock();
	
	//Storing data in the dat file specified in the command line.
	if (argc > 1)
	{
		ofstream data;

		data.open(argv[1], ios::trunc);

		data << "#x value" << "	" << "u(x) numerical" << "	" << "u(x) analytical" << endl;
		data << 0 << "		" << 0 << "		" << 0 << endl;
		for (i = 0; i < N; i++) data << x[i] << "	" << b[i] << "		" << analytic[i] << endl;
		data << 1 << "		" << 0 << "		" << 0 << endl;

		data << endl << "#maximum error(log) = " << errmax << endl;

		data.close();

		cout << endl << "Values of x and u(x) (both numerical and analytical) and relative error have been stored in " << argv[1] << endl;
	}
	
	cout << "maximum error(log) = " << errmax << endl;
	cout << endl << "time taken = " << ((float)(stop - start)) / CLOCKS_PER_SEC << " s" << endl;
	
	return 0;
}
/*ERROR AS A FUNCTION OF h 
NOTE: Maximum error computed for each n value has been printed in the dat file. You can see the behaviour of the error in it.

Total error is a sum of approximation error and roundoff (precision) error. For larger h, approximation error is higher and roundoff error is less, 
because we are using less number of points. For smaller h, the approximation error reduces, but since we are summing over many terms, the roundoff error
(machine precision limit) cascades with every point for which we solve the equations. 
The overall behaviour is expected to decrease, reach a minimum, then increase again.*/

/*TIME COMPARISON WITH LU DECOMPOSITION (see the jpg file attached)
For a 1000x1000 matrix, the tridiagonal solver takes only ~1 ms (shows 0 s in photo as it is actually <1 ms) but the LU decomposition code 
(see pois_with_lud.cpp) takes ~2.5 to 3 s, almost 2500 times longer. This shows the efficiency of solving a tridiagonal matrix. Even for n=10^5 it is
able to solve the equations in about 0.03 s*/