//A program to calculate the inverse of an NxN matrix (N<=M)

#include<iostream>
#include<cmath>
using namespace std;

//Determinant function. For explanation, kindly see my matdet.cpp program.
float det(float **A, int N)
{
	if (N == 1) return A[0][0];

	float** B, res = 0; int i, j, k;

	B = new float* [N-1];
	for (i = 0;i < N-1;i++) B[i] = new float[N-1];

	for (i = 0;i < N;i++)
	{
		for (j = 1;j < N;j++)
		{
			for (k = 0;k < N;k++)
			{
				if (k == i) continue;
				if (k < i) B[j - 1][k] = A[j][k];
				else B[j - 1][k - 1] = A[j][k];
			}
		}
		res += pow(-1, i) * A[0][i] * det(B, N - 1);
	}
	return res;
}


int main()
{
	float **A, **B, **cof; int i, j, k, l, N;

	cout << "enter dimension of matrix: ";
	cin >> N;

	A = new float* [N];
	for (i = 0;i < N;i++) A[i] = new float[N];

	B = new float* [N - 1];
	for (i = 0;i < N - 1;i++) B[i] = new float[N - 1];

	cof = new float* [N];
	for (i = 0;i < N;i++) cof[i] = new float[N];
	
	cout << "enter matrix" << endl;
	for (i = 0;i < N;i++)
	{
		for (j = 0;j < N;j++) cin >> A[i][j];
	}
	
	if (N == 1) 
	{
		cout << 1 / A[0][0] << endl; 
		return 0;
	}

	//The main part of the program, that calculates the cofactor matrix
	for (i = 0;i < N;i++)
	{
		for (j = 0;j < N;j++)
		{
			for (k = 0;k < N;k++)
			{
				if (k == i) continue;
				for (l = 0;l < N;l++)
				{
					//For every A[k][l], elements in the same row as A[k][l] are skipped using the
					//continue command. Others are stored in the temporary matrix B (the minor)
					if (l == j) continue;
					if (k < i && l < j) B[k][l] = A[k][l];
					if (k < i && l > j) B[k][l - 1] = A[k][l];
					if (k > i && l < j) B[k - 1][l] = A[k][l];
					if (k > i && l > j) B[k - 1][l - 1] = A[k][l];
				}
			}
			//Here the cofactor matrix element cof[i][j] is defined as the determinant of B (with appropriate sign)
			cof[i][j] = pow(-1, (i + j)) * det(B,N-1);
		}
	}

	//Prints the adjoint, which is the transpose of the cofactor, hence it prints cof[j][i]
	cout << endl << "adjoint is" << endl;
	for (i = 0;i < N;i++)
	{
		for (j = 0;j < N;j++)
		{
			//Note: this if condition has been put to overcome the floating point error. Sometimes, due to error in
			//calculation, the computer prints something like 2.23543 x 10^-38 instead of 0. Hence it is asked to
			//print 0 if the number is below the precision (same logic as taught in class)
			if (1 + cof[j][i] == 1) cout << 0 << " ";
			else cout << cof[j][i]  << " ";
		}
		cout << endl;
	}
	
	//Prints the inverse which is adj(A)/det(A). 
	cout << endl << "inverse is:" << endl;
	for (i = 0;i < N;i++)
	{
		for (j = 0;j < N;j++)
		{
			//if statement is present for the same reason as the case of printing adjoint.
			if (1 + cof[j][i] == 1) cout << 0 << " ";
			else cout << cof[j][i]/det(A,N) << " ";
		}
		cout << endl;
	}

	for (i = 0;i < N;i++) delete[] A[i];
	delete[] A;

	for (i = 0;i < N-1;i++) delete[] B[i];
	delete[] B;
	
	for (i = 0;i < N;i++) delete[] cof[i];
	delete[] cof;

	return 0;
}



	
