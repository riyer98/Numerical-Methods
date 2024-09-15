//Program to multiply arbitrary MxN matrices

#include<iostream>
using namespace std;

int main()
{
	int M1, N1, M2, N2, i, j, k;
	cout << "enter dimensions of matrix 1(M1xN1) as M<space>N:";
	cin >> M1 >> N1;
	cout << "enter dimension of matrix 2(M2xN2) in the same way:";
	cin >> M2 >> N2;
	
	//To multiply 2 matrices, the number of columns of the first should equal the number of
	//rows in the second. If it is not, the program returns.
	if (N1 != M2)
	{
		cout << "Error, cannot multiply as M2 must be equal to N1" << endl;
		return 0;
	}

	float **A, **B;
	A = new float*[M1];
	for(i = 0;i < M1; i++) A[i] = new float[N1];

	B = new float*[M2];
	for(i = 0;i < M2; i++) B[i] = new float[N2];
	
	//Input the two matrices.
	cout << "Enter first matrix:" << endl;
	for (i = 0;i < M1;i++) 
	{
		for (j = 0;j < N1;j++) 
		{
			cin>>A[i][j];
		}
	}
	
	cout << endl;
	
	cout << "Enter second matrix:" << endl;
	for (i = 0;i < M2;i++) 
	{
		for (j = 0;j < N2;j++) 
		{
			cin >> B[i][j];
		}
	}
	
	cout << endl;

	cout << "Result is:" << endl;
	//Printing the result: we are summing A[i][j]*B[j][k] to get the answer C[i][k]
	//Instead of defining another matrix C i have directly printed the answer.
	for (i = 0;i < M1;i++)
	{
		for (k = 0;k < N2;k++) 
		{
			float sum = 0;
			
			for (j = 0;j < N1;j++) 
			{
				sum += A[i][j] * B[j][k];
			}
			
			cout << sum << " ";
		}
		
		cout << endl;
	}

	return 0;
}
