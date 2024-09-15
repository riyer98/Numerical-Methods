//LU DECOMPOSITION CODE
#include<iostream>
using namespace std;

//Function that does the decomposition.
void decomp(int N, float **A, float **L, float **U, float *b, float *x)
{
	int i, j, k; float *y, sum;
	y = new float[N];

	for (i = 0; i < N; i++) L[i][i] = 1; 

	//Setting the upper elements of L and lower elements of U equal to 0.
	for(i = 0; i < N; i++)
	{
		for (j = 0; j < i; j++) 
		{
			U[i][j] = 0;
			L[j][i] = 0;
		}
	}

	for(i = 0; i < N; i++)
	{
		for(j = i; j < N; j++)
		{
			sum = 0;
			//U[i][j] = A[i][j] - Sum(k=0 to i-1)L[i][k]*U[k][j]
			for(k = 0; k < i; k++) sum+= L[i][k] * U[k][j];
			U[i][j] = A[i][j] - sum;
			sum = 0;
			//L[j][i] = 1/U[i][i] (A[j][i] - Sum(k=0 to i-1)L[j][k]*U[k][i])
			for (k = 0; k < i; k++)  sum += L[j][k] * U[k][i];
			L[j][i] = (A[j][i] - sum) / U[i][i];	
		}
	}
	
	//Solving Ly=b by forward substitution
	for (i = 0; i < N; i++) 
	{
		sum = 0;
		for (j = 0; j < i; j++)  sum += L[i][j] * y[j];
		y[i] = b[i] - sum;
	}
	
	//Solvinng Ux=y by backward substituion.
	for(i = N - 1;i >= 0;i--)
	{
		sum = 0;
		for(j = i + 1;j < N; j++) sum+= U[i][j] * x[j];
		x[i] = (y[i] - sum) / U[i][i];
	}
}
/*COMMENTS ON ORDER OF OPERATIONS
We are solving for nxn elements of the matrix A=n^2 elements.But solving for each element requires computing a sum of 
O(n) elements (see comments on line 28 and 32). Therefore, similar to Gauss elimination, LU decomposition requires O(n^3) operations.*/

int main()
{
	float **A, **L, **U, *b, *x; int N, i, j, temp;

	cout<<"Enter no of variables - ";
	cin >> N;
	
	A = new float*[N];
	for (i = 0; i < N; i++) A[i] = new float[N];

	L = new float*[N];
	for (i = 0; i < N; i++) L[i] = new float[N];

	U = new float*[N];
	for (i = 0; i < N; i++) U[i] = new float[N];

	b = new float[N];
	x = new float[N];

	cout<<endl<<"Enter coefficient matrix"<<endl;
	for (i = 0; i < N; i++) 
	{
		for (j = 0; j < N; j++) cin>>A[i][j];
	}
	
	cout << endl << "Enter rhs" << endl;
	for (i = 0;i < N;i++) cin >> b[i];

	//Pivoting done here only if the first element of A is close to 0. This is because U[0][0]=A[0][0] and diagonal elements of U should be non zero.
	if (1 + A[0][0] == 1) 
	{
		for (i = 1; i < N; i++) 
		{
			if (A[i][0] == 0) continue;
			for (j = 0; j < N; j++) 
			{
				temp = A[0][j];
				A[0][j] = A[i][j];
				A[i][j] = temp;
			}
			temp = b[0];
			b[0] = b[i];
			b[i] = temp;
			break;
		}
		
		cout<<endl<<"Matrix has been pivoted. New matrix:"<<endl;
		for (i = 0;i < N;i++) 
		{
			for (j = 0;j < N;j++) cout << A[i][j] << "	";
			cout << endl;
		}
	}				
	
	decomp(N, A, L, U, b, x);
	
	cout << endl << "L matrix" << endl;
	for (i = 0; i < N; i++) 
		{
			for(j = 0; j < N; j++) cout<<L[i][j]<<"	";
			cout << endl;
		}
	
	cout << endl << "U matrix" << endl;
	for (i = 0; i < N; i++) 
		{
		for (j = 0; j < N; j++) cout << U[i][j] << "	";
		cout << endl;
		}
	
	cout << endl << "solutions are" << endl;
	for (i = 0; i < N; i++) cout << "x" << i + 1 << " = " << x[i] << endl;
	
	return 0;
}