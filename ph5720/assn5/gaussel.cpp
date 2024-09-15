//GAUSS ELIMINATION WITHOUT PIVOTING
#include<iostream>
#include "time.h"
using namespace std;

//The function that does the gauss elimination.
void gauss(int n, float** lhsmat, float* rhs)
{
	int i, k, l; float sum;
	for (i = 0;i < n;i++) 
	{
		for (k = i + 1;k < n;k++) 
		{
			float temp = lhsmat[k][i];
			//Rhs column is like the last column of the augmented matrix. I am performing row operations here.
			rhs[k] -= rhs[i] * temp / lhsmat[i][i];	
			for (l = 0;l < n;l++) 
			{
				//Row operations on the matrix elements to make the lower triangular elements 0. Hence the matrix will be an upper triangular matrix after this.
				lhsmat[k][l] -= lhsmat[i][l] * temp / lhsmat[i][i];
			}
		}
	}
	//Now. back substitution, going from x(n-1) to x(0).
	for(i = n - 1; i >= 0; i--)
	{
		sum = 0;
		//Defined the variable sum. x[i] = rhs[i] - Sum(k = 0 to i-1)lhsmat[i][k]*x[k];
		//Note: I have not defined an extra array to store the answers. I am directly operating on the rhs and storing it as the answer.
		for(k = i + 1; k < n; k++) sum+= lhsmat[i][k] * rhs[k];
		rhs[i] = (rhs[i] - sum) / lhsmat[i][i];
	}
}
	
int main()
{
	int n, i, j; float* rhs, ** lhsmat;
	cout << "how many variables? - ";
	cin >> n;
	
	rhs = new float[n];
	lhsmat = new float* [n];
	for (i = 0; i < n; i++) lhsmat[i] = new float[n];

	cout<<"enter coefficient matrix"<<endl;
	for (i = 0; i < n; i++) 
	{	
		for (j = 0; j < n; j++) cin >> lhsmat[i][j];
	}
	
	cout<<endl<<"enter rhs"<<endl;
	for (i = 0; i < n; i++) cin >> rhs[i];
	
	gauss(n, lhsmat, rhs);
	
	cout<<endl<<"Triangular matrix after row operations:"<<endl;
	for (i = 0; i < n; i++) 
	{
		for (j = 0; j < n; j++) cout << lhsmat[i][j] << "	";
		cout << endl;
	}
	
	cout << endl << "The solutions are:" << endl;
	for (i = 0; i < n; i++) cout << "x" << i + 1 << " = " << rhs[i] << endl;

	return 0;
}