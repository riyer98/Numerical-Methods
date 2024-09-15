//A program to calculate determinant of an NxN matrix (N<=M)

#include<iostream>
#include<cmath>
using namespace std;

//The determinant function that returns the value of the determinant. It is a recursive function
//that calculates the determinant of smaller blocks in the matrix and adds them to give the result.
float det(float **A, int N)
{
	//for a 1x1 matrix the function obviously returns the number itself.
	if (N == 1) return A[0][0];
	
	//Here a temporary matrix B is defined that will store the minor matrix for a particular element.
	float** B, res = 0; int i,j,k;
	
	B = new float* [N-1];
	for (i = 0;i < N-1;i++) B[i] = new float[N-1];


	for (i = 0;i < N;i++)
	{
		for (j = 1;j < N;j++)
		{
			for (k = 0;k < N;k++)
			{
				//Here, the elements of B are added in such a way that all elements that are in
				//the same row or column as A[0][i] are skipped. The others are added.
				if (k == i) continue;
				if (k < i) B[j - 1][k] = A[j][k];
				else B[j - 1][k - 1] = A[j][k];
			}
		}
		//Here is the recursion: the function calls on itself to calculate determinant of B, then multiplies it with the
		//element A[0][i] and adds it to the result with the appropriate sign.
		res += pow(-1, i) * A[0][i] * det(B, N - 1);
	}
	return res;
}


int main()
{
	int i,j,N; float **A; 

	cout << "enter dimension of square matrix(N): ";
	cin >> N;

	A = new float* [N];
	for (i = 0;i < N;i++) A[i] = new float[N];

	//The matrix A[i][j] is entered here
	cout << "enter matrix:" << endl;
	for (i = 0;i < N;i++) 
	{
		for (j = 0;j < N;j++) 
		{
			cin>>A[i][j];
		}
	}
	//Final answer is printed here.
	cout << endl << "determinant is:" << det(A, N) << endl;
	
	for (i = 0;i < N;i++) delete[] A[i];
	delete[] A;
	
	return 0;
}
