//GAUSS ELIMINATION WITH PIVOTING
#include<iostream>
using namespace std;

//Almost the same function as the non pivoting case. For additional comments, see gaussel.cpp
void gauss(int n, float** lhsmat, float* rhs)
{
	int i, k, l; float sum, temp;
	for (i = 0;i < n;i++) 
	{
		//Pivoting occurs here. I am checking if a diagonal element is 0.
		if (1 + lhsmat[i][i] == 1) 
		{
			for(k = i + 1; k < n; k++)
			{
				//Here I am checking for the first row below the current row which has a non-zero element below the diagonal element.
				if (lhsmat[k][i] == 0) continue;
				//The rows are then swapped.
				for(l = 0; l < n; l++)
				{
					temp = lhsmat[i][l];
					lhsmat[i][l] = lhsmat[k][l];
					lhsmat[k][l] = temp;
				}
				temp = rhs[i];
				rhs[i] = rhs[k];
				rhs[k] = temp;
				break;
			}
		}				
		for (k = i + 1;k < n;k++) 
		{
			temp = lhsmat[k][i];
			rhs[k] = rhs[k] - rhs[i] * temp / lhsmat[i][i];
			for (l = 0;l < n;l++) 
			{
				lhsmat[k][l] = lhsmat[k][l] - lhsmat[i][l] * temp / lhsmat[i][i];
			}
		}
	}
	for(i = n - 1;i >= 0; i--)
	{
		sum = 0;
		for(k = i + 1; k < n; k++) sum+= lhsmat[i][k] * rhs[k];
		rhs[i] = (rhs[i] - sum) / lhsmat[i][i];
	}
}
/*COMMENTS ON ORDER OF OPERATIONS
For the 1st row we perform n-1 row operations (all rows below it). For the 2nd, row, n-2 operations, and so on. 
But for the 1st row, we are doing n floating point operations with each row operation, as each row has n elements. 
For the 2nd row, we are doing n-1 floating point operations with each row operations (2nd row has n-1 non-zero
elements except the one which was made 0 by the 1st row's operations) and so on.
Thus the number of operations go like n(n-1) + (n-1)(n-2) + ... +2x1, which is O(n^3)

Here I have used only partial pivoting; I swap the rows only when the diagonal element is close to 0. It will only contribute O(n^2) at most.*/

int main()
{
	int n,i,j; float *rhs,**lhsmat;
	cout<<"how many variables? - ";
	cin>>n;
	
	rhs=new float[n];
	lhsmat= new float*[n];
	for (i = 0;i < n;i++) lhsmat[i] = new float[n];

	cout<<"enter coefficient matrix"<<endl;
	for (i = 0;i < n;i++) 
	{	
		for (j = 0;j < n;j++) cin >> lhsmat[i][j];
	}
	
	cout<<endl<<"enter rhs"<<endl;
	for (i = 0;i < n;i++) cin >> rhs[i];
	
	gauss(n,lhsmat,rhs);
	
	cout<<endl<<"Triangular matrix after row operations (including pivoting):"<<endl;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++) cout<<lhsmat[i][j]<<"	";
		cout<<endl;
	}
	
	cout<<endl<<"The solutions are:"<<endl;
	for(i = 0;i < n; i++) cout<<"x"<<i+1<<" = "<<rhs[i]<<endl;

	return 0;
}

	
