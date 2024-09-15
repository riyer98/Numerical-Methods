/*=======================================================================
Program: Describe Goal of the program here 
To read a file and interpolate the data in that file.

Author : ROHAN IYER
Roll No: PH16B007
Date   : 17/03/2020
*/
// NOTE: Always try to explain the method(functions or syntax) by putting comment
//       for better understanding to the examiner.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <cmath>
#include<fstream>
using namespace std;

//Lagrange interpolation function
double lagrange(double* x, double* y, double xp, int n)
{
	//For n points, we approximate it as an n-1 degree polynomial.
	//f(x) = Sum over i[ yi* product over j!=i { (x-xj)/(xi-xj) }]
	//where yi = f(xi)
	double yp = 0;
	for (int i = 0;i < n;i++)
	{
		double p = 1;
		for (int j = 0;j < n;j++)
		{
			if (i != j)
			{
				p = p * (xp - x[j]) / (x[i] - x[j]);
			}
		}
		yp += p * y[i];
	}
	return yp;
}

//Main Function
int main()
{
	char line[70]; FILE* fpr;
	//Opening the file for reading
	fpr = fopen("sample_data.dat", "r");
	fgets(line, sizeof(line), fpr);
	
	int n = 0, i = 0;
	//This is to determine how many lines there are in the file, i.e how many data points.
	//The number n will increment until it hits NULL, thus n becomes the no of data points.
	while (fgets(line, sizeof(line), fpr) != NULL) n++;
	
	double* x_old = new double[n];
	double* f_old = new double[n];
	
	//We need to reopen the file because we had read it previously to determine n.
	fpr = fopen("sample_data.dat", "r");
	fgets(line, sizeof(line), fpr);

	while (fgets(line, sizeof(line), fpr) != NULL)
	{
		sscanf(line, "%lf	%lf", &x_old[i], &f_old[i]);
		i++;
	}
	int N = int((x_old[n - 1] - x_old[0]) / 2) + 1;
	double* x_new = new double[N];
	double* f_new1 = new double[N];
	double* f_new2 = new double[N];
	
	int j = 0;
	for (i = 0; i<N; i++)
	{
		//Interpolation using all 9 points
		//x_new is defined for 2 MeV intervals
		x_new[i] = 2 * i;
		f_new1[i] = lagrange(x_old, f_old, x_new[i], n);
		
		if (x_new[i] > x_old[j + 2]) j += 2;
		
		//Interpolation using 3 points at a time
		double xol[3] = { x_old[j], x_old[j + 1], x_old[j + 2] }, fol[3] = { f_old[j], f_old[j + 1], f_old[j + 2] }; 
		f_new2[i] = lagrange(xol, fol, x_new[i], 3);
	}

	//Storing data in a dat file.
	ofstream xvsf;
	xvsf.open("interpolated_data.dat", ios::trunc);
	xvsf << "#x" << "\t" << "interpolation using all 9 points" << "\t" << "Interpolation 3 at a time" << endl;
	for (i = 0; i<N; i++)
		xvsf << x_new[i] << "\t\t\t" << f_new1[i] << "\t\t\t" << f_new2[i] << endl;

	xvsf.close();
	cout << "Interpolated data has been stored in 'interpolated_data.dat'" << endl;

	return 0;
}//End of Main Function
