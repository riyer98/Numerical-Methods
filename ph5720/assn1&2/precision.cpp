//Program to calculate precision of float and double

#include<iostream>
using namespace std;

int main()
{
	float eps1 = 1, preeps1=eps1; double eps2 = 1, preeps2=eps2;
	
	//This is for float. The while loop will keep making epsilon smaller until 1+eps is no longer distinguishable from 1. That is the precision of the computer
	while (1+eps1 != 1)
	{
		preeps1 = eps1;
		eps1 /= 2;
	}

	//This is for double.
	while (1+eps2 != 1)
	{
		preeps2 = eps2;
		eps2 /= 2;
	}
	
	cout << "precision of float is " << preeps1 << endl;
	cout << "precision of double is " << preeps2 << endl;

	return 0;
}