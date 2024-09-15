//A program to print the smallest missing positive integer from an array

#include<iostream>
using namespace std;

int main()
{
	//m is the result we want. I have inititalized it to 1.
	int N, i, m = 1;
	
	cout << "enter length of array:";
	cin >> N;

	int A[N];
	cout << "enter array of integers(put space between nos):";
	for (i = 0;i < N;i++) cin >> A[i];
	
	i = 0;
	
	//This loop starts with m=1 and checks if it is equal to any number in the array.
	//If it is, then m is increased by 1 as that number is not missing. The loop is run again.
	//Once we have found an m value that is not equal to any number in the array, the loop is
	//broken and the result is printed.
	while (i != N)
	{
		//i is the counter which starts from 0. The below 'for' loop breaks if m is equal to any
		//array element and i doesn't reach N. If it does, it means m was not equal to any number
		//in the array. Then the while loop is broken.
		for (i = 0;i < N;i++)
		{
			if (m == A[i])
			{
				m++; 
				break;
			}
		}
	}
	cout << "smallest positive missing integer:" << m << endl;
	
	return 0;

}