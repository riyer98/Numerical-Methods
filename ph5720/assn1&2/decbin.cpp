//This program converts both decimal to binary and vice versa. 
//Note: these programs only work for small integers, not for floating numbers or integers larger than about 1000 
//I tried doing the float part but there was an error while the program took in the values.
//eg. when i inputed 2.01, the computer took it as something like 2.0100004, which gave an error in the final result.

#include<iostream>
#include<cmath>
using namespace std;

//First program: to convert decimal to binary.
int dectobin(int x)
{
	//Here I have initialised a local variable k to the input value x.
	int i = 0, k = x, ans = 0;

	//Here I am dividing the input continuously by 2 and collecting the remainders (k%2) which will give me the number in binary
	//e.g. 10 in binary is 1010 =  1*2^3 + 0 * 2^2 + 1*2^1 + 0*2^0
	//the digits are basically remainders when 10 is divided continuosly by 2.
	//but I want to display it as 1010, so I expand in powers of 10:
	//1*10^3 + 0*10^2 + 1*10^1 + 0*10^0    (that's what the pow(10,i++) does)
	//I am displaying it as a "fake decimal integer", because 2 in binary is nothing but 10.
	while (k!=0)
	{
		ans+= (k%2)* pow(10,i++);
		k/=2;
	}

	return ans;
}

//Second program: converts binary to decimal
int bintodec(int x)
{
	int i = 0, k = x, ans = 0;
	
	//It's the opposite here: divide by 10 and get the remainders
	//We expand not in powers of 10, but powers of 2 (which is what pow(2,i++) does), because the number is in binary
	while (k!=0)
	{
		ans+= (k%10)* pow(2,i++);
		k/=10;
	}
	
	return ans;
}


int main()
{
	int n,x;

	cout<<"Enter: 1 for decimal to binary"<<endl;
	cout<<"       2 for binary to decimal"<<endl;
	
	cin>>n;
	
	//I have created 2 cases: if we enter 1 for the input, the program will call decimal to binary function
	//if 2 is entered, the program will call the binary to decimal function. The default will exit from the program
	switch(n)
	{
		case(1): 
		{
			cout<<"enter number:";
			cin>>x;
			cout<<dectobin(x)<<endl;
			break;
		}
		case(2):
		{
			cout<<"enter number:";
			cin>>x;
			cout<<bintodec(x)<<endl;
			break;
		}
		default: return 0;
	}
}