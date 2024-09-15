//A program to calculate e^-x using 3 different methods

#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//For the power series we need n!, so I have defined a recursive function that 
//returns n(n-1)! until it reaches 0!=1, so it gives n(n-1)(n-2)...1
//I have defined it as double, otherwise it reaches NaN pretty quickly as factorial grows rapidly.
double factorial(double m)
{
	if (m == 0) return 1;
	else return m * factorial(m - 1);
}


//First method: adding (-x)^n/n! to a separate variable and returning the answer.
double method1(float x,int N)
{
	double n, ans = 0;
	
	for (n = 0;n <= N;n++)
	{
		ans += pow(-x, n) / factorial(n);
	}
	return ans;
}
//Comments on large N behaviour and variation with x value
//For small x values, as N inceases, the series converges to the actual value, so the error decreases to zero. 
//However, for large x, the series does not converge fast enough as each term x^n/n! becomes very large, 
//and due to the limited precision of the computer it becomes a constant value and if increases further, gives "NaN" due to the factorial function blowing up..


//Second method: using the sequence s(0)=1 and s(n)= -s(n-1)*x/n, we can get s(n) = (-x)^n/n!
//We sum this series to get e^-x
double method2(double x,int N)
{
	int n; double s = 1, sum = 1;
	for (n = 1;n <=N;n++)
	{
		s = -s * x / n;
		sum += s;
	}
	return sum;
}
//Comments on large N behaviour and variation with x value
//Similar to method 1, the error goes to zero for small x values but goes to a non-zero
//constant for large x. Unlike method 1 however, it does not give "NaN" as it does not compute n factorial but uses recursion.



//Third method: I calculate e^x using s(0)=1 and s(n)=s(n-1)*x/n then take reciprocal to get e^-x
double method3(double x, int N)
{
	int n; float p = 1, res = 1;

	for (n = 1;n <= N;n++)
	{
		p = p * x / n;
		res += p;
	}
	return 1/res;
}
//Comments on large N behaviour and variation with x value
//This method is different from the other two because it is actually computing e^x and returning the recciprocal.
//For small x the error goes to zero, for large x, the function returns 0 (as e^x -> inf, e^-x -> 0), so unlike the other 
//two methods, the error doesn't blow up but becomes 1 
//( err = |res(cal) - res|/res); res(cal) = 0 hence err = 1)


int main()
{
	double x; int N, m; ofstream x0;

	//Value of the exponent is taken as input
	cout << "Enter x value:";
	cin >> x;
	

	//the below command (commented out) is used to create a .dat file
	//x0.open("exp_mx_data.dat", ios::trunc);

	//Just printing the heading for each column. x0 prints it in the dat file whereas cout prints it on the screen
	x0 << "#Truncation Value N"<<"      "<<"Method 1" << "   " << "Method2"<<"   "<<"   "<<"Method3"<<"    "<<"Relative Error1"<<"   "<<"Relative Error2"<<"    "<<"Relative Error3"<<endl;
	cout<< "Truncation Value N" << "      " << "Method 1" << "   " << "Method2" << "   " << "   " << "Method3" << "    " << "Relative Error1" << "   " << "Relative Error2" << "    " << "Relative Error3" << endl;
	
	
	//This is a loop that calculates e^-x and the relative error for different truncation(N) values ranging from 10 to 300 in steps of 10.
	//There are 7 columns printed. First the N value, then the next 3 are the value of e^-x using the 3 different methods, and finally the errors due to all 3.
	for (N = 10;N <= 300;N+=10)
	{
		x0 << "       " << N << "             " << method1(x, N) << "    " << method2(x, N) << "     " <<"      "<<method3(x,N)<<"     "<< abs((method1(x, N) - exp(-x)) / exp(-x)) << "   " << abs((method2(x, N) - exp(-x)) / exp(-x)) << "    "<<abs((method3(x,N)-exp(-x))/exp(-x))<<endl;
		cout << "       " << N << "             " << method1(x, N) << "    " << method2(x, N) << "     " << "      " << method3(x, N) << "     " << abs((method1(x, N) - exp(-x)) / exp(-x)) << "   " << abs((method2(x, N) - exp(-x)) / exp(-x)) << "    " << abs((method3(x, N) - exp(-x)) / exp(-x)) << endl;
	}
	
	//The below command closes the file (commented out as file has already been created.)
	//x0.close();
	
	return 0;
}
