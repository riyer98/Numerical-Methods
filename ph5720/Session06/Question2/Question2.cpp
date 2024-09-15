/*=======================================================================
Program: Describe Goal of the program here 
A program to integrate 1/2+x^2 from 0 to 3 using trapezoid method.
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

//member function declaration
double trapezoid(double lower, double upper, int subinterval);

double f(double x) { return 1 / (2 + x * x); }

//Main Function
int main()
{
    int * N, i; double upper, lower, * res, * err, exact;

    cout << "Enter limits of integration: ";
    cin >> lower >> upper;

    exact = atan(3. / sqrt(2)) / sqrt(2);
    
    //I have initialised the arrays to size 15 because I will take 15 N values, starting from 2^0 (1) to 2^14 (16384), 
    //doubling the N value in each step.
    N = new int[15];
    res = new double[15];
    err = new double[15];
   
    for (i = 0; i < 15; i++)
    {
        N[i] = pow(2, i);
        res[i] = trapezoid(lower, upper, N[i]);
        err[i] = abs((res[i] - exact) / exact);
    }

    //Storing data in valanderr.dat
    ofstream data;
    data.open("valanderr.dat",ios::trunc);
    
    data << "#N(steps)" << "\t" << "Value (trapezoid)" << "\t" << "Exact value" << "\t" << "Error" << endl;
    for (i = 0; i < 15; i++) data << N[i] << "\t\t" << res[i] << "\t\t" << exact << "\t\t" << err[i] << endl;
    data.close();

    return 0;
}//End of Main Function

// member function definition
double trapezoid(double lower, double upper, int subinterval)
{
  //For comments on the trapezoid method, see Question 3.
    double stepSize = (upper - lower) / subinterval;
    double integration = f(lower) + f(upper);
  
    for (int i = 1; i <= subinterval - 1; i++) 
    {
        double k = lower + i * stepSize;
        integration = integration + 2 * (f(k));
    }
  
    integration = integration * stepSize / 2;
    return integration;
}

/*COMMENTS ON ERROR
Error decreases as a power of the number of steps N (err proportional to N^k for some constant k)
Hence in the log log plot, it comes out as a straight line.*/