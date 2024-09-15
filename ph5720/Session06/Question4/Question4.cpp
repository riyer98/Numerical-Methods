/*=======================================================================
Program: Describe Goal of the program here 
Program to integrate exp(x) cos(x) from 0 to 1 using Trapezoid and Simpson Method (implemented via class).
Author : ROHAN IYER 
Roll No: PH16B007
Date   : 17/03/2020
*/
// NOTE: Always try to explain the method(functions or syntax) by putting comment
//       for better understanding to the examiner.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <cmath>
using namespace std;

//The function
double f(double x) {
	return exp(x) * cos(x);
}
// Write Your Class here

class integrate
{
public:
    //Setting limits of integration and subintervals
    void setval()
    {
        cout << "Enter limits of integration: ";
        cin >> lower >> upper;
        cout << "Enter no of subintervals: ";
        cin >> subinterval;
    }

    //Exact value of integral
    void exactval()
    {
        double exact = 0.5 * (exp(upper) * (cos(upper) + sin(upper)) - exp(lower) * (cos(lower) + sin(lower)));
        cout << endl << "Exact value of integral is: " << exact << endl;
    }

    void trapezoid()
    {
        //For comments, see Q3.
        double stepSize = (upper - lower) / subinterval;
        double integration = f(lower) + f(upper);

        for (int i = 1; i <= subinterval - 1; i++)
        {
            double k = lower + i * stepSize;
            integration = integration + 2 * (f(k));
        }

        integration = integration * stepSize / 2;

        cout << endl << "Required value of integration (trapezoid method) is: " << integration << endl;
    }

    void simpson()
    {
        //For comments, see Q3.
        double stepSize = (upper - lower) / subinterval;
        double integration = f(lower) + f(upper);

        for (int i = 1; i <= subinterval - 1; i++)
        {
            double k = lower + i * stepSize;

            if (i % 2 == 0)
            {
                integration = integration + 2 * (f(k));
            }
            else
            {
                integration = integration + 4 * (f(k));
            }

        }

        integration = integration * stepSize / 3;

        cout << endl << "Required value of integration (simpson method) is: " << integration << endl;
    }

    
private:
    //In private, I have defined the limits of integration and the no of steps.
    double lower;
    double upper;
    double subinterval;
};

//Main Function
int main()
{
//Declare objects from here
    integrate f1;
    f1.setval();
    f1.exactval();
    f1.trapezoid();
    f1.simpson();

    return 0;
}//End of Main Function
