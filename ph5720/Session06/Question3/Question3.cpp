/*=======================================================================
Program: Describe Goal of the program here
This is a double precision program to integrate an arbitrary function using thr Trapezoid, Simpson and
Gauss Quadrature Method.

Author : ROHAN IYER 
Roll No: PH16B007
Date   : 17/03/2020
*/
// NOTE: Always try to explain the method(functions or syntax) by putting comment
//       for better understanding to the examiner.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

//member function declaration
double trapezoid(double lower, double upper, int subinterval);
double simpson(double lower, double upper, int subinterval);
double gaussquad(double lower, double upper, int subinterval);
//double lagrange( double *x, double *y, double xp);

//Note: I have commented out the lagrange interpolation function because it is not required for
//the current program. However, in cases where functional form is not given and we have only data points,
//the interpolation can be used to approximate the function.

//The function. If you want to integrate a different function, change the function here.
double f(double x) { return exp(-x); };

//Main Function
//Using command line arguments to implement the 3 methods separately, as you will see below.
int main(int argc, char** argv)
{
    double lower, upper, exact; int subinterval = 2;

    cout << "Enter limits of integration: ";
    cin >> lower >> upper;
  
  /* call the function here with provided argument */
    
    //The exact value for integral of e^(-x). For different functions, it needs to be changed accordingly.
    exact = exp(-lower) - exp(-upper);

    if (argc > 1)
    {
        int cases = atoi(argv[1]);
        ofstream data;
        switch (cases)
        {
        case(1):
            //Uses the trapezoid method and stores numerical value, exact value and error in a dat file. 
            //To implenent, enter in the command line:  
            // ./a.out 1
            data.open("trapezoid_intgl.dat", ios::trunc);
            data << "#Subintls" << "\t" << "Numerical Value" << "\t\t" << "Exact Value" << "\t\t" << "Error" << endl;
            data << subinterval << "\t\t" << trapezoid(lower, upper, subinterval) << "\t\t" << exact << "\t\t" << abs((trapezoid(lower, upper, subinterval) - exact) / exact) << endl;
            for(subinterval=10; subinterval<=160; subinterval*=2)
                data << subinterval << "\t\t" << trapezoid(lower, upper, subinterval) << "\t\t" << exact << "\t\t" << abs((trapezoid(lower, upper, subinterval) - exact) / exact) << endl;

            data.close();
            cout << endl << "Data table is stored in trapezoid_intgl.dat" << endl;
            break;

        case(2):
            //Uses the Simpson method and stores numerical value, exact value and error in a dat file. 
           //To implenent, enter in the command line:  
           // ./a.out 2
            data.open("simpson_intgl.dat", ios::trunc);
            data << "#Subintls" << "\t" << "Numerical Value" << "\t\t" << "Exact Value" << "\t\t" << "Error" << endl;
            data << subinterval << "\t\t" << simpson(lower, upper, subinterval) << "\t\t" << exact << "\t\t" << abs((simpson(lower, upper, subinterval) - exact) / exact) << endl;
            for (subinterval = 10; subinterval <= 160; subinterval *= 2)
                data << subinterval << "\t\t" << simpson(lower, upper, subinterval) << "\t\t" << exact << "\t\t" << abs((simpson(lower, upper, subinterval) - exact) / exact) << endl;

            data.close();
            cout << endl << "Data table is stored in simpson_intgl.dat" << endl;
            break;

        case(3):
            //Uses the Gauss quadrature method and stores numerical value, exact value and error in a dat file. 
            //To implenent, enter in the command line:  
            // ./a.out 3
            data.open("gaussquad_intgl.dat", ios::trunc);
            data << "#Subintls" << "\t" << "Numerical Value" << "\t\t" << "Exact Value" << "\t\t" << "Error" << endl;
            data << subinterval << "\t\t" << gaussquad(lower, upper, subinterval) << "\t\t" << exact << "\t\t" << abs((gaussquad(lower, upper, subinterval) - exact) / exact) << endl;
            for (subinterval = 10; subinterval <= 160; subinterval *= 2)
                data << subinterval << "\t\t" << gaussquad(lower, upper, subinterval) << "\t\t" << exact << "\t\t" << abs((gaussquad(lower, upper, subinterval) - exact) / exact) << endl;

            data.close();
            cout << endl << "Data table is stored in gaussquad_intgl.dat" << endl;
            break;

        default: break;
        }
    }

    return 0;
}//End of Main Function


// member function definition
double trapezoid(double lower, double upper, int subinterval)
{
    //For small intervals, the area under the curve can be approximated as a trapezium:
    //Area = 0.5*h*(f(a) + f(a+h)) 
    //If we sum over all the steps:
    //Area = 0.5*h*(f(a) + f(a+h) + f(a+h) + f(a+2h) + ... + f(b-h) + f(b))
    //     = 0.5*h*(f(a) + 2f(a+h) + 2f(a+2h) + ... + 2f(b-h) + f(b))
  double stepSize = (upper - lower)/subinterval;
  double integration = f(lower) + f(upper);
  
  for(int i=1; i<= subinterval-1; i++)
    {
      double k = lower + i*stepSize;
      integration = integration + 2 * (f(k));
    }
  
  integration = integration * stepSize/2;
  return integration;
}


double simpson(double lower, double upper, int subinterval)
{
    //This method approximates the function between the 2 points as a parabola, by taking the midpoint.
    //Area between 2 points = (b-a)/6 * (f(a) + 4(f((a+b)/2) + f(b))
    //                      = h/3*(f(a) + 4(f((a+b)/2) + f(b))    (Since we take the midpoint, step size becomes half)

    //For several subintervals:     (Let h = twice the step size = b-a for any 2 points a and b)
    //Area = h/6*(f(a) + 4f(a+0.5h) + f(a+h) + f(a+h) + 4f(a + 1.5h) + f(a+2h) + ... + f(b-h) + 4f(b-0.5h) + f(b))
    
    //Let h'= h/2, the actual step size:
    //Area = h'/3*(f(a) + 4f(a+h') + 2f(a+2h') + ... + 4f(b-h') + f(b))
    
    //It is evident that if there is an odd multiple of h' (except the endpoints) we add 4f.
    //If an even multiple of h' (i%2==0) appears, we add 2f.

  double stepSize = (upper - lower)/subinterval;
  double integration = f(lower) + f(upper);
  
  for (int i = 1; i <= subinterval - 1; i++) 
    {
      double k = lower + i*stepSize;
      
      if(i%2==0)
	{
	  integration = integration + 2 * (f(k));
	}
      else
	{
	  integration = integration + 4 * (f(k));
	}
      
    }
  
  integration = integration * stepSize/3;
  return integration;
}

/*void lagrange( double *x, double *y, double xp)
{
    
  for(int i=1;i<=n;i++)
    {
      double p=1;
      for(int j=1;j<=n;j++)
	{
	  if(i!=j)
	    {
	      p = p* (xp - x[j])/(x[i] - x[j]);
	    }
	}
     double  yp = yp + p * y[i];
    }
  cout<< endl<<"Interpolated value at "<< xp<< " is "<< yp;
}*/

double gaussquad(double lower, double upper, int subinterval)
{
    //Two point gauss quadrature formula:
    //Area = (b-a)/2 * (f((a+b)/2 + (b-a)/(2*sqrt(3))) + f((a+b)/2 - (b-a)/(2*sqrt(3))))
    
    //For several intervals:
    //Area = h/2 * Sum( f(xi + 0.5h + 0.5h/sqrt(3)) + f(xi + 0.5h - 0.5h/sqrt(3))) where xi run from a to b-h.

    double stepSize = (upper - lower) / subinterval;
    double integration = 0;

    for (int i = 0; i < subinterval; i++)
    {
        double k = lower + i * stepSize;
        integration += (f(k + stepSize / 2 + stepSize / 2 / sqrt(3)) + f(k + stepSize / 2 - stepSize / 2 / sqrt(3)));
    }
    integration *= stepSize / 2;

    return integration;
}


/*COMMENTS ON RELATIVE ERROR
The trapezoid method approximates the function as a straight line between 2 points so it is less accurate.
It can be seen that the error decreases somewhat slowly with decreasing step size.

Simpson method approximates the function as a parabola so it is much more accurate. The error decreases fast with decreasing step size.

Gauss quadrature method is even more accurate with a cubic approximation ( n - point quadrature takes a 2n-1 degree polynomial approximation, here n=2) 
The error decreases rapidly even if there are only a few steps (larger step size).*/