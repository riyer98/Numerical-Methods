#include<iostream>
#include<math.h>
using namespace std;

/* Define function here  like this*/
#define f(x) 1/(1+pow(x,2))

//member function declearation
void trapezoid(double lower, double upper, int subinterval);
void simpson(double lower, double upper, int subinterval);
void lagrange( double *x, double *y, double xp);


// main function
int main()
{
  
  /* call the function here with provided argument */
  
}

// member function defination
void trapezoid(double lower, double upper, int subinterval)
{
  
  double stepSize = (upper - lower)/subinterval;
  double integration = f(lower) + f(upper);
  
  for(int i=1; i<= subinterval-1; i++)
    {
      double k = lower + i*stepSize;
      integration = integration + 2 * (f(k));
    }
  
  integration = integration * stepSize/2;
  
  cout<< endl<<"Required value of integration is: "<< integration<<endl;
  
}


void simpson(double lower, double upper, int subinterval)
{
  double stepSize = (upper - lower)/subinterval;
  double integration = f(lower) + f(upper);
  
  for(int i =1; i<= subinterval-1; i++)
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
  
  cout<< endl <<"Required value of integration is: "<< integration<<endl;
}

void lagrange( double *x, double *y, double xp)
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
}
