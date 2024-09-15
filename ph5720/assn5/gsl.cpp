//SOLVING POISSON'S EQUATION MATRIX BY GSL'S LU DECOMPOSITION
#include<iostream>
#include<fstream>
#include<cmath>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

using namespace std;

float f(float x){ return (3.*x+(x*x))*exp(x);}

float u(float x){ return x*(1.-x)*exp(x);}  

int main(int argc, char** argv) 
{
	clock_t start,stop;
    
	int N;

	cout << "Enter no of steps: " ;
	cin >> N;

    start = clock();

    //Assigning and initializing the tridiagonal matrix.
    gsl_matrix *A;
    gsl_permutation *I;
    I = gsl_permutation_calloc(N);
    A = gsl_matrix_alloc(N,N);
    
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(i==j) 
            {
                gsl_matrix_set(A,i,j,2);
            }
            else if(i==j-1||i==j+1) gsl_matrix_set(A,i,j,-1);
            else
            {
                gsl_matrix_set(A,i,j,0);
            }
        }
    }

    float h = 1. / (N + 1.);
    int* signum;
    signum = new int[1];
    *signum = 0;

    gsl_vector *vec,*X;
    vec = gsl_vector_alloc(N);
    X = gsl_vector_alloc(N);

    for(int i=0;i<N;i++)
    {
       gsl_vector_set(vec,i,h*h*f((i+1)*h));
    } 

    gsl_linalg_LU_decomp(A,I,signum);
    gsl_linalg_LU_solve(A,I,vec,X);

    stop = clock();

    if (argc > 1)
    {
        ofstream ofile;
        ofile.open(argv[1], ios::trunc);
        ofile << "x \t u(x) Numerical \t u(x) Analytic \n";
        for (int k = 0;k < N;k++) ofile << h * (k + 1) << "\t" << gsl_vector_get(X, k) << "\t" << u((k + 1) * h) << "\t" << abs(gsl_vector_get(X, k) - u((k + 1) * h)) / (u((k + 1) * h)) << "\n";
        ofile.close();
    }
	
    cout << "Result is printed to " << argv[1] << endl;

	cout << "Time taken = " << float(stop-start)/CLOCKS_PER_SEC << " s" << endl;
	//printf("x \t solution \t error \n");
	//for(int k=0;k<N;k++) cout << h*(k+1) << "\t" << gsl_vector_get(X,k) << "\t" << abs(gsl_vector_get(X,k)-u((k+1)*h))/(u((k+1)*h)) << "\n";
	//cout<<endl;

	gsl_matrix_free(A);
	gsl_vector_free(X);
	gsl_vector_free(vec);
    gsl_permutation_free(I);

	return 0;

}
/*COMMENTS ON THE ERRORS
As expected, error decreases with decrease in step size, as can be seen in the plots. The log of the error varies linearly with x.
The gsl LU decomposition function is notably faster than the programs that we have written as it is more efficient.*/