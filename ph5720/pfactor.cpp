#include<iostream>
#include<cmath>
#include<ctime>

using namespace std;

unsigned int getprimes(unsigned long int n, unsigned int * primes)  //gets a list of primes <=sqrt(n)
{
    primes[0] = 2; primes[1] = 3; 
    unsigned int j, i, k=2; 
    
        for (i=5; i<=(unsigned int)sqrt(n); i+=2)           //check every odd no i <= sqrt(n)
        {
            for(j=1; i%primes[j]!=0; j++)                   //see if i has any prime factors <= sqrt(i)
            {
                if (primes[j]>(unsigned int)sqrt(i))
                {
                    primes[k]=i;                            //if not, then i is prime and is added to the list
                    k++;
                    break;
                }
            }
        }
    return k;
}

int main()
{
    unsigned long int n; unsigned int * primes; 
    unsigned int j, n_prime; clock_t timereq;
    
    cout<<"\nEnter integer >=2: ";
    cin>>n;     //enter integer
    
    if(n<2) {cout<<"invalid entry\n";
    return 0;}
    
    timereq = clock();

    primes = new unsigned int[(unsigned int)(2.5*sqrt(n)/log(n))];

    n_prime = getprimes(n,primes);

    cout<<"\nFactors: ";
    for(j=0; n!=1 && j<n_prime; j++)
    {
        while (n%primes[j]==0)
        {
            cout<<primes[j]<<" ";       //prime factors are printed as many times as they divide n
            n/=primes[j];
        }
    }  
    if(n!=1) cout<<n; 

cout<<"\n\n";
cout<<"time elapsed = "<<clock()-timereq<<" microseconds\n\n";

return 0;
}