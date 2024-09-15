##ASSIGNMENT 1
##FINDING ENERGIES OF BOUND STATES IN A SQUARE WELL POTENTIAL
##USING UNITS WHERE MASS M=1 AND hbar = 1
import sys
import numpy as np;
import matplotlib.pyplot as plt;


def f(x,x0,parity):
    #The function whose roots we want to find. This is obtained by
    #applying boundary conditions to the wavefunction. If parity is
    #1 (even), Then it returns the first output. If it is -1 (odd),
    #it returns the second output.
    if((1+parity)/2):
        return (x*np.sin(x)- np.sqrt(x0**2 - x**2)*np.cos(x));
    else:
        return (x*np.cos(x)+ np.sqrt(x0**2 - x**2)*np.sin(x));



#getting sub-intervals which have roots
def bruteforce(a,b,N,x0,parity):
    '''Enter the interval [a,b] you want to scan. N is the number of
sub-intervals we will divide [a,b] into. The program will check for
a root in each of the sub-intervals. 
parity = +1 or -1 only. It will choose one of the 2 return outputs.
x0 is related to V0. it is defined in the sqwellenergies function.'''
    h=(b-a)/N; solns=[];
    for i in range(N):
        if(f(a+i*h,x0,parity)*f(a+(i+1)*h,x0,parity)<=0):
            solns.append([a+i*h,a+(i+1)*h]);
    return solns;



def bisection(a,b,tolerance,x0,parity):
    '''enter interval [a,b]. Note: a must be < b.
parity = +1 or -1 only. It will choose one of the 2 return outputs.
x0 is related to V0. it is defined in the sqwellenergies function.'''

    #This is a recursive function that calls upon itself every
    #time the interval is bisected.
    
    if (f(a,x0,parity)*f(b,x0,parity)>0):
        return ('Inconclusive'); #This method works only if
                                                #f(a)*f(b) <= 0.
    if (f(a,x0,parity)==0): return a;
    if (f(b,x0,parity)==0): return b;

    x=a; y=b;

    while(abs(x-y)>=tolerance):
        c=(x+y)/2;
        #bisection step happens here.
        if (f(c,x0,parity)==0): break;
        if(f(x,x0,parity)*f(c,x0,parity)<0):
            y=c;
        else:
            x=c;

    return c;



#Determining energy eigenvalues. 
def sqwellenergies(V0,a,tolerance,N):
    '''V0= depth of potential well in eV
a=half the width of the well in Angstrom
N=we will divide the interval [0,V0] into N subintervals and
use brute force to get the subintervals which contain roots.'''
    x0= np.sqrt(2*V0)*a;
    evenparity=[]; oddparity=[];
    
    intls=bruteforce(0,x0,N,x0,1);
    for z in intls:
        E=(bisection(z[0],z[1],tolerance,x0,1)/a)**2 /2;
        evenparity.append('%1.3f'%E);   #creating list of even parity energies

    intls=bruteforce(0,x0,N,x0,-1);
    for z in intls:
        E=(bisection(z[0],z[1],tolerance,x0,-1)/a)**2/2;
        if(E==0): continue; #for odd parity we get a trivial root which
                            #we want to avoid.
        oddparity.append('%1.3f'%E);    #creating list of odd parity energies

    print('Even parity energies (in eV) are: ',evenparity);
    print('Odd parity energies (in eV) are: ',oddparity);




#For plotting wavefunctions
def wfplotter(E,V0,a,parity,plttitle):
    '''E=energy eigenvalue in eV
V0=depth of potential well in eV
a=half width of well in Angstrom
parity = +1 or -1 => decides if eigenfunction is even or odd
plttitle=a string that will become the title of the plot.
Please enter plttitle in inverted commas.'''
    xvals=np.linspace(-2*a,2*a,201); yvals=[];
    k=np.sqrt(2*E); K=np.sqrt(2*(V0-E));
    if((1+parity)/2):
        A=1/np.sqrt((np.cos(k*a))**2 /K + a + np.sin(2*k*a)/(2*k));
        for x in xvals:
            if(abs(x)>a):
                yvals.append(A*np.cos(k*a)*np.exp(K*(a-abs(x))));
            else:
                yvals.append(A*np.cos(k*x));
    else:
        A=1/np.sqrt((np.sin(k*a))**2 /K + a - np.sin(2*k*a)/(2*k));
        for x in xvals:
            if(abs(x)>a):
                yvals.append(abs(x)/x * A*np.sin(k*a)*np.exp(K*(a-abs(x))));
            else:
                yvals.append(A*np.sin(k*x));

    plt.plot(xvals,yvals);
    plt.xlabel('x');
    plt.ylabel('psi(x)');
    plt.title('%s'%(plttitle));
    plt.show();
    

                        
