##ASSIGNMENT 5.1 AND 5.2
##MONTE CARLO AND METROPOLIS INTEGRATION 

from numpy import *
from matplotlib.pyplot import *


#Q 5.1 (von Neumann sampling)
def f(x):
    return 1/(1+x**2);

def w(x):
    '''the weight function'''
    return (4-2*x)/3;
    

def wpr():
    '''w'(x), used for von Neumann acceptance/rejection'''
    return 4/3;



def mc_intg(xmin,xmax,N,case):
    '''xmin=lower limit, xmax=upper limit,
N=no of points sampled, case -> enter 0(without weight) or 1(with weight)'''
    x=random.uniform(xmin,xmax,N);
    if (case==0):           #uniform sampling
##        hist(x,bins=50,density='True');
##        xlabel('x');
##        ylabel('w(x)');
##        title('Points sampled uniformly');
##        show();

        x=f(x);
        xintg=(xmax-xmin)*mean(x);
        #print('Value of integral is %f'%(xintg));
        return xintg;

    elif (case==1):         #sampling with weights.
        y=[];
        for i in range(N):
            if (random.random()<w(x[i])/wpr()):
                y.append(x[i]);

##        hist(y,bins=50,density='True');
##        xlabel('x');
##        ylabel('w(x)');
##        title('Points sampled according to w(x) density');
##        show();

        y=array(y);
        y=f(y)/w(y);
##        print('no of points accepted = %d'%(len(y)));
##        print('Value of integral is %f'%(mean(y)));

        return mean(y);
    else: return 0;



def varcomp(xmin,xmax,N,nint):
    '''xmin=lower limit, xmax=upper limit
N=no of sampling points, nint=no of times integral is calculated
This function plots and compares the deviation in the computation'''
    int0=[]; int1=[];
    for n in range(nint):
        int0.append(mc_intg(xmin,xmax,N,0));
        int1.append(mc_intg(xmin,xmax,N,1));

    hist(int0,bins=50,density='True');
    hist(int1,bins=20,density='True');
    xlabel('x')
    ylabel('Normalized Frequency');
    title('Comparing Deviations,no of samples=%d'%(N));
    legend(('Without w(x)','With w(x)'));
    show();

    print('mean value of integral without w(x)= ',mean(int0));
    print('standard deviation= ',std(int0));
    print('mean value of integral with w(x)= ',mean(int1));
    print('standard deviation= ',std(int1));
    
        



#Q 5.2 (Metropolis sampling)
def f2(x):
    return x**2 * exp(-x**2/2);

def w2(x):
    return exp(-x**2/2)/(2*pi)**0.5;


def C(k,x):
    '''Computes autocorrelation C(k) for array x and distance k'''
    s=0;
    for i in range(k,len(x)):
        s+=x[i]*x[i-k];

    s=(s/(len(x)-k)-(mean(x))**2)/var(x);
    return s;


def metropolis(x0,N,delta,freq,therm):
    '''x0=starting point, N=no of sampling points, delta=step size'''
    x=x0;  xvals=[]; acc=0;
    for i in range(freq*N+therm):
        if(i>=therm):
            if ((i-therm)%freq==0):
                xvals.append(x);
        #trial step
        temp=x+random.uniform(-delta,delta);
        a=w2(temp)/w2(x);
        #accept/reject happens here
        if(a>1):
            x=temp;
            acc+=1;
        elif(a>random.random()):
            x=temp;
            acc+=1;

    return xvals;


def corr_plot(x0,N,delta,freq,therm):
    '''x0=starting point, N=no of sampling points, delta=step size
Plots autocorrelation as function of k for given step size delta.'''
    x=array(metropolis(x0,N,delta,freq,therm));
    
    corr=[];
    for k in arange(0,100,1):
        corr.append(C(k,x));
    corr=array(corr);
    
    savetxt('correlation,delta=%1.4f.dat'%(delta),corr.T,fmt='%f',header='Autocorrelation, step size=%1.4f \n'%(delta));
    print('data saved in txt file');
    plot(arange(0,100*freq,freq),corr);
    xlabel('k');
    ylabel('C(k)');
    title('Autocorrelation, step size=%1.4f'%(delta));
    show();



def metro_intg(x0,N,delta,freq,therm):
    '''x0=starting point, N=no of sampling points, delta=step size
Evaluates the integral of the f2(x) using the points sampled from
metropolis algorithm.'''
    x=array(metropolis(x0,N,delta,freq,therm));

    xmin=floor(min(x));
    xmax=ceil(max(x));
    bs=4*(xmax-xmin)+2;
    hist(x,bins=linspace(xmin-0.125,xmax+.125,bs),density='True');
    xlabel('x');
    ylabel('w(x)');
    title('w(x) sample, no of pts=%d, start pt=%1.3f'%(len(x),x0));
    show();

    x=f2(x)/w2(x);
    print('Value of integral x exp(-x**2/2) from -inf to inf is ',mean(x));
