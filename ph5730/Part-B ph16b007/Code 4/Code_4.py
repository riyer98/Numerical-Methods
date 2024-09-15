##CODE 4
##RANDOM WALK AND INTEGRATION


#I use the Metropolis Algorithm to sample points according to
#i) the Exponential distribution exp(-|x|)
#ii) Gaussian distribution with variance 2 (exp(-x^2/4))

#I will use both of these to integrate x^2 exp(-x^2/4) from -inf to inf.

#I use 100000 sample points. Step size delta is chosen such that the acceptance
#ratio is close to 1/2.
# for Exponential, delta=3
# For Gaussian, delta=4

#Sampling frequency is 10, that is, I sample at every 10th step.

#Different thermalization steps and starting points are used and I have evaluated
#the integral for those values.

from numpy import *
from matplotlib.pyplot import *


#2 cases for w(x)
def w(x,case):
    if (case==1): return 0.5*exp(-abs(x));
    if (case==2): return 1/sqrt(4*pi) *exp(-x**2/4);



#The metropolis algorithm
def metropolis(x0,N,delta,therm,freq,case):
    '''x0=starting point, N=no of sampling points, delta=step size
therm=thermalization, freq=sampling freq, case=1 or 2 (exp/gauss)'''
    #start from some x0
    x=x0;  xvals=[]; acc=0;
    for i in range(freq*N+therm):
        if(i>=therm):
            if ((i-therm)%freq==0):
                xvals.append(x);
        #trial step
        temp=x+random.uniform(-delta,delta);
        a=w(temp,case)/w(x,case);
        #accept/reject happens here
        if(a>1):
            x=temp;
            acc+=1;
        elif(a>random.random()):
            x=temp;
            acc+=1;
                
    return xvals;



#Plotting the distributions
def distplot(N,therm,freq,delta,x0,case):
     '''x0=starting point, N=no of sampling points, delta=step size
therm=thermalization, freq=sampling freq, case=1 or 2 (exp/gauss)'''
    x=metropolis(x0,N,delta,therm,freq,case);

    if (case==1): dist='Exponential exp(-|x|)';
    if (case==2): dist='Gaussian (variance=2)';

    nmin=floor(min(x))-.125;
    nmax=ceil(max(x))+.125;
    binvals=linspace(nmin,nmax,(nmax-nmin)/.25 +1);
    hist(x,bins=binvals,density='True');
    xlim([-6.5,6.5]);
    title('%s, No of points=%d, start pt= %1.3f'%(dist,N,x0));
    xlabel('x');
    ylabel('w(x)');
    show();




def f(x):
    return x**2 * exp(-x**2/4);



#Performing integration
def intg(N,therm,freq,delta,x0,case):
     '''x0=starting point, N=no of sampling points, delta=step size
therm=thermalization, freq=sampling freq, case=1 or 2 (exp/gauss)'''
    if (case==1): dist='Exponential';
    if (case==2): dist='Gaussian';

    #Generate ensemble
    x=array(metropolis(x0,N,delta,therm,freq,case));
    #We need f(x)/w(x) average.
    x=f(x)/w(x,case);

    #The integral
    I=mean(x);
    #Error = std devn
    sigma=sqrt(1/N * var(x));

    print('value of integral using %s distribution = %f'%(dist,I));
    print('error = ',sigma);

    file=open('%s_intgr.txt'%(dist),'a');
    file.write('%d\t\t%.2f\t\t%.5f\t\t%.5f\n'%(therm,x0,I,sigma));
    file.close();
    


##COMMENTS
    #The code is in fact pretty accurate and even for less thermalization and
    #far starting points gives quite accurate values.

    #Nevertheless, we can see that as we go away from the origin (which is the
    #ideal starting point because probability is maximum) eg. In Gaussian for
    #x0>50, the error begins to increase.

    #We also see that this error reduces as we increase the thermalization steps
    #as the system will return to equilibrium.

    #If we use the exponential distribution for integration, the sensitivity is
    #even lesser and the error stays more or less constant.
    #This is because the integrand x^2 exp(-x^2/4) falls off faster than exp(-x)
    #allowing the integral to converge even for far starting points.
