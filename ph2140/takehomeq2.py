import numpy as np
import pylab as plt

def unif_randwalk(T,N):
    """
    arguments:
    T: the number of steps of the walker,i.e, the time
    N: the number of times the random walk is generated.

    this function has 2 parts:
    PART I)
    here two lists x and y are initialised as [0]. then with equal probability,
    either x increases or decreases by unity or y increases or decreases by
    unity. the probability of each movement here is 1/4.
    the values are then appended to x and y, and their path is plotted.

    PART II)
    there are a 4 arrays of zeroes. the random walk is carried out N times,
    each time with T steps. the sums of square of distance (x**2+y**2)are
    assigned to sumrsquare, sum of distances to sumr, the standard deviation for
    each T to stddev. the distance of the endpoints (position of each of the N
    walker after T steps) is appended to a list rN.

    the mean and std deviation are displayed by using np.mean(rN) and np.std(rN).
    the normalised distribution of r is plotted. next, <r>, sigma and <r^2> are
    all plotted in one plot as a function of T.. <r^2> in red, <r> in green and
    sigma in blue.

    the mean of r^2 = <r^2>, in fact, comes out to be close to T for large N.
    thus <r^2> =D*t^a, where a=1. this shows that r^2 for a random walk is
    proportional to t. 
    since <r^2> and <r>^2 are both proportional to T, sigma comes out to be
    proportional to T**0.5.
    D here is the length of each step taken, i.e.= 1. if x or y increased or
    decreased by 2, for instance, then D would be equal to 2.
    """
    x=[0]
    y=[0]
    for i in range(T):
        j=np.random.randint(0,4)
        if j==0:
            x.append(x[i]+1)
            y.append(y[i])
        elif j==1:
            x.append(x[i]-1)
            y.append(y[i])
        elif j==2:
            y.append(y[i]+1)
            x.append(x[i])
        else:
            y.append(y[i]-1)
            x.append(x[i])
    plt.plot(x,y)
    plt.show()

    sumr=np.zeros(T)
    sumrsquare=np.zeros(T)
    rN=[]
    stddev=np.zeros(T)
    for k in range(N):
        x=np.zeros(T)
        y=np.zeros(T)

        for i in range(T-1):
            j=np.random.randint(0,4)
            if j==0:
                x[i+1]=x[i]+1
                y[i+1]=y[i]
            elif j==1:
                x[i+1]=x[i]-1
                y[i+1]=y[i]
            elif j==2:
                y[i+1]=y[i]+1
                x[i+1]=x[i]
            else:
                y[i+1]=y[i]-1
                x[i+1]=x[i]
        rN.append((x[T-1]**2 + y[T-1]**2)**.5)
        sumrsquare+=x**2 + y**2
        sumr+=(x**2+y**2)**0.5
        stddev=(sumrsquare/N - (sumr/N)**2)**0.5
        
        
    print('mean of r = <r> =',np.mean(rN))
    print('std deviation = sigma =',np.std(rN))
    
    plt.hist(rN,bins=100,normed='True')
    plt.xlabel('r')
    plt.ylabel('p(r)')
    plt.title('distribution of r')
    plt.show()
    plt.plot(sumrsquare/N,'r')
    plt.plot(sumr/N,'g')
    plt.plot(range(T),'y')
    plt.plot(stddev,'b')
    plt.xlabel('time T')
    plt.ylabel('<r^2>,<r> and sigma')
    plt.title('variation of <r^2>,<r> and sigma with T')
    plt.show()
    


def nonunif_randwalk(T,N):
    """
    arguments:
    T: the number of steps of the walker,i.e, the time
    N: the number of times the random walk is generated.

    this function has 2 parts:
    PART I)
    here two arrays x and y are initialised as [0]. then it increments either
    x[i] or y[i] by a number from the standard normal distribution, with equal
    probability. the values are then appended to lists x and y, and their path is
    plotted.

    PART II)
    the random walk is done N times here. here r = distance from origin after T
    steps = sqrt(x[T]^2 + y[T]^2). the mean and standard deviation of the r's is
    then printed.

    the mean of r^2 = <r^2>, in fact, comes out to be close to T for large N.
    thus <r^2> =D*t^a, where a=1. this shows that r^2 for a random walk is
    proportional to t.
    """
    x=[0]
    y=[0]
    for i in range(T):
        j=np.random.randint(0,2)
        if j==0:
            x.append(x[i]+np.random.standard_normal())
            y.append(y[i])
        else:
            y.append(y[i]+np.random.standard_normal())
            x.append(x[i])
    plt.plot(x,y)
    plt.show()
    
    sumr=np.zeros(T)
    sumrsquare=np.zeros(T)
    stddev=np.zeros(T)
    rN=[]
    for j in range(N):
        x=np.zeros(T)
        y=np.zeros(T)
        for i in range(T-1):
            k=np.random.randint(0,2)
            if k==0:
                x[i+1]=x[i]+np.random.standard_normal()
                y[i+1]=y[i]
            else:
                y[i+1]=y[i]+np.random.standard_normal()
                x[i+1]=x[i]
        
        sumr+=(x**2 + y**2)**0.5
        sumrsquare+=x**2 + y**2
        rN.append((x[T-1]**2 + y[T-1]**2)**0.5)
        stddev=(sumrsquare/N - (sumr/N)**2)**.5

    print('mean of r = <r> =',np.mean(rN))
    print('std deviation = sigma =',np.std(rN))
    plt.hist(rN,bins=100,normed='True')
    plt.xlabel('r')
    plt.ylabel('p(r)')
    plt.title('normalised distribution of r')
    plt.show()
    plt.plot(sumrsquare/N,'r')
    plt.plot(sumr/N,'g')
    plt.plot(stddev,'b')
    plt.plot(range(T),'y')
    plt.xlabel('time T')
    plt.ylabel('<r^2>,<r>,sigma')
    plt.title('variation of <r^2>,<r> and sigma with T')
    plt.show()

def tg_unif_rwalk(T):
    """
    uniform random walk with turtle grapics."""
    import turtle
    for i in range(T):
        j=np.random.randint(0,4)
        if j==0:
            turtle.right(90)
            turtle.forward(20)
            turtle.left(90)
        elif j==1:
            turtle.left(90)
            turtle.forward(20)
            turtle.right(90)
        elif j==2:
            turtle.forward(20)
        else:
            turtle.backward(20)
            
def tg_nonunif_rwalk(T):
    """
    non uniform random walk with turtle graphics"""
    import turtle
    for i in range(T):
        j=np.random.randint(0,2)
        if j==0:
            turtle.right(90)
            turtle.forward(20*np.random.standard_normal())
            turtle.left(90)
        else:
            turtle.forward(20*np.random.standard_normal())
