import numpy as np
import pylab as plt

def logistic_dist(r,N,K,nbins):
    """arguments:
    r=real number betweeen 0 and 4
    N=number of iterations of the logistic function
    K=number of times a random x0, the initial value, is picked from (0,1)
    What the function does:
    1) creates an empty array x=[].
    2) generates a random number x0 K times in a for loop.
    3) for each for loop, it appends x0 to x and runs the logistic map N times.
    4) gives an orbit plot for the first x0 generated.
    5) it gives normalised distibutions of xi's for the first 6 x0's. (which was
       asked in part a) of the questions).
    6) plots the normalised distribution for the entire ensemble of xi's that
       were generated and saves it.
    """
    
    #PART b) and c)
    x=[]
    m=0
    for i in range(K):
        x0=np.random.random()

        for j in range(N):
            x.append(x0)
            x0=r*x0*(1-x0)

        if i==0:
            plt.plot(x)
            plt.title('orbit plot for any one x0 (the first)')
            plt.show()
        if m<6:
            plt.hist(x[m*N:],bins=nbins,normed='True')
            plt.title('normalised distribution for 6 different x0')
            plt.xlabel('xi')
            plt.ylabel('p(xi)')
            plt.show()
            m+=1

    plt.hist(x,bins=nbins,normed='True')
    if r==4:
        h=np.linspace(0,1,100)
        k=1/(np.pi*(h*(1-h))**.5)
        plt.plot(h,k)
    plt.xlabel('xi')
    plt.ylabel('p(xi)')
    plt.title('normalised distribution for ensemble of iterations')
    plt.savefig('logistic_ensemble.pdf')
    plt.show()



#PART d)
def l_vs_r(rval,N,k):
    """
    arguments:
    rval: the number of values of r between 0 and 4
    N: nuumber of iterations of the logistic map
    k: number of transient states o be discarded.

    this function creates a linspace of r between 0 and 4 with rval number of
    values in between.
    it generates a random number x0 between 0 and 1.
    for each r in the linspace it generates a list x using the logistic map.
    it then discards the k transient elements and appends the remaining list to
    another list y.
    this means y is now a list of lists, with each list containing the
    non-transient elements.
    y is then plotted against r, and the plot is saved as a pdf.
    """
    y=[]
    rrange=np.linspace(0,4,rval)
    x0=np.random.random()
    for r in rrange:
        x=[x0]
        for i in range(N):
            x.append(r*x[i]*(1-x[i]))
        y.append(x[k:])

    plt.plot(rrange,y)
    plt.xlabel('r values')
    plt.ylabel('logistic map after N iterations')
    plt.title('nontransient points vs r')
    plt.savefig('logistic_rval.pdf')
    plt.show()
   
