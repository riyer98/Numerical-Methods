                #PH5500 ASSIGNMENT 1
                #BY ROHAN IYER - PH16B007


import numpy as np;
import matplotlib.pyplot as plt;



#Function that returns period of an array of numbers

def pertraj(x):
    for i in range(1,len(x)):
        if(round(x[0],8)==round(x[i],8)):   #rounding off numbers as machine sometimes gives wrong result due to limited precision
            return i;
    return len(x);



#Code to generate the bifurcation diagram

def bifdiag(mui,muf,mustep,N):
    '''mustep here is the number of mu values in between mui and muf.
        N is the number of iterations of the map for each mu.'''

    x0=np.random.random(); #initial condition is a random no.
    xeq=[];mu=[];

    for r in np.linspace(mui,muf,mustep+1):
        x=[];t=x0; 

        for i in range(N): #run the map N times until the values settle into its stable orbit.
            x.append(t);    
            t=r*t*(1-t);

        mu.append(r);   #creating a list of mu values
        xeq.append(x[-500:]);   #the list of numbers forming the orbit for each mu(taking final 500 iterates)

    #plotting
    plt.plot(mu,xeq,'kp',ms=0.01);
    plt.xlabel("mu");
    plt.ylabel("x");
    plt.title("Bifurcation diagram Total");
    plt.show();


#Code to calculate the first Feigenbaum constant(alpha=4.6692...)
#Code is not perfect. It comes to within 10% of the actual value.
#It gives different values when run multiple times as the initial condition is a random number. 
def feigalpha(mui,muf,mustep,N):

    '''mustep here is the number of mu values in between mui and muf.
        N is the number of iterations of the map for each mu.
        Gives optimum result for mu range 3.55 to 3.6'''

    x0=np.random.random();
    mu=[]; period=[]; feig=[];
    
    for r in np.linspace(mui,muf,mustep+1):
        x=[];t=x0; 

        for i in range(N):  #iterate the map N times
            x.append(t);
            t=r*t*(1-t);
            
        mu.append(r);   #list of mu values
        period.append(pertraj(x[-500:]));   #list of the period for each mu.


    n=64;
    while(n!=8):    #scanning the period list to see when a particular period is appearing for the first time.
                    #that corresponding mu value is inserted into another list 'feig'.
                    #I am taking mu values for period 64,32 and 16.
        for i in range(mustep+1):
            if(period[i]==n):
                feig.append(mu[i]);
                break;
        n/=2;

    print("Feigenbaum const (alpha) is ",(feig[1]-feig[2])/(feig[0]-feig[1]));



#Code to calculate second Feigenbaum constant (delta = -2.506...)

def feigdelta(mui,muf,mustep,N):

    '''mustep here is the number of mu values in between mui and muf.
        N is the number of iterations of the map for each mu.
        Gives optimum result for mu range 3.55 to 3.6'''

    x0=np.random.random(); delta=[];
    for r in np.linspace(mui,muf,mustep+1):
        x=[];t=x0; 

        for i in range(N):
            x.append(t);
            t=r*t*(1-t);
            
        T=pertraj(x[-500:]);

        if(T==64 or T==32): #scanning mu values with period 32 or 64
            eqvals=x[-500:][:T];
            for i in eqvals:
                if(abs(i-0.5)<0.0001):  #checking if any of the stable values is close to 0.5
                                        #(not exact as machine precision is finite)
                    eqvals.remove(i);
                    delta.append((-1)**(np.log2(T)+1)*min(abs(np.array(eqvals)-i)));    #taking the minimum of the distance between 0.5 and other fixed points
                    break;

    print("Feigenbaum const(delta) is",delta[0]/delta[1]);



#Code to calculate the Lyapunov exponent for a given mu value

def lyapunov(mu,N):

    '''N=no of iterations of the map'''
    
    x=np.random.random();
    lexps=0;
    for i in range(N):
        x=mu*x*(1-x);
        lexps+=np.log(abs(mu*(1-2*x)));

    return lexps/N;
    



#Code to plot Lyapunov exponent for a given range of mu

def lyaplot(mui,muf,mustep,N):

    '''mustep here is the number of mu values in between mui and muf.
        N is the number of iterations of the map for each mu.'''
    
    lyaexp=[]; mu=[];

    for r in np.linspace(mui,muf,mustep+1):
        mu.append(r);
        lyaexp.append(lyapunov(r,N));

    plt.plot(mu,lyaexp);
    plt.grid();
    plt.xlabel("mu value");
    plt.ylabel("Lyapunov exponent");
    plt.title("Lyapunov exponent Plot");
    plt.show()




##COMMENTS ON THE BEHAVIOUR OF THE MAP
##    1)  For mu<1, x=0 is a stable fixed point and the map goes to 0.
##
##    2)  at mu=1, transcritical bifurcation happens and the new fixed point is 1-1/mu. It is a stable fixed point and
##          the Lyapunov exponent is negative.
##
##    3)  mu=2 is a super-stable trajectory and we can see that the Lyapunov exponent goes to minus infinity.
##
##    4)  at mu=3, pitchfork bifurcation (period doubling) begins to takes place.
##
##    5)  until mu~3.57, we have periodic stable trajectories with super-stable trajectories in between (marked by
##        dips in the Lyapunov exponent). The periods double with smaller and smaller gaps
##
##    6)  After the Feigenbaum point(~3.57) the period goes to infinity and the Lyapunov exponent becomes positive.
##        Thus this is the onset of chaos.
##
##    7)  Now we see bands of x values that merge with increasing mu instead of increasing mu.
##
##    8)  There are windows of stable periodic behaviour. The order in which they appear (from the first period doubling
##        to the last band merging) form the Sarkovski sequence.
##
##    9)  Period 3 is the last period that appears in the bifurcation diagram.
##
##    10) The bifurcation diagram shows self similarities. I have zoomed into a part of the plot where bifurcation takes
##        place and we can see copies of the entire map present within the map.

    
