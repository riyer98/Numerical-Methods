##ASSIGNMENT 5.3
##SOLVING ISING MODEL USING MONTE CARLO METHODS

#I will be using Metropolis Algorithm to generate microstates and use them
#to compute thermodynamic quantities and plot them versus the coupling J


from numpy import *
from numpy.random import *
from scipy.special import *
from matplotlib.pyplot import *


#The function that performs the sweep
#returns change in change in energy and magnetization

def sweep(Nx,Ny,J,h,spins,r,ech):
    '''h=magnetic field, r=list of acceptance ratios,
ech= list of change in energy'''
    dE=0; dm=0;

    #the sweep
    for j in range(Nx):
        for k in range(Ny):
            s=spins[j][k];
            f=spins[j-1][k]+spins[(j+1)%Nx][k]+spins[j][k-1]+spins[j][(k+1)%Ny];
            #f is the sum of neighbouring spins. I have used
            #periodic boundary conditions.
            
            if (ech[(f,s)]<0):
                dE+=ech[(f,s)];
                dm+=-2*s;
                spins[j][k]*=-1;
            elif(random()<r[(f,spins[j][k])]):
                dE+=ech[(f,s)];
                dm+=-2*s;
                spins[j][k]*= -1;

    return (dE,dm);


#returns the energy of the lattice
def energy(spins,Nx,Ny,J,h):
    E=0;
    for i in range(Nx):
        for j in range(Ny):
            s0=spins[i][j];
            s1=spins[i-1][j];
            s2=spins[i][j-1];
            E+=-J*s0*(s1+s2)-h*s0;
    return E;


#returns magnetization of the lattice
def magnetization(spins,Nx,Ny):
    m=0;
    for i in range(Nx):
        for j in range(Ny): m+=spins[i][j];
    return m;

    

#The main function which runs the sweeps
def ising(Nx,Ny,h,J,Nsw,Ntherm):
    '''Nsw=no of sweeps
Ntherm= no of thermalization steps'''

    #initializing the lattice as random spins
    spins=choice([-1,1],(Nx,Ny));

    #listing all the 10 possible r values
    r={}; ech={};
    for i in arange(-4,5,2):
        for j in [-1,1]:
            ech[(i,j)]=2*J*i*j+2*j*h;
            r[(i,j)]=exp(-ech[(i,j)]);


    #thermalization
    for n in range(Ntherm):
        sweep(Nx,Ny,J,h,spins,r,ech);

    evals=[energy(spins,Nx,Ny,J,h)];
    mvals=[magnetization(spins,Nx,Ny)];


    #running sweeps
    for n in range(Nsw):
        dE,dm=sweep(Nx,Ny,J,h,spins,r,ech);
        evals.append(evals[-1]+dE);
        mvals.append(mvals[-1]+dm);

    evals=array(evals); 
    mvals=array(mvals);

    #averaging the quantities over all sweeps
    E=mean(evals)/Nx/Ny;
    m=abs(mean(mvals))/Nx/Ny;
    Cb=var(evals)/Nx/Ny;
    X=var(mvals)/Nx/Ny;
    
##    f1=open('therm_quants,lattice=%dx%d,sweeps=%d.txt'%(Nx,Ny,Nsw),'a');
##    f1.write("%1.1f\t\t%1.4f\t\t%1.4f\t\t%1.4f\t\t%1.4f\n"%(J,E,m,Cb,X));
##    f1.close();

    return (E,m,Cb,X);



#Here I plot the quantities and compare with the analytical expression
def plotter():
    J,E,M,C,X=loadtxt('therm_quants,b=0.txt',unpack='True');
    J1=J[1:];
    k=2*sinh(2*J1)/(cosh(2*J1))**2
    k1=2*(tanh(2*J1))**2-1
    z=exp(-2*J1);
    K1=ellipk(k)
    E1=ellipe(k);

    Ea=-J1/tanh(2*J1) * (1+2/pi *k1*K1);

    Ma=zeros(len(J1));
    for i in range (len(J1)):
        Ma[i]= (1+z[i]**2)**0.25 * (max(0,1-6*z[i]**2+z[i]**4))**0.125 / (1-z[i]**2)**0.5

    Ca= 2/pi *(J1/tanh(2*J1))**2 * (2*K1 - 2*E1 - (1-k1)*(pi/2 + k1*K1));

    plot(J1,Ea,J,E);
    xlabel('J');
    ylabel('E per particle');
    title('E vs J, 20x20 lattice, B=0, sweeps=5000');
    legend(('analytic','numerical'));
    show();

    plot(J1,Ma,J,M);
    xlabel('J');
    ylabel('M per particle');
    title('M vs J, 20x20 lattice, B=0, sweeps=5000');
    legend(('analytic','numerical'));
    show();

    plot(J1,Ca,J,C);
    xlabel('J');
    ylabel('Cb per particle');
    title('Cb vs J, 20x20 lattice, B=0, sweeps=5000');
    legend(('analytic','numerical'));
    show();

    plot(J,X);
    xlabel('J');
    ylabel('Susc per particle');
    title('Susc vs J, 20x20 lattice, B=0, sweeps=5000');
    show();




##COMMENTS
    #Since it is not practical to generate all the 2**N microstates, We use
    #the metropolis algorithm to generate microstate with probability equal
    #to the Boltzmann weight.

    #The analytical expression are for an infinite lattice, using a finite
    #lattice for computation introduces an error. We can see that the numerical
    #plots deviate slightly from the analytical plots

    #To ignore edge/boundary effects, we have used periodic boundary conditions.


    #Note: We can see that the specific heat and susceptibility have a peak
    #approximately at J=Jc=0.44. This indicates that the system undergoes a 2nd
    #order phase transition at this point.

    #For J<Jc M=0 (for B=0), so we have the paramagnetic
    #phase. For J>Jc, we have the ferromagnetic phase as we have non-zero
    #magnetization for B=0.
