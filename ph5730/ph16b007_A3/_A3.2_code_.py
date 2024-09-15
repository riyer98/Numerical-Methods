##ASSIGNMENT 3.2
##FINITE DIFFERENCE METHODS FOR HEAT EQUATION

#This code contains the explicit, implicit and Crank-Nicolson method
#to solve the heat equation. Note, this code can only solve for Dirichlet
#boundary conditions. For Neumann conditions, the code will have to be
#modified.

from numpy import *
from matplotlib.pyplot import *


#First boundary condition at x = initial point (for all t)
#If you want to change the boundary condition, change this function
#with suitable arguments and outputs.
def bdrycond1(t):
    return 0;


#Second boundary condition at x = final point  (for all t)
#If you want to change the boundary condition, change this function
#with suitable arguments and outputs.
def bdrycond2(t):
    return 0;


#Initial condition for all x but t = initial time
#If you want to change initial conditions, change the function here.
def inicond(x):
    return 4*x*(1-x);



#Explicit method - takes in array of the function at all x values
#and returns the array at the next time step.
def explicit(l,u0,t):
    '''l=D(delt)/(delx)^2
u0 = array at time t
t = time'''
    u1=zeros(len(u0));

    #Dirichlet conditions
    u1[0]=bdrycond1(t);
    u1[-1]=bdrycond2(t);
    
    for i in range(1,len(u0)-1):
        u1[i]=l*u0[i-1]+(1-2*l)*u0[i]+l*u0[i+1];

    return u1;



#Explicit method - takes in array of the function at all x values
#and returns the array at the next time step.
#Solves the tridiagonal matrix to obtain array at next t step.
def implicit(l,u0,t):
    '''l=D(delt)/(delx)^2
u0 = array at time t
t = time'''
    u1=u0;
    #Adjusting for Dirichlet conditions
    #If you are using Neumann conditions, comment out these 4 lines.
    u1[0]=bdrycond1(t);
    u1[-1]=bdrycond2(t);
    u1[1]+=l*bdrycond1(t);
    u1[-2]+=l*bdrycond2(t);
    
    A=[];
    A.append(array([1+l,-l]));
    for i in range(3,len(u0)-1):
        A.append(array([-l,1+2*l,-l]));
    A.append(array([-l,1+l]));
    A=array(A);

    u1[2]-=u1[1]*A[1][0]/A[0][0];
    A[1][:2]-=A[0]*A[1][0]/A[0][0];
    for i in range(1,len(u0)-3):
        u1[i+2]-=u1[i+1]*A[i+1][0]/A[i][1];
        A[i+1][:2]-=A[i][1:]*A[i+1][0]/A[i][1];

    for i in arange(len(u0)-3,1,-1):
        u1[i]-=u1[i+1]*A[i-1][2]/A[i][1];
        A[i-1][2]=0;
    u1[1]-=u1[2]*A[0][1]/A[1][1];
    A[0][1]=0;

    u1[1]/=A[0][0];
    for i in range(1,len(u0)-2):
        u1[i+1]/=A[i][1];

    #Adjusting for Neumann conditions (derivative=0).
    #When using Dirichlet conditions, comment out these 2 lines.
##    u1[0]=u1[1];
##    u1[-1]=u1[-2];

    return u1;



#Crank-Nicolson method - takes in array of the function at all x
#values and returns the array at the next time step.
#Average of explicit and implicit methods.
def cranknicolson(l,u0,t):
    '''l=D(delt)/(delx)^2
u0 = array at time t
t = time'''
    u1=zeros(len(u0));

    #Dirichlet conditions
    u1[0]=bdrycond1(t);
    u1[-1]=bdrycond2(t);
    for i in range(1,len(u0)-1):
        u1[i]=l/2*u0[i-1]+(1-l)*u0[i]+l/2*u0[i+1];
    u1[1]+=l/2*bdrycond1(t);
    u1[-2]+=l/2*bdrycond2(t);

    A=[];
    A.append(array([1+l,-l/2]));
    for i in range(3,len(u0)-1):
        A.append(array([-l/2,1+l,-l/2]));
    A.append(array([-l/2,1+l]));
    A=array(A);
   
    u1[2]-=u1[1]*A[1][0]/A[0][0];
    A[1][:2]-=A[0]*A[1][0]/A[0][0];
    for i in range(1,len(u0)-3):
        u1[i+2]-=u1[i+1]*A[i+1][0]/A[i][1];
        A[i+1][:2]-=A[i][1:]*A[i+1][0]/A[i][1];

    for i in arange(len(u0)-3,1,-1):
        u1[i]-=u1[i+1]*A[i-1][2]/A[i][1];
        A[i-1][2]=0;
    u1[1]-=u1[2]*A[0][1]/A[1][1];
    A[0][1]=0;

    u1[1]/=A[0][0];
    for i in range(1,len(u0)-2):
        u1[i+1]/=A[i][1];

    return u1;



#The main function that plots the contours.
def plotter(D,ti,tf,xi,xf,delt,delx,method):
    '''D=diffusion coefficient, ti=initial time, tf=final time,
xi=initial x, xf=final x, delt=time step, delx=space step
method -> Chooses method. Please enter one of the following strings:
"Explicit", "Implicit", "Crank-Nicolson" '''
    l=D*delt/delx**2;
    u=zeros((int((tf-ti)/delt)+1,int((xf-xi)/delx)+1));

    time=[ti]; t=ti;
    space=linspace(xi,xf,int((xf-xi)/delx)+1);
    u[0]=inicond(space);
    
    for i in range(1,int((tf-ti)/delt)+1):
        t+=delt;
        time.append(t);
        if (method=='Explicit'):
            u[i]=explicit(l,u[i-1],t);
        elif (method=='Implicit'):
            u[i]=implicit(l,u[i-1],t);
        elif (method=='Crank-Nicolson'):
            u[i]=cranknicolson(l,u[i-1],t);
        else: return 'Invalid method name';

    time=array(time);
    cp=contourf(space,time,u,levels=100,cmap=cm.jet);
    colorbar(cp);
    xlabel('x');
    ylabel('t');
    title('Contour Plot, %s Method'%(method)); 
    show();
    
    
