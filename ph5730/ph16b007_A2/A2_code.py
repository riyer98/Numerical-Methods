##ASSIGNMENT 2 PART 1
##METHODS FOR SOLVING 1st ORDER ODEs.

import numpy as np;
import matplotlib.pyplot as plt;

#I am using 3 different functions to solve.
#the DE is y'=dy/dt=f(t,y).
def f1(t,y): return (-t*y);
def f2(t,y): return (y+t);
def f3(t,y): return (y-y**2);


#These programs choose 1 of the above 3 ODEs and solve them. 
#The plotter function saves and plots the data simultaneously
##for all 4 methods. 

#The codes for the methods are given below. Note: They only return
#the array of solution values and not the plots. The plots are done
#by the plotter function only.

#Euler Method
def euler(fn_no,h,y0,t0,tmax):
    '''fn_no=enter either 1,2, or 3 only.
    It chooses which DE/f(t,y) to solve.
    h=Step size
    y0=initial condition y(t0)
    t0=initial time
    tmax=final time
    Note: tmax should be > t0.'''
    y=[y0]; t=t0;

    if (fn_no==1):
        while(t<tmax):
            y.append(y[-1]+f1(t,y[-1])*h);
            t+=h;
    if (fn_no==2):
        while(t<tmax):
            y.append(y[-1]+f2(t,y[-1])*h);
            t+=h;
    if (fn_no==3):
        while(t<tmax):
            y.append(y[-1]+f3(t,y[-1])*h);
            t+=h;
    return y;


#Midpoint Method
def midpt(fn_no,h,y0,t0,tmax):
    '''fn_no=enter either 1,2, or 3 only.
    It chooses which DE/f(t,y) to solve.
    h=Step size
    y0=initial condition y(t0)
    t0=initial time
    tmax=final time
    Note: tmax should be > t0.'''
    y=[y0]; t=t0;

    if (fn_no==1):
        while(t<tmax):
            y.append(y[-1]+h*f1(t+h/2,y[-1]+f1(t,y[-1])*h/2));
            t+=h;
    if (fn_no==2):
        while(t<tmax):
            y.append(y[-1]+h*f2(t+h/2,y[-1]+f2(t,y[-1])*h/2));
            t+=h;
    if (fn_no==3):
        while(t<tmax):
            y.append(y[-1]+h*f3(t+h/2,y[-1]+f3(t,y[-1])*h/2));
            t+=h;
    return y;


#Average of Slopes Method
def avg(fn_no,h,y0,t0,tmax):
    '''fn_no=enter either 1,2, or 3 only.
    It chooses which DE/f(t,y) to solve.
    h=Step size
    y0=initial condition y(t0)
    t0=initial time
    tmax=final time
    Note: tmax should be > t0.'''
    y=[y0]; t=t0;

    if (fn_no==1):
        while(t<tmax):
            y.append(y[-1]+h*(f1(t,y[-1])+f1(t+h,y[-1]+h*f1(t,y[-1])))/2);
            t+=h;
    if (fn_no==2):
        while(t<tmax):
            y.append(y[-1]+h*(f2(t,y[-1])+f2(t+h,y[-1]+h*f2(t,y[-1])))/2);
            t+=h;
    if (fn_no==3):
        while(t<tmax):
            y.append(y[-1]+h*(f3(t,y[-1])+f3(t+h,y[-1]+h*f3(t,y[-1])))/2);
            t+=h;
    return y;


#4th Order Runge-Kutta (RK4) Method
def rk4(fn_no,h,y0,t0,tmax):
    '''fn_no=enter either 1,2, or 3 only.
    It chooses which DE/f(t,y) to solve.
    h=Step size
    y0=initial condition y(t0)
    t0=initial time
    tmax=final time
    Note: tmax should be > t0.'''
    y=[y0]; t=t0;

    if (fn_no==1):
        while(t<tmax):
            k0=f1(t,y[-1]);
            k1=f1(t+h/2,y[-1]+h*k0/2);
            k2=f1(t+h/2,y[-1]+h*k1/2);
            k3=f1(t+h,y[-1]+h*k2);
            y.append(y[-1]+h*(k0+2*k1+2*k2+k3)/6);
            t+=h;

    if (fn_no==2):
        while(t<tmax):
            k0=f2(t,y[-1]);
            k1=f2(t+h/2,y[-1]+h*k0/2);
            k2=f2(t+h/2,y[-1]+h*k1/2);
            k3=f2(t+h,y[-1]+h*k2);
            y.append(y[-1]+h*(k0+2*k1+2*k2+k3)/6);
            t+=h;

    if (fn_no==3):
        while(t<tmax):
            k0=f3(t,y[-1]);
            k1=f3(t+h/2,y[-1]+h*k0/2);
            k2=f3(t+h/2,y[-1]+h*k1/2);
            k3=f3(t+h,y[-1]+h*k2);
            y.append(y[-1]+h*(k0+2*k1+2*k2+k3)/6);
            t+=h;
        
    return y;


#The Plotter
def plotter(fn_no,h,y0,t0,tmax):
    '''fn_no=enter either 1,2, or 3 only.
    It chooses which DE/f(t,y) to solve.
    h=Step size
    y0=initial condition y(t0)
    t0=initial time
    tmax=final time
    Note: tmax should be > t0.'''
    t=[t0]; yanalytic=[];

    #The analytic solution array is created here.
    tanalytic=np.linspace(t0,tmax,int((tmax-t0)*1000)+1);
    while (t[-1]<tmax): t.append(t[-1]+h);
    if (fn_no==1):
        for ta in tanalytic:
            yanalytic.append(y0*np.exp(-(ta**2-t0**2)/2));
    if (fn_no==2):
        for ta in tanalytic:
            yanalytic.append(-ta-1+(t0+y0+1)*np.exp(ta-t0));
    if(fn_no==3):
        for ta in tanalytic:
            yanalytic.append(y0/(y0+(1-y0)*np.exp(t0-ta)));

    yeu=euler(fn_no,h,y0,t0,tmax);
    ymid=midpt(fn_no,h,y0,t0,tmax);
    yavg=avg(fn_no,h,y0,t0,tmax);
    yrk=rk4(fn_no,h,y0,t0,tmax);

    if(t[-1]>tmax+h/2):
        t=t[:-1];
        yeu=yeu[:-1];
        ymid=ymid[:-1];
        yavg=yavg[:-1];
        yrk=yrk[:-1];

    np.savetxt('f%d_h=%1.2f.dat'%(fn_no,h),np.array([t,yeu,ymid,yavg,yrk]).T,'%1.2f\t%1.6f\t%1.6f\t%1.6f\t%1.6f');
    np.savetxt('f%d_analytic.dat'%(fn_no),np.array([tanalytic,yanalytic]).T,'%1.4f\t%1.6f');
    plt.plot(t,yeu);
    plt.plot(t,ymid);
    plt.plot(t,yavg);
    plt.plot(t,yrk);
    plt.plot(tanalytic,yanalytic);
    plt.xlabel('t');
    plt.ylabel('y(t)');
    plt.legend(('Euler','Midpoint','Average','RK4','Analytic Soln'));
    plt.title('Plots for DE %d, y(%1.2f)=%1.2f, h=%1.3f'%(fn_no,t0,y0,h));
    plt.show();



##COMMENTS
    #I have plotted the solution to each DE for 2 different h values, 0.2
    #and 0.4. The initial conditions used for each DE are as follows:
    #1) y1(0)=1      2)y2(0)=0        3)y3(0)=1/2
    
    
    #It is obvious that with decreasing step size h, the accuracy of the
    #methods increase.

    #Euler method is the most inaccurate as it only takes the linear
    #order term. It deviates quite significantly from the actual solution.

    #The midpoint and average methods have O(h^2) accuracy and both more
    #or less have the same accuracy.

    #RK4 has O(h^4) accuracy and is accurate even for large step sizes.
    #Even though the interpolation between the points  for RK4 deviates,
    #The points themselves are very close to the actual solution.
