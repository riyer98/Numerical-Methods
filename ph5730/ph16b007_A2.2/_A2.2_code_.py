##ASSIGNMENT 2.2
##SOLVING HARMONIC OSCILLATOR AS COUPLED DEs USING RK4
##AND MIDPOINT METHOD

##I have added damping factor as well to generalise it.

import numpy as np
import matplotlib.pyplot as plt

#f is a vector function that returns the DE in coupled
#vector form.
def f(x,w,y):
    '''x=[x1,x2], a vector where x1=x(t) and x2=dx/dt=v.
x must be inputed as a list as shown above.
w=Angular frequency of the oscillator
y=damping parameter'''
    return np.array([x[1],-(w**2)*x[0]-y*x[1]]);


#Midpoint Method
def midpt(x,w,y,h):
    '''x=[x1,x2], a vector where x1=x(t) and x2=dx/dt=v.
x must be inputed as a list as shown above.
w=Angular frequency of the oscillator
y=damping parameter
h=step size
Uses the midpoint method to return the vector x=[x1,x2] in
the next step.'''
    return x+h*f(x+h/2*f(x,w,y),w,y);


#RK4 Method
def rk4(x,w,y,h):
    '''x=[x1,x2], a vector where x1=x(t) and x2=dx/dt=v.
x must be inputed as a list as shown above.
w=Angular frequency of the oscillator
y=damping parameter
h=step size
Uses the RK4 method to return the vector x=[x1,x2] in
the next step.'''
    k0=f(x,w,y);
    k1=f(x+h/2*k0,w,y);
    k2=f(x+h/2*k1,w,y);
    k3=f(x+h*k2,w,y);
    return x+h/6*(k0+2*k1+2*k2+k3);


#The Main Program. Use this function to solve and plot.
def solver(w,y,x0,v0,t0,tf,h):
    '''w=Angluar frequency = sqrt(k/m) => (I have reduced 2 parameters to 1 parameter.)
y=damping parameter, x0=x(t0), v0=v(t0), t0=initial time, tf=final time, h=step size
Solves the SHO using midpoint and RK4 methods for the given
initial condition and plots them along with the analytic solution.'''
    xvecmid=[[x0,v0]]; xvecrk=[[x0,v0]]; t=[t0];

    while(t[-1]<tf-h/2):
        xvecmid.append(midpt(xvecmid[-1],w,y,h));
        xvecrk.append(rk4(xvecrk[-1],w,y,h));
        t.append(t[-1]+h);

    tanalytic=np.linspace(t0,tf,int((tf-t0)*1000)+1);
    xanalytic=[];

    if (w**2>(y/2)**2):
        w1=np.sqrt(w**2 - (y/2)**2);
        for ta in tanalytic:
            xanalytic.append(np.exp(-y/2*(ta-t0))*(x0*np.cos(w1*(ta-t0))+(v0+x0*y/2)/w1*np.sin(w1*(ta-t0))));

    elif (w**2==(y/2)**2):
        for ta in tanalytic:
            xanalytic.append(np.exp(-w*(ta-t0))*(x0+(v0+w*x0)*(ta-t0)));

    else:
        y1=np.sqrt((y/2)**2-w**2);
        A=(x0*(1+y/(2*y1))+v0/y1)/2;
        B=(x0*(1-y/(2*y1))-v0/y1)/2;
        for ta in tanalytic:
            xanalytic.append(A*np.exp((-y/2+y1)*(ta-t0))+B*np.exp((-y/2-y1)*(ta-t0)));

    np.savetxt('x(t),w=%1.1f,y=%1.1f,h=%1.1f.dat'%(w,y,h),np.array([t,np.array(xvecmid).T[0],np.array(xvecrk).T[0]]).T,'%.2f\t%.6f\t%.6f',header='t\tx(t) Midpt\t x(t)RK4');
    np.savetxt('x(t)_analytic,w=%1.1f,y=%1.1f.dat'%(w,y),np.array([tanalytic,xanalytic]).T,'%1.3f\t%1.6f',header='Analytic Solution\nt\tx(t)'); 
    plt.plot(tanalytic,xanalytic);
    plt.plot(t,np.array(xvecmid).T[0]);
    plt.plot(t,np.array(xvecrk).T[0]);
    plt.legend(('Analytic Solution','Midpoint Method','RK4'))
    plt.xlabel('t');
    plt.ylabel('x(t)');
    plt.title('x(t) vs t, x(0)=%1.3f, v(0)=%1.3f, w=%1.3f, y=%1.3f, h=%1.3f'%(x0,v0,w,y,h));
    plt.show();
    



##COMMENTS
    #It is evident that the accuracy of the numerical method decreases
    #with increasing step size.

    #In general, the RK4 method is much more accurate than midpoint method
    #and gives values close to the actual value. It is often seen that the
    #midpoint method's troughs and peaks of the graph do not match
    #with the analytic solution but are closer together (higher frequency
    #than the actual frequency).

    #For the undamped case, the midpoint method's amplitude grows with time.
    #This effect is faster for larger h values. This shows that our
    #approximation is eventually unstable and will diverge. 

    #For the underdamped case, we can also see that for higher values of h,
    #the damping is swamped by the error and the amplitude does not get damped
    #but stays almost constant.

    #For the overdamped case, if the step size is large, the midpoint method
    #is highly unstable and diverges exponentially to infinity.

    #RK4 method, on the other hand, is much more stable, and for the given h
    #values, does not diverge. It may, however, diverge if we take very large
    #step size like h=1 or h=2.
