##CODE 2
##SOLVING GIVEN DE WITH DIRECT FINITE DIFFERENCE METHOD

#The DE given is  y'' +3y'-5y =7x  with boundary conditions
#y(0)=-20, y(1)=100

#We use central difference for the x derivatives. We discretize the
#unit intverval with step dt. If the index for the steps be i.
#Then the equation is

#(y(i+1)-2y(i)+y(i-1))/(dt^2) + 3(y(i+1)-y(i-1))/(2dt) - 5 y(i) = 7x(i)
#Multiplying by dt^2 and rearranging gives

#(1-3dt/2)y(i-1) - (2+5dt^2)y(i) + (1+3dt/2)y(i+1) = 7dt^2 x(i)

#The boundary terms are y(xi=0)=-20 and y(xi=1)=100

#For the middle terms, this can be written as a matrix equation A.y=b where

#b=7dt^2 x

#A is a tridiagonal matrix with diagonal elements -(2+5dt^2) and the other 2
#being the off diagonal elements. All other elements are 0.

#A=( -(2+5dt^2)     (1+3dt/2)
#    (1-3dt/2)      -(2+5dt^2)      (1+3dt/2)
#                   (1-3dt/2)      -(2+5dt^2)      (1+3dt/2)
#                                   (1-3dt/2)      -(2+5dt^2)      (1+3dt/2) 
#       ...             ...             ...             ...          ...
#                                                                 -(2+5dt^2)      (1+3dt/2) )

#We can perform row operations on the matrix to reduce it to the identity
#matrix/upper triangular matrix. Then we can directly use it to determine
#y(i).

#Central difference derivative gives accuracy of O(dt^2) at each step. If
#we integrate the equations, the global error will be O(dt). Hence, if we
#want accuracy to 4 significant digits, we need dt<= 10^(-4).

#I have solved the equations for dt=10^(-4) and 10^(-5). I have plotted it
#with the analytic expressions as well.


from numpy import *
from matplotlib.pyplot import *


#The code that performs row ops on the tridiagonal matrix A.  
def tridiag(A,x,y,b):
    #These are the row ops: R(n+1) = R(n+1) - R(n) * A[n+1][n]/A[n][n];
    #This is to bring the off diagonal elements to 0 and solve the matrix.
    b[1]-=A[1][0]*b[0]/A[0][0];
    A[1][:2]-=A[0]*A[1][0]/A[0][0];

    for i in range(1,len(b)-1):
        b[i+1]-=b[i]*A[i+1][0]/A[i][1];
        A[i+1][:2]-=A[i][1:]*A[i+1][0]/A[i][1];


    #Now we repeat the procedure backwards for the second off-diagonal.
    for i in arange(len(b)-1,1,-1):
        b[i-1]-=b[i]*A[i-1][2]/A[i][1];
        A[i-1][2]=0;

    b[0]-=b[1]*A[0][1]/A[1][1];
    A[0][1]=0;

    for i in range(1,len(b)):
        y[i+1]=b[i]/A[i][1];
    y[1]=b[0]/A[0][0];
        



#The main code that stores data and plots the functions.
def diffsolve(xmin,xmax,yi,yf,dt):
    '''xmin=initial x, xmax=final x
       yi=initial y, yf=final y
       dt=time step.'''
    N=int(1./dt)+1;

    #setting x values
    x=linspace(xmin,xmax,N);
    y=zeros(N);

    #boundary values
    y[0]=yi; y[-1]=yf;

    #Initialize b
    b=x[1:N-1]*7*dt**2; 

    b[0]-=(1-3*dt/2)*yi;
    b[-1]-=(1+3*dt/2)*yf;

    #Initializing A
    A=[array([-(2+5*dt**2),1+3*dt/2])];

    for i in range(1,N-3):
        A.append(array([1-3*dt/2,-(2+5*dt**2),1+3*dt/2]));
    A.append(array([1-3*dt/2,-(2+5*dt**2)]));

    

    A=array(A);
    #Call the function to solve for y
    tridiag(A,x,y,b);

    #Analytic expression for the solution.
    l=0.5*(3+sqrt(29));
    c1=(102.24*exp(l)+19.16)/(exp(sqrt(29))-1);
    c2=-19.16-c1;
    ya=-21/25 -7/5 * x + exp(-l*x) * (c1*exp(x*sqrt(29))+c2);
    
    data=array([x,y,ya]).T
    savetxt('y(x),dt=%f.txt'%(dt),data,fmt='%.5f\t\t%.5f\t\t%.5f',header='y(x) values\n x\t\ty_numerical\t\ty_analytic');

    plot(x,ya,x,y);
    xlabel('x');
    ylabel('y(x)');
    title('y(x),dt=%f,'%(dt));
    legend(('analytical','numerical'));
    show();



##RESULTS
    #The analytic and numerical graphs overlap perfectly. We can see in the
    #txt file that the numbers match at least up to 4 significant digits.
