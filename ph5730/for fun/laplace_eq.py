import numpy as np;
import sys;
import matplotlib.pyplot as plt;

def xbdrymin(y):
    ##return 10*np.sin(np.pi*y/(y[-1]-y[0]));
    return np.sin(2*np.pi*(y-y[0])/(y[-1]-y[0]));


def ybdrymin(x):
    return -2*(x-x[0])*(x[-1]-x[0]);

def ybdrymax(x):
    return 2*(x-x[0])*(x[-1]-x[0]);

def xbdrymax(y):
    return 2*(2*y-(y[0]+y[-1]))/(y[-1]-y[0]);



def sdfj(xmin,xmax,ymin,ymax,h):
    Nx= int((xmax-xmin)/h); Ny= int((ymax-ymin)/h); 
    u=np.zeros((Ny+1,Nx+1));
    x=np.arange(xmin,xmax+h,h); y=np.arange(ymin,ymax+h,h);

    u.T[0]=xbdrymin(y);
    u.T[-1]=xbdrymax(y);
    u[0]=ybdrymin(x);
    u[-1]=ybdrymax(x);
    
    u[1,1:-1] += u[0,1:-1];
    u[-2,1:-1]+= u[-1,1:-1];
    u.T[1,1:-1]+=u.T[0,1:-1];
    u.T[-2,1:-1]+=u.T[-1,1:-1];

    matrix=np.zeros(((Nx-1)*(Ny-1),2*Nx-1));
    for i in range(len(matrix)):
        matrix[i][0]=-1; matrix[i][-1]=-1; 
        if(i%(Nx-1)!=0): matrix[i][Nx-2]=-1; 
        if(i%(Nx-1)!= Nx-2): matrix[i][Nx]=-1;
        matrix[i][Nx-1]=4;

    for i in range(len(matrix)):
        u[int(i/(Nx-1))+1,i%(Nx-1)+1] /= matrix[i,Nx-1];
        matrix[i,Nx-1:] /= matrix[i][Nx-1];

        for j in range(1,Nx):
            if (i+j>=len(matrix)): break;
            if (matrix[i+j,Nx-1-j]==0): continue;
            else: 
                u[int((i+j)/(Nx-1))+1,(i+j)%(Nx-1)+1] -= matrix[i+j,Nx-1-j] * u[int(i/(Nx-1))+1,i%(Nx-1)+1];
                matrix[i+j,Nx-1-j:-j] -= matrix[i+j,Nx-1-j] * matrix[i,Nx-1:];
    
    for i in range(len(matrix)-1,0,-1):
        for j in range(1,Nx):
            if (i-j<0): break;
            if (matrix[i-j,Nx-1+j]==0): continue;
            else: 
                u[int((i-j)/(Nx-1))+1,(i-j)%(Nx-1)+1] -= matrix[i-j,Nx-1+j] * u[int(i/(Nx-1))+1,i%(Nx-1)+1];
    
    solplot=plt.contourf(x,y,u,levels=50,cmap='jet');
    plt.colorbar(solplot);
    plt.contour(x,y,u,levels=50);
    plt.show();


if (__name__=='__main__'):
    if (len(sys.argv)!=6):
        print("Format: xmin xmax ymin ymax h");
        sys.exit(1);
    xmin = float(sys.argv[1]);
    xmax = float(sys.argv[2]);
    ymin = float(sys.argv[3]);
    ymax = float(sys.argv[4]);
    h = float(sys.argv[5]);
    
    sdfj(xmin,xmax,ymin,ymax,h);

    
