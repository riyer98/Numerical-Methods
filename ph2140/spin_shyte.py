import numpy as np;
import pylab as plt;

sx=np.array([[0,1],[1,0]]);
sy=np.array([[0,-1j],[1j,0]]);
sz=np.array([[1,0],[0,-1]]);
I=np.eye(2);

def spin1(J,hx,hy,hz):
    H=J*np.kron(sz,sz) + hz*(np.kron(sz,I)+np.kron(I,sz)) + hx*(np.kron(sx,I)+np.kron(I,sx)) + hy*(np.kron(sy,I)+np.kron(I,sy));
    print("eigenvalues",np.linalg.eigh(H)[0]);
    print('eigenvectors',np.linalg.eigh(H)[1]);

def spin2():
    hval=[0.5,1,1.5]; 
    for h in hval:
        theta=[]; evals=[];
        for th in np.linspace(0,np.pi/2,51):
            hx=h*np.sin(th); hz=h*np.cos(th);
            H=np.kron(sz,sz) + hz*(np.kron(sz,I)+np.kron(I,sz)) + hx*(np.kron(sx,I)+np.kron(I,sx));
            evals.append(list(np.linalg.eigvalsh(H)));
            theta.append(th);

        eigv=[theta,evals];
##        np.savetxt("eigv.dat",eigv,fmt="%2.5f");
        plt.plot(theta,evals);
        plt.xlabel("angle");
        plt.ylabel("eigenvalues");

        if (h==0.5):
            plt.title("h=0.5");
            plt.savefig("plotfig_h=0.5.pdf");
        if (h==1):
            plt.title("h=1");
            plt.savefig("plotfig_h=1.pdf");
        if (h==1.5):
            plt.title("h=1.5");
            plt.savefig("plotfig_h=1.5.pdf");
        plt.show();
            
            
