import numpy as np
import pylab as plt
sx=np.array([[0,1],[1,0]])
sy=np.array([[0,-1j],[1j,0]])
sz=np.array([[1,0],[0,-1]])
I=np.identity(2)
    
def hamil_eig(J,hx,hy,hz):
    """
    this program calculates the hamiltonian matrix for input values of J, hx, hy
    and hz. then it displays the eigenvalues and the eigenvectors.
    """
    H=J*np.kron(sz,sz)+hx*(np.kron(I,sx)+np.kron(sx,I))+hy*(np.kron(I,sy)+np.kron(sy,I))+hz*(np.kron(I,sz)+np.kron(sz,I))
    E=np.linalg.eigh(H)
    print('eigenvalues\n',E[0])
    print('eigenvectors\n',E[1])

def H_theta(h):
    """
    this program takes input as h, i.e. sqrt(hx**2 + hz**2), (here hy=0). for
    theta values between 0 and Pi/2, it appends the 4 different eigenvalues to 4
    different lists. it saves the eigenvalue list as a dat file.
    then it plots the eigenvalues vs theta in the same plot and saves it as pdf.
    """
    e1=[]
    e2=[]
    e3=[]
    e4=[]
    theta=[]
    H_eval=[]
    for t in np.linspace(0,(np.pi)/2,200):
        H=np.kron(sz,sz)+h*np.cos(t)*(np.kron(I,sx)+np.kron(sx,I))+h*np.sin(t)*(np.kron(I,sz)+np.kron(sz,I))
        Heig=np.linalg.eigvalsh(H)
        e1.append(Heig[0])
        e2.append(Heig[1])
        e3.append(Heig[2])
        e4.append(Heig[3])
        theta.append(t)
        H_eval.append([t,Heig[0],Heig[1],Heig[2],Heig[3]])
    np.savetxt('eig(theta).dat',H_eval,fmt='%2.5f')
    plt.plot(theta,e1,'b')
    plt.plot(theta,e2,'g')
    plt.plot(theta,e3,'y')
    plt.plot(theta,e4,'r')
    plt.xlabel('theta')
    plt.ylabel('eigenvalues')
    plt.title('eigenvalues vs theta')
    plt.savefig('eval_vs_theta.pdf')
    plt.show()

def eigintensity(theta):
    B=[]
    for h in np.linspace(0.5,1.5):
        H=np.kron(sz,sz)+h*np.sin(theta)*(np.kron(I,sx)+np.kron(sx,I))+h*np.cos(theta)*(np.kron(I,sz)+np.kron(sz,I))
        A=np.linalg.eigh(H)
        D=[]
        for i in A[1]:
            D.append((np.linalg.norm(i))**2)
        B.append(D)

    plt.plot(np.linspace(0.5,1.5),B)
    plt.xlabel('h')
    plt.ylabel('eigenvector intensities')
    plt.title('eigenvector intensities vs h for a particular theta')
    plt.savefig('eigintensity.pdf')
    plt.show()


    
