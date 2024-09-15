##ASSIGNMENT 6
##MOLECULAR DYNAMICS USING LJ POTENTIAL


#We simulate a fluid where the particles interact with a cut-off Lennard-Jones
#Potential. The system starts out from a 4x4x4 fcc lattice arrangement with
#periodic boundary conditions.
#Values of parameters used:

    #No of steps=2000
    #density (rescaled) = 0.55
    #Initial temperature = T0 = 0.2
    #Time step = dt = 0.002
    #Cutoff r value = rc = 2.5

from numpy import *
from matplotlib.pyplot import *


#This sets the initial position for the particles as an fcc lattice
def initpos(ncells,a):
    pos=zeros((4*ncells**3,3));

    for i in range(4*ncells**3):
        nz=int(i/(2*ncells**2));
        nx=int((i%(2*ncells**2))/ncells);
        ny= i%ncells;
        pos[i][2]=nz*a/2;
        pos[i][0]=nx*a/2;
        pos[i][1]=ny*a + ((nx+nz)%2)*a/2;
        
    return pos;


#Sets initial velocities according to a Gaussian (normal) distribution.
#I also adjust the values of the velocities such that total sum of velocities
#(and hence total momentum) is 0.
def initv(N,T):
    vvals=zeros((N,3));
    while(True):
        vf=zeros(3);
        for i in range(N-1):
            vvals[i]=random.normal(0,sqrt(T),3); #normal distribution with
                                                 #variance as temperature, T
            vf-=vvals[i];

        vvals[-1]=vf;       #adjust velocity so that total momentum is 0.
        if(abs(var(vvals)-T)<0.05): return vvals;




#Calculates the force on a particle.
#I have shifted the force by a constant value so that F is continuous and 0
#at r=rc.
def force(x,y,z,rc):
    r=(x**2+y**2+z**2)**0.5;
    if (r>=rc): return zeros(3);
    
    f= (-24*(1/r**7-1/rc**7) + 48*(1/r**13-1/rc**13))/r;
    return array([f*x,f*y,f*z]);




#Returns the potential energy between 2 particles and r*grad(V). I have
#interpolated the potential such that it is continuous and differentiable
#at r=rc. I have added extra terms of the form -6*r/rc**7 and so on.
def pe(r,rc):
    if (r>=rc): return (0.,0.);
    
    V= 4*(-1/r**6 -6*r/rc**7 + 7/rc**6 + 1/r**12 + 12*r/rc**13 - 13/rc**12);
    pr= -24*(1/r**7-1/rc**7) + 48*(1/r**13-1/rc**13);
    return array([V,r*pr]);




#Main code. 
def simulation(steps,ncells,rho,T,dt,rc):
    '''steps=no of simulation steps;  ncells=no of fcc cells along 1 dimension
rho=density;  T=temperature;   dt=time step;   rc=potential cutoff distance'''
    #no of particles and cell size
    N= 4*ncells**3; a=(4/rho)**(1/3);


    #initializing positions and velocities.
    pos0=initpos(ncells,a);
    vel=initv(N,T);
    pos1=zeros((N,3));
    pos2=zeros((N,3));
    pevals=[]; kevals=[]; tevals=[]; pvals=[];


    #The first time step. Position is calculated using usual Taylor series
    K,U,p=0,0,0;
    for i in range(N):
        f=zeros(3);
        
        for j in range(N):
            if(i==j): continue;

            r=zeros(3);
            #Calculating pairwise force
            for k in range(3):
                if(abs(pos0[i][k]-pos0[j][k])<ncells*a/2):
                    r[k]= pos0[i][k]-pos0[j][k];
                elif(pos0[j][k]>pos0[i][k]):
                    r[k]=ncells*a + pos0[i][k] - pos0[j][k];
                else: r[k]=pos0[i][k]-pos0[j][k]-ncells*a;

            f+=force(r[0],r[1],r[2],rc);

            #Calculating potential energy, kinetic energy and pressure
            if (j>i):
                d=(r[0]**2+r[1]**2+r[2]**2)**0.5;
                U+= pe(d,rc)[0];
                p+= pe(d,rc)[1];
        K+=((vel[i][0])**2 +(vel[i][1])**2 +(vel[i][2])**2)/2;
            
        for k in range(3):
            pos1[i][k]=(pos0[i][k]+ vel[i][k]*dt + 1/2 * f[k]* dt**2)%(ncells*a);


    K/= N; U/= N;
    p=(p/(3*N)+K*2/3)*rho; 
    pevals.append(U);
    kevals.append(K);
    tevals.append(K+U);
    pvals.append(p);


    #2nd step onwards, the Verlet algorithm is used for position.
    #Velocity is calculated using central difference method.
    for n in range(1,steps):
        K,U,p=0,0,0;
        for i in range(N):
            f=zeros(3);
            
            for j in range(N):
                if(i==j): continue;

                r=zeros(3);
                for k in range(3):
                    if(abs(pos1[i][k]-pos1[j][k])<ncells*a/2):
                        r[k]= pos1[i][k]-pos1[j][k];
                    elif(pos1[j][k]>pos1[i][k]):
                        r[k]=ncells*a+pos1[i][k]-pos1[j][k];
                    else: r[k]=pos1[i][k]-pos1[j][k]-ncells*a;

                f+=force(r[0],r[1],r[2],rc);

                if (j>i):
                    d=(r[0]**2+r[1]**2+r[2]**2)**0.5;
                    U+= pe(d,rc)[0];
                    p+= pe(d,rc)[1];

            #Calculating new positions using Verlet algorithm
            for k in range(3):
                if (abs(pos1[i][k]-pos0[i][k])>a):
                    if(pos0[i][k]>pos1[i][k]):
                        pos2[i][k]=(2*pos1[i][k]-pos0[i][k]+ncells*a + f[k]*dt**2)%(ncells*a);
                    else:
                        pos2[i][k]=(2*pos1[i][k]-pos0[i][k]-ncells*a + f[k]*dt**2)%(ncells*a);
                else:
                    pos2[i][k]=(2*pos1[i][k]-pos0[i][k] + f[k]*dt**2)%(ncells*a);

                if (abs(pos0[i][k]-pos2[i][k])>a):
                    if(pos0[i][k]>pos2[i][k]):
                        vel[i][k]=(pos2[i][k]-pos0[i][k]+ncells*a)/(2*dt);
                    else:
                        vel[i][k]=(pos2[i][k]-pos0[i][k]-ncells*a)/(2*dt);
                else:
                    vel[i][k]=(pos2[i][k]-pos0[i][k])/(2*dt);
                K+= (vel[i][k])**2;

                pos0[i][k]=pos1[i][k];
                pos1[i][k]=pos2[i][k];

        K/= (2*N); U/= N;
        p=(p/(3*N**2)+K*2/3)*rho;
        pevals.append(U);
        kevals.append(K);
        tevals.append(K+U);
        pvals.append(p);


    data=array([arange(0,steps,1)*dt,array(kevals),array(pevals),array(tevals),array(pvals)]);

    savetxt('avg_vals.txt',data.T,'%1.3f\t\t%1.4f\t\t%1.4f\t\t%1.4f\t\t%1.4f',header='Averages varying with time steps\n time\t\tKE\t\tPE\t\tTE\t\tPr');




def plotter(ncells,rho,rc):
    t,K,U,E,P=loadtxt('avg_vals.txt',unpack='True');

    plot(500*t,array([K,U,E,P]).T);
    xlabel('time steps');
    ylabel('Quantities');
    title('Thermodynamic Quantities, %dx%dx%d fcc lattice, rc=%.1f, rho=%.2f'%(ncells,ncells,ncells,rc,rho));
    legend(('KE/N','PE/N','Total E/N','Pressure'));
    show();

    avgK= mean(K[len(K)-500:]);
    avgU= mean(U[len(U)-500:]);
    avgP=mean(P[len(P)-500:]);

    print('Average KE per particle =',avgK);
    print('Average temperature =',2/3*avgK);
    print('Average PE per particle =',avgU);
    print('Total E per particle =',avgK+avgU);
    print('Average Pressure =',avgP);

##COMMENTS
    #I have plotted the kinetic energy, potential energy, total energy and
    #the pressure as a function of time.
    #I did not plot the temperature as it is simply proportional to the kinetic
    #energy and would give a similar plot.

    #We see that the quantities gradually stabilize and reach equilibrium.
    
    #The code is not perfect as the energy begins to vary over large periods
    #of time. This could be due to the local error in the forces that
    #accummulates over time. Also, the potential is not completely smooth as
    #my interpolated potential has a discontinuity in the second derivative.

