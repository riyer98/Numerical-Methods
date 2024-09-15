##PROJECT
##GROUND STATE ENERGY OF H2 MOLECULE.

#We use Variational and Path Integral Monte Carlo to arrive at the
#the ground state energy of the H2 molecule.

#It involves defining a trial wavefunction

from numpy import *
from matplotlib.pyplot import *


e0=27.211386;


#trial wf for single electron
def phi0(rl,rr,a):
    return exp(-rl/a)+exp(-rr/a);


#gradient of trial wavefunction
def phi01(pos,rl,rr,s,a):
    posl=pos+array([0,0,s/2]);
    posr=pos-array([0,0,-s/2]);
    return -1/a *(exp(-rl/a)/rl *posl + exp(-rr/a)/rr *posr);


#divergence of trial wf
def phi02(rl,rr,a):
    return phi0(rl,rr,a)/a**2 - 2/a * (exp(-rl/a)/rl+exp(-rr/a)/rr);


#correlation f12
def f(r,beta):
    return exp(r/(2*(1+beta*r)));

#correlation gradient
def f1(pos,r,beta):
    return f(r,beta)/(2*r*(1+beta*r)**2) * pos;

#correlation divergence
def f2(r,beta):
    return f(r,beta)/(1+beta*r)**2 * (1/r + 1/(4*(1+beta*r)**2) -1/(1+beta*r)); 



#The joint trial wavefunction
def phi(r1l,r1r,r2l,r2r,r12,a,beta):
    return phi0(r1l,r1r,a)*phi0(r2l,r2r,a)*f(r12,beta);



#The local variational energy for a given set of positions
def vare(pos1,pos2,beta,a,s):
    r1l=sqrt(sum((pos1+array([0.,0.,s/2]))**2));
    r1r=sqrt(sum((pos1-array([0.,0.,s/2]))**2));

    r2l=sqrt(sum((pos2+array([0.,0.,s/2]))**2));
    r2r=sqrt(sum((pos2-array([0.,0.,s/2]))**2));

    r12=sqrt(sum((pos1-pos2)**2));

    ke1 = -1/2 * (phi02(r1l,r1r,a)/phi0(r1l,r1r,a) + 2*sum(f1(pos1-pos2,r12,beta)*\
           phi01(pos1,r1l,r1r,s,a))/(phi0(r1l,r1r,a)*f(r12,beta)) + f2(r12,beta)/f(r12,beta));

    ke2= -1/2 * (phi02(r2l,r2r,a)/phi0(r2l,r2r,a) + 2*sum(f1(pos2-pos1,r12,beta)*\
           phi01(pos2,r2l,r2r,s,a))/(phi0(r2l,r2r,a)*f(r12,beta)) + f2(r12,beta)/f(r12,beta));

    if (s==0):
        pe= -1/r1l - 1/r1r - 1/r2l - 1/r2r + 1/r12;
    else: pe= -1/r1l - 1/r1r - 1/r2l - 1/r2r + 1/r12 + 1/s;
    return e0*(ke1+ke2+pe);
    



#Samples points in space according to |phi(r)|^2 for calculating averages
#I use the metropolis algorithm to accept/reject points
def sampler(beta,s,a,N):
    posvals=[]; pos=random.normal(0,1,6);

    for i in range(N+100):
        newpos = pos+random.uniform(-1/2,1/2,6);
        r1lold=sqrt(sum((pos[:3]+array([0.,0.,s/2]))**2));
        r1rold=sqrt(sum((pos[:3]-array([0.,0.,s/2]))**2));

        r2lold=sqrt(sum((pos[3:]+array([0.,0.,s/2]))**2));
        r2rold=sqrt(sum((pos[3:]-array([0.,0.,s/2]))**2));

        r12old=sqrt(sum((pos[3:]-pos[:3])**2));


        r1lnew=sqrt(sum((newpos[:3]+array([0.,0.,s/2]))**2));
        r1rnew=sqrt(sum((newpos[:3]-array([0.,0.,s/2]))**2));

        r2lnew=sqrt(sum((newpos[3:]+array([0.,0.,s/2]))**2));
        r2rnew=sqrt(sum((newpos[3:]-array([0.,0.,s/2]))**2));

        r12new=sqrt(sum((newpos[3:]-newpos[:3])**2));

        acc=(phi(r1lnew,r1rnew,r2lnew,r2rnew,r12new,a,beta))**2/(phi(r1lold,r1rold,r2lold,r2rold,r12old,a,beta))**2;

        if (acc>1): pos=newpos;
        elif (acc>random.random()): pos=newpos;

        if (i>=100): posvals.append(pos);

    return posvals;



#The code that implements the variational monte carlo for the trial
#wavefunction
def variational(bmin,bmax,s,N):
    betavals=linspace(bmin,bmax,50); grnden=[];
    a=1;
    for i in range(10):
        a=1/(1+exp(-s/a));


    #I sample points for each beta value to and calculate the average of
    #the variational energy
    for beta in betavals:
        posvals=sampler(beta,s,a,N);
        evals=zeros(N);
        for i in range(N):
            pos1=posvals[i][:3];
            pos2=posvals[i][3:];
            evals[i]=vare(pos1,pos2,beta,a,s);

        grnden.append(mean(evals));

    #Plotting the the energy obtained vs beta
    grnden=array(grnden);
    plot(betavals,grnden);
    xlabel('beta (units of 1/a0)');
    ylabel('Energy (eV)');
    title('Energy vs beta, s=%1.2f, sample points =%d'%(s,N));
    show();

    #Prints the output and saves it in a txt file.
    for i in range(50):
        if (grnden[i]==min(grnden)):
            file=open('grndenergy.txt','a');
            file.write('%1.2f\t\t%1.4f\t\t%1.4f\n'%(s,betavals[i],grnden[i]));
            file.close();
            print('beta for which energy is minimum is ',betavals[i])
            print('Ground state energy is ',grnden[i],' eV');
            return;



#The drift function that governs behaviour of G(r,t)
def D(pos1,pos2,a,beta,s):
    r1l=sqrt(sum((pos1+array([0.,0.,s/2]))**2));
    r1r=sqrt(sum((pos1-array([0.,0.,s/2]))**2));
    r2l=sqrt(sum((pos2+array([0.,0.,s/2]))**2));
    r2r=sqrt(sum((pos1-array([0.,0.,s/2]))**2));
    r12=sqrt(sum((pos1-pos2)**2));
    drift=zeros(6);
    for i in range(3):
        drift[i]=f1(pos1-pos2,r12,beta)[i]/f(r12,beta) + phi01(pos1,r1l,r1r,s,a)[i]/phi0(r1l,r1r,a);
        drift[i+3]=f1(pos2-pos1,r12,beta)[i]/f(r12,beta) + phi01(pos2,r2l,r2r,s,a)[i]/phi0(r2l,r2r,a);

    return drift;




#Path integral monte carlo 
def pathint(beta,dt,s,N,npath):
    #setting a, the wavefunction length scale
    a=1;
    for i in range(10):
        a=1/(1+exp(-s/a));
    posvals=sampler(beta,s,a,N); newpos=zeros((N,6));

    pathen=[];

    for n in range(npath):
        evals=[];
        for i in range(N):
            #Calculating variational energy and drift at each point.
            pos1=posvals[i][:3];
            pos2=posvals[i][3:];
            
            evals.append(vare(pos1,pos2,beta,a,s));
            drift=D(pos1,pos2,a,beta,s);

            #New positions are initialized using the drifts.
            for j in range(6):
                newpos[i][j]=random.normal(posvals[i][j]+drift[j]*dt,sqrt(dt));

        et=mean(evals);
        pathen.append(et);
        
        #New positions are accepted/rejected
        for i in range(N):
            pos1=newpos[i][:3]; pos2=newpos[i][3:];
            er=vare(pos1,pos2,beta,a,s);
            #the weight is exponentially related to the local energy.
            acc=exp(-(er-et)*dt);

            if(acc>1): posvals[i]=newpos[i];
            elif(acc>random.random()): posvals[i]=newpos[i];


    print('Ground state energy is ',pathen[-1],'eV');
        
        
        


##COMMENTS
    #I have calculated ground state energy for proton distance=0 which is a
    #Helium atom. In variational Monte Carlo the value comes out to be -78.3 eV
    #Which is close to the actual value of -79.02 eV.

    #For the H2 molecule, the minimum energy is found near s= 1.4 a0 where a0
    #is the Bohr radius. Numerically the value is 0.74 Angstrom.

    #Ground state energy of H2 molecule using variational MC came out to
    #be -31.75 eV for a separation of s= 1.35 a0 = .7 Angstrom.

    #I have plotted beta, the variational parameter vs the energies. The plot
    #is quite irregular because of the randomness of the ensemble that is
    #generated for each run. However, the final results are reasonably
    #accurate.

    #The path integral MC code is not perfect as I made my own code instead
    #of translating the textbook, but if we take a large number of points the
    #results come out accurate.

    #The binding energy of H2 molecule = Ground state of H2 - Ground state of
    #two H atoms = -31.75 - (-27.2) = -4.55 eV

    
