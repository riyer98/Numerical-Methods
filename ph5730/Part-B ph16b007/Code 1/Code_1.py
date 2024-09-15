##CODE 1
#Finding roots and minima of Lennard-Jones Potential.

#The LJ potential is given by V(r) = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)
#Here I work in reduced units V -> V/epsilon    and   r -> r/sigma.


##COMMENTS
#The LJ potential has a root at r=sigma. in reduced units: r=1.
#The force is balanced when F(r)=-V'(r)=0 which gives
#r= 2**(1/6) sigma= 1.12246 sigma.

#The minimum energy required to separate 2 atoms in the bound state is
#simply the negative of the value of the potential at the minimum which is
#E=epsilon.


#Defining the potential here.
def V(r):
    return 4*(1/r**12 - 1/r**6);


# Force F(r)= -V'(r) is defined here.
def F(r):
    return 24*(2/r**13 - 1/r**7);



#Brute force root checking for V(r) 
def Vbruteforce(a,b,N):
    '''Enter the interval [a,b] you want to scan. N is the number of
sub-intervals we will divide [a,b] into. The program will check for
a root in each of the sub-intervals. '''
    h=(b-a)/N; solns=[];
    for i in range(N):
        if(V(a+i*h)*V(a+(i+1)*h)<=0):
            solns.append([a+i*h,a+(i+1)*h]);

    print ('Intervals in the range in which V(r)=0 are',solns);



#Brute force root checking for F(r)
def Fbruteforce(a,b,N):
    '''Enter the interval [a,b] you want to scan. N is the number of
sub-intervals we will divide [a,b] into. The program will check for
a root in each of the sub-intervals. '''
    h=(b-a)/N; solns=[];
    for i in range(N):
        if(F(a+i*h)*F(a+(i+1)*h)<=0):
            solns.append([a+i*h,a+(i+1)*h]);

    print ('Intervals in the range in which F(r)=0 are',solns);



#Bisection method for V(r)
def Vrootbisect(a,b,tolerance):
    '''enter interval [a,b]. Note: a must be < b.'''

    
    if (V(a)*V(b)>0):
        return ('Inconclusive'); #This method works only if
                                                #f(a)*f(b) <= 0.
    #if the root is a or b, we directly return a or b.
    if (V(a)==0):
        print('Value of  r for which V(r)=0 is r0 = ',a,'sigma');
        return;
    if (V(b)==0):
        print('Value of  r for which V(r)=0 is r0 = ',b,'sigma');
        return;

    x=a; y=b;

    while(abs(x-y)>=tolerance):
        c=(x+y)/2;
        #bisection step happens here.
        if (V(c)==0): break;
        if(V(x)*V(c)<0):  y=c;
        else:  x=c;

    print('Value of  r for which V(r)=0 is r0 = ',c,'sigma');



#Bisection method for F(r)
def Frootbisect(a,b,tolerance):
    '''enter interval [a,b]. Note: a must be < b.'''

    
    if (F(a)*F(b)>0):
        return ('Inconclusive'); #This method works only if
                                                #f(a)*f(b) <= 0.
    
    if (F(a)==0):
        print('Value of  r for which V(r)=0 is r0 = ',a,'sigma');
        print('Min energy reqd to separate the atoms is -V(r0) = ',-V(a),'epsilon');
        return;
    if (F(b)==0):
        print('Value of  r for which V(r)=0 is r0 = ',b,'sigma');
        print('Min energy reqd to separate the atoms is -V(r0) = ',-V(b),'epsilon');
        return;

    x=a; y=b;

    while(abs(x-y)>=tolerance):
        c=(x+y)/2;
        #bisection step happens here.
        if (F(c)==0): break;
        if(F(x)*F(c)<0):  y=c;
        else:  x=c;

    print('Value of r for which F(r)=0 is r0=',c,'sigma');
    print('Min energy reqd to separate the atoms is -V(r0) = ',-V(c),'epsilon');




##RESULTS
    #I scanned the interval (0.5,10) to search for multiple roots (in case
    #there are more than one, but in this case we have only 1).

    #After zeroing in on the correct interval, I used the bisection method,
    #with tolerance 10^(-6). It confirmed the analytical results.

