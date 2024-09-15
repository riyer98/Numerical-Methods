##import numpy as np;

def primef(n):
    if (n<=1):
        print("fuck off");
        return;
    if(n==2):
        print("prime");
        return;

    p=[2];
    for i in range(2,int(n/2)+1):
        for j in p:
            if (i%j==0):
                break;
            if (j>i**0.5):
                p.append(i);
                break;

    m=n;

    for k in p:
        while(m%k==0):
            print(k);
            m=m/k;
        if (m==1):
            break;

    if(m==n):
        print("prime");
        return;
    

def prime(n):
    if (n<=1):
        print("fuck off");
        return;

    q=[]; 

    for k in p:
        if(n%k==0):
            q.append(k);
        if(k>=n):
            break;

    return q;


def piest():
    a=0;
    for j in range(1000):
        n1=np.random.randint(1,1000001);
        n2=np.random.randint(1,1000001);

        if (n1==1 or n2==1):
            a+=1;
            continue;

        p1=prime(n1);
        p2=prime(n2);
        k=0;

        for i in p1:
            if ((i in p2)==True):
                k=1;
                break;
        if(k==0):
            a+=1;

    return ((6000/a)**0.5)
