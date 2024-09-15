import pylab as plt;

def totient(n):
    if (n<=1):
        return 0;
    if n==2:
        return 1;
    q=0; p=[2];
    for i in range(2,n-1):
        k=0;
        for j in p:
            if (i%j==0 and n%j==0):
                k=1;
                break;
            
            if (j>i**0.5):
                if ((i in p)==False):
                    p.append(i);
            
            
        if (k==0):
            q+=1;
    return(q+2);

def plot(m):
    r=[];
    for i in range(1,m+1):
        r.append(totient(i));
    plt.plot(range(1,m+1),r,'r');
    plt.show();
