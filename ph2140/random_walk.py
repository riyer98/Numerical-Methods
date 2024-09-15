import turtle as tl;
import numpy as np;
import pylab as plt;
import statistics as st;

def rwalk_gauss(T,N):
    xens=[]; yens=[]; 
    for i in range(N):
        x=0;y=0;
        xt=[];yt=[] ; 
        for j in range(T):
            x1=np.random.standard_normal();
            y1=np.random.standard_normal();
            xt.append(x);
            yt.append(y);
            x+=x1;
            y+=y1;
            
            if (i==0):
                if (x1>=0):
                    tl.left(180/np.pi*np.arctan(y1/x1));
                    tl.forward(30*(x1**2 + y1**2)**0.5);
                    tl.right(180/np.pi*np.arctan(y1/x1));
                else:
                    tl.left(180/np.pi*np.arctan(y1/x1)+180);
                    tl.forward(30*(x1**2 + y1**2)**0.5);
                    tl.right(180/np.pi*np.arctan(y1/x1)+180);
        xens.append(xt);
        yens.append(yt);
        
    xens=np.array(xens).T.tolist();
    yens=np.array(yens).T.tolist();
    rT=((np.array(xens[-1])**2 + np.array(yens[-1])**2)**0.5).tolist();

    meanr=[];
    varofr=[];
    for l in range(T):
        meanr.append((st.mean(xens[l])**2 + st.mean(xens[2])**2)**0.5);
        varofr.append(st.variance(xens[l])+st.variance(yens[l]));

    plt.plot(varofr);
    plt.plot((np.array(varofr))**0.5)
    plt.plot(meanr);
    plt.plot(2*np.linspace(0,T,T+1));
    plt.xlabel('time');
    plt.ylabel("<r^2>,<r> etc");
    plt.title('expectation values vs time');
    plt.show();
    plt.hist(rT,bins=50,normed='True');
    plt.show();
    

def rwalk_cauchy(T,N):
    xens=[]; yens=[]; 
    for i in range(N):
        x=0;y=0;
        xt=[];yt=[] ; 
        for j in range(T):
            x1=np.random.standard_cauchy();
            y1=np.random.standard_cauchy();
            xt.append(x);
            yt.append(y);
            x+=x1;
            y+=y1;
            
            if (i==0):
                if (x1>=0):
                    tl.left(180/np.pi*np.arctan(y1/x1));
                    tl.forward(10*(x1**2 + y1**2)**0.5);
                    tl.right(180/np.pi*np.arctan(y1/x1));
                else:
                    tl.left(180/np.pi*np.arctan(y1/x1)+180);
                    tl.forward(10*(x1**2 + y1**2)**0.5);
                    tl.right(180/np.pi*np.arctan(y1/x1)+180);
        xens.append(xt);
        yens.append(yt);
        
    xens=np.array(xens).T.tolist();
    yens=np.array(yens).T.tolist();
    rT=((np.array(xens[-1])**2 + np.array(yens[-1])**2)**0.5).tolist();

    meanr=[];
    varofr=[];
    for l in range(T):
        meanr.append((st.mean(xens[l])**2 + st.mean(xens[2])**2)**0.5);
        varofr.append(st.variance(xens[l])+st.variance(yens[l]));

    plt.plot(varofr);
##    plt.plot((np.array(varofr))**0.5)
##    plt.plot(meanr);
##    plt.plot(2*np.linspace(0,T,T+1));
    plt.xlabel('time');
    plt.ylabel("<r^2>");
    plt.title('<r^2> vs time');
    plt.show();
##    plt.hist(rT,bins=50,normed='True');
##    plt.show();
