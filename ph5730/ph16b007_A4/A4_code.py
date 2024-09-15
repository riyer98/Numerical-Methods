##ASSIGNMENT 4
##FFT OF EXAMPLE SIGNALS

from numpy import *
from matplotlib.pyplot import *
import scipy.signal as sg



#I have made 8 cases for the 8 functions.
#I have plotted the modulus of the Fourier transform for each case.

def ftplot(fn,N):
    '''fn=any integer between 1-8 for the 8 cases.
N=No of points taken in discretization.
Plots the time signal and its Fourier transform for each case.'''
    
    if(fn==1):
        f=5;
        t=linspace(0,1,N+1)[:-1];
        x=sin(2*pi*f*t);
        fx=fft.fft(x,4*N)[:2*N];
        w=fft.fftfreq(4*N,1/N)[:2*N];
        plot(t,x);
        title('Sine Function, freq=%1.2f'%(f));
        xlabel('time');
        ylabel('f(t)');
        show();
        plot(w,abs(fx));
        xlim([0,50]);
        title('Sine Wave FT, no of sampling=%d'%(2*N));
        xlabel('frequency');
        ylabel('FT function');
        show();

    elif(fn==2):
        f=5;
        t=linspace(0,1,N+1)[:-1];
        x=cos(2*pi*f*t);
        fx=fft.fft(x,4*N)[:2*N];
        w=fft.fftfreq(4*N,1/N)[:2*N];
        plot(t,x);
        title('Cosine Function, freq=%1.2f'%(f));
        xlabel('time');
        ylabel('f(t)');
        show();
        plot(w,abs(fx));
        xlim([0,50]);
        title('Cosine Wave FT, no of sampling=%d'%(2*N));
        xlabel('frequency');
        ylabel('FT function');
        show();

    elif(fn==3):
        f=5; phi=pi/3;
        t=linspace(0,1,N+1)[:-1];
        x=cos(2*pi*f*t + phi);
        fx=fft.fft(x,4*N)[:2*N];
        w=fft.fftfreq(4*N,1/N)[:2*N];
        plot(t,x);
        title('Cos Wave with Phase shift, freq=%1.2f, phase shift=%1.2f'%(f,phi));
        xlabel('time');
        ylabel('f(t)');
        show();
        plot(w,abs(fx));
        xlim([0,50]);
        title('Square Wave FT, no of sampling=%d'%(2*N));
        xlabel('frequency');
        ylabel('FT function');
        show();

    elif(fn==4):
        f=5; 
        t=linspace(0,1,N+1)[:-1];
        x=sg.square(2*pi*f*t);
        fx=fft.fft(x,4*N)[:2*N];
        w=fft.fftfreq(4*N,1/N)[:2*N];
        plot(t,x);
        title('Square pulse, freq=%1.2f'%(f));
        xlabel('time');
        ylabel('f(t)');
        show();
        plot(w,abs(fx));
        xlim([0,100]);
        title('Square Wave FT, no of sampling=%d'%(2*N));
        xlabel('frequency');
        ylabel('FT function');
        show();
        

    elif(fn==5):
        d=0.2;
        t=linspace(-.5,.5,N+1)[:-1];
        x=zeros(N);
        for i in range(N):
            if (abs(t[i])<=d/2): x[i]=1;

        fx=fft.fft(x,4*N)[:2*N];
        w=fft.fftfreq(4*N,1/N)[:2*N];
        plot(t,x);
        title('Rectangular Pulse, width=%1.2f'%(d));
        xlabel('time');
        ylabel('f(t)');
        show();
        plot(w,abs(fx));
        xlim([0,100]);
        title('Rectangular Pulse FT, no of sampling=%d'%(2*N));
        xlabel('frequency');
        ylabel('FT function');
        show();

    elif(fn==6):
        sigma=0.1;
        t=linspace(-.5,.5,N+1)[:-1];
        x= 1/(sqrt(2*pi)*sigma) * exp(-t**2/(2*sigma**2));
        fx=fft.fft(x,4*N)[:2*N];
        w=fft.fftfreq(4*N,1/N)[:2*N];
        plot(t,x);
        title('Gaussian Pulse, std deviation=%1.2f'%(sigma));
        xlabel('time');
        ylabel('f(t)');
        show();
        plot(w,abs(fx));
        xlim([0,10]);
        title('Gaussian Pulse FT, no of sampling=%d'%(2*N));
        xlabel('frequency');
        ylabel('FT function');
        show();

    elif(fn==7):
        t=linspace(0,1,N+1)[:-1];
        tau=0.2; A=2;
        x=A*exp(-t/tau);
        fx=fft.fft(x,4*N)[:2*N];
        w=fft.fftfreq(4*N,1/N)[:2*N];
        plot(t,x);
        title('Exponential Decay, tau=%1.2f, amp=%1.2f'%(tau,A));
        xlabel('time');
        ylabel('f(t)');
        show();
        plot(w,abs(fx));
        xlim([0,50]);
        title('Exponential Decay FT, no of sampling=%d'%(2*N));
        xlabel('frequency');
        ylabel('FT function');
        show();

    elif(fn==8):
        t=linspace(0,1,N+1)[:-1];
        x=sg.chirp(t,0,0.5,N/32);
        fx=fft.fft(x,4*N)[:2*N];
        w=fft.fftfreq(4*N,1/N)[:2*N];
        plot(t,x);
        title('Chirp Signal');
        xlabel('time');
        ylabel('f(t)');
        show();
        plot(w,abs(fx));
        xlim([0,100]);
        title('Chirp Signal FT, no of sampling=%d'%(2*N));
        xlabel('frequency');
        ylabel('FT function');
        show();

    else: return 0;





##COMMENTS
    #Since we are taking only a limited interval and not the full real line,
    #and also because we are discretizing the signal, the resulting FT plot
    #we get has some errors/fluctuations in it.


    #Ideally, FT of sine and cosine function must be a delta-function as they
    #have only one frequency component. However, the discretization and taking
    #limited interval results in fluctuating components that go as sin(f)/f or
    #the sinc function. These fluctuations depend on the phase of the wave,
    #as we can see that they are different for sine, cos and the shifted cos.

    
    #The square wave FT has peaks at odd multiples of the frequency: that is,
    #at (2n-1)f for n=1,2,3....
    #The amplitude of these peaks is proportional to 1/(2n-1).


    #The FT of a rectangular pulse is the sinc function sin(2*pi*f)/f. I have 
    #taken the modulus, so the negative parts are plotted as positive.


    #The FT of a Gaussian is also a Gaussian, but the width of the FT is the
    #reciprocal (1/sigma) of the original function. That is, if the wave is
    #spread out in the time domain, the frequency domain has a sharp peak and
    #vice-versa.


    #The FT of an exponential decay e^(-l*t) is a Lorentzian curve 1/(l^2+f^2)
    #where l is the half width at half maximum.


    #The chirp signal has varying frequency between a given range. Hence the
    #FT has non-zero values for a range of f and then goes to 0.
