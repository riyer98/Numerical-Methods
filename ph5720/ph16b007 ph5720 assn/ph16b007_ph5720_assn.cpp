//CODE TO SOLVE FOR DIFFERENTIAL EQUATION OF A FORCED, DAMPED, SIMPLE PENDULUM

//Comments about the four methods are written in a separate file, so I shall not
//be including them here

#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//Euler Predictor-Corrector Method
void eulerpred(double* theta, double* omega, double k, double A, double h, double W) 
{
	int i;
	for (i = 0; i < 30; i++)
	{
		omega[i + 1] = (omega[i] * (1 - k * h / 2 - h * h / 4) - h * theta[i] + h * A / 2 * (cos(W * i * h) + cos(W * (i + 1) * h))) / (1 + k * h / 2 + h * h / 4);
		theta[i + 1] = theta[i] + h / 2 * (omega[i] + omega[i + 1]);
	}
}

//Taylor-Series Method
void taylor(double* theta, double* omega, double k, double A, double h, double W)
{
	int i;
	for (i = 0; i < 30; i++)
	{
		omega[i + 1] = omega[i] + h * (-theta[i] - k * omega[i] + A * cos(W * i * h)) + h * h / 2 * (-omega[i] + k * theta[i] + k * k * omega[i] - k * A * cos(W * i * h) - A * W * sin(W * i * h));
		theta[i + 1] = theta[i] + h * omega[i] + h * h / 2 * (-theta[i] - k * omega[i] + A * cos(W * i * h));
	}
}

//this takes in the arguments u=(theta,omega) and returns p=du/dt
void f(double* u, double* p, double t, double k, double A, double W)
{
	p[0] = u[1];
	p[1] = -u[0] - k * u[1] + A * cos(W * t);
}

//2nd Order Runge-Kutta (Heun's) Method 
void rk2(double* theta, double* omega, double k, double A, double h, double W)
{
	int i; double* k1, * k2, * u, * temp;
	k1 = new double[2];
	k2 = new double[2];
	u = new double[2];
	temp = new double[2];
	for (i = 0; i < 30; i++)
	{
		u[0] = theta[i];
		u[1] = omega[i];
		f(u, k1, i * h, k, A, W);
		temp[0] = u[0] + h * k1[0];
		temp[1] = u[1] + h * k1[1];
		f(temp, k2, (i + 1) * h, k, A, W);
		omega[i + 1] = omega[i] + h / 2 * (k1[1] + k2[1]);
		theta[i + 1] = theta[i] + h / 2 * (k1[0] + k2[0]);
	}
}

//4th Order Runge-Kutta Method
void rk4(double* theta, double* omega, double k, double A, double h, double W)
{
	int i; double* u, * temp, *k1, * k2, * k3, * k4;
	k1 = new double[2];
	k2 = new double[2];
	k3 = new double[2];
	k4 = new double[2];
	u = new double[2];
	temp = new double[2];

	for (i = 0; i < 30; i++)
	{
		u[0] = theta[i];
		u[1] = omega[i];
		f(u, k1, i * h, k, A, W);
		temp[0] = u[0] + h / 2 * k1[0];
		temp[1] = u[1] + h / 2 * k1[1];
		f(temp, k2, i * h + h / 2, k, A, W);
		temp[0] = u[0] + h / 2 * k2[0];
		temp[1] = u[1] + h / 2 * k2[1];
		f(temp, k3, i * h + h / 2, k, A, W);
		temp[0] = u[0] + h * k3[0];
		temp[1] = u[1] + h * k3[1];
		f(temp, k4, (i + 1) * h, k, A, W);
		omega[i + 1] = omega[i] + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
		theta[i + 1] = theta[i] + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
	}
}

//When we run the program, we must enter the file name in which we want to save the data
int main(int argc, char** argv)
{
	double h, k, A, W, * theta, * omega; int n, i;
	
	cout << "Enter the following values." << endl;
	cout << "k = damping, A = driven force amplitude, W = driven force frequency, h = step size" << endl;
	cout << "k = ";
	cin >> k;
	cout << "h = ";
	cin >> h;
	cout << "A = ";
	cin >> A;
	cout << "W = ";
	cin >> W;
	
	cout << "What method of solving to be used? Enter any of the following numbers:" << endl;
	cout << " 1 - Euler Predictor Corrector\n 2 - Taylor Series\n 3 - 2nd Order Runge-Kutta (Heun's Method)\n 4 - 4th Order Runge-Kutta" << endl;
	cin >> n;

	theta = new double[31];
	omega = new double[31];
	//Initial Conditions
	theta[0] = 0.2; omega[0] = 0;

	ofstream data;
	switch (n)
	{
	case(1):
		eulerpred(theta, omega, k, A, h, W);
		data.open(argv[1], ios::trunc);
		data << "#Time\t\tTheta(t)" << endl;
		for (i = 0; i < 31; i++) data << i * h << "\t\t" << theta[i] << endl;
		data.close();
		cout << "data theta(t) saved in " << argv[1] << endl;
		break;

	case(2):
		taylor(theta, omega, k, A, h, W);
		data.open(argv[1], ios::trunc);
		data << "#Time\t\tTheta(t)" << endl;
		for (i = 0; i < 31; i++) data << i * h << "\t\t" << theta[i] << endl;
		data.close();
		cout << "data theta(t) saved in " << argv[1] << endl;
		break;

	case(3):
		rk2(theta, omega, k, A, h, W);
		data.open(argv[1], ios::trunc);
		data << "#Time\t\tTheta(t)" << endl;
		for (i = 0; i < 31; i++) data << i * h << "\t\t" << theta[i] << endl;
		data.close();
		cout << "data theta(t) saved in " << argv[1] << endl;
		break;

	case(4):
		rk4(theta, omega, k, A, h, W);
		data.open(argv[1], ios::trunc);
		data << "#Time\t\tTheta(t)" << endl;
		for (i = 0; i < 31; i++) data << i * h << "\t\t" << theta[i] << endl;
		data.close();
		cout << "data theta(t) saved in " << argv[1] << endl;
		break;

	default: break;
	}

	return 0;
}