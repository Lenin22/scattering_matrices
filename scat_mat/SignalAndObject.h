#pragma once
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES 
#include <cmath>
#include <fftw3.h>
#define LIGHTSPEED 299792458.0
#ifndef M_PI
	#define M_PI 3.1415
#endif // !M_PI

using namespace std;

class MatPoint {
	
public:
	double az, el, R, azV, elV, V;
	MatPoint(double az0, double el0, double R0,
			 double azV0, double elV0, double V0);
	double getTau(double t);
	~MatPoint();
};

class RLS {
private:
	double BW, f0, Ti, Tp, Fs, Tnp;
public:
	RLS(double BW0, double f00, double Ti0,
		double Tp0, double Fs0, double Tnp0);
	vector<vector<double>> signal(MatPoint& point);
	vector<vector<double>> fft(vector<vector<double>>& signal);
	vector<vector<double>> ifft(vector<vector<double>>& signal);
	~RLS();
};