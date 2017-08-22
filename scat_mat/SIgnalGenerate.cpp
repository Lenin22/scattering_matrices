// SIgnalGenerate.cpp: определяет точку входа для консольного приложения.
//
#include "stdafx.h"
#include "SignalAndObject.h"
#include "volumeobject.h"
#include <fstream>


void dumpVector(string str, vector<vector<double>>& vec){

	ofstream file;
	file.open(str);
	if (file.is_open()){
		for (int i = 0; i < vec.size(); i++){
			file << vec[i][0] << "\t" << vec[i][1] << "\n";
		}
	}
	else{
		std::cout << str.c_str() << " not opened";
	}
	file.close();
}

int main()
{
	double BW = 1e6, f0 = 1e9, Ti = 0.5e-3, Tp = Ti, Fs = 4 * BW, Tnp = 2e-3;

	double az = 0, el = M_PI / 2, R = LIGHTSPEED*(Tnp+Ti/8),
		  azV = 0, elV = 0, V = 5000;

	MatPoint pt(az, el, R, azV, elV, V);
	RLS rls(BW, f0, Ti, Tp, Fs, Tnp);

	vector<vector<double>> signal = rls.signal(pt);
	vector<vector<double>> spectrum = rls.fft(signal);
	spectrum = rls.ifft(spectrum);

	dumpVector("signal.txt", signal);
	dumpVector("spectrum.txt", spectrum);

	return 0;
}
