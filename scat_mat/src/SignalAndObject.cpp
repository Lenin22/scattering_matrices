#include "stdafx.h"
#include "SignalAndObject.h"
#include <valarray>
#include <complex>

inline double min(double a, double b) { return(a > b ? b : a); };
inline double max(double a, double b) { return(a > b ? a : b); };

MatPoint::MatPoint(double az0, double el0, double R0,
				   double azV0, double elV0, double V0) :
	az(az0), el(el0), R(R0), azV(azV0), elV(elV0), V(V0) {}

double MatPoint::getTau(double t)
{
	vector<double> R_vec = { R*cos(el)*cos(az),
							R*cos(el)*sin(az),
							R*sin(el) };

	vector<double> V_vec = { V*cos(elV)*cos(azV),
							V*cos(elV)*sin(azV),
							V*sin(elV) };
	double RV = 0;
	if (!(R==0.0 || V==0.0)){
		RV = (R_vec[0] * V_vec[0] + R_vec[1] * V_vec[1] + R_vec[2] * V_vec[2]) / (R*V);
	}
	
	double a = pow(LIGHTSPEED, 2) - pow(V, 2);
	double b = -2 * (pow(LIGHTSPEED, 2)*t + RV);
	double c = pow(LIGHTSPEED, 2)*pow(t,2) - pow(R, 2);
	double D = b*b - 4 * a*c;
	double tau = (-b - sqrt(D)) / (2 * a);
	return tau;
}

MatPoint::~MatPoint()
{
}

RLS::RLS(double BW0, double f00, double Ti0, 
	     double Tp0, double Fs0, double Tnp0):
	BW(BW0), f0(f00), Ti(Ti0), Tp(Tp0), Fs(Fs0), Tnp(Tnp0){}

vector<vector<double>> RLS::signal(MatPoint& point)
{
	double alpha = BW/(2*Ti);
	int LT = ceil(Tp*Fs);
	vector<vector<double>> signal(LT, vector<double>(2, 0));
	
	double tns = point.R / LIGHTSPEED;
	double gamma = (point.getTau(Ti) - point.getTau(0)) / Ti;
	double newTi = Ti / gamma;
	
	if ((tns + newTi >= Tnp) && (tns <= Tnp + Tp)) {
		
		double tks = min(tns+newTi, Tnp+Tp);
		tns = Tnp + ceil((max(Tnp, tns)-Tnp)*Fs)/Fs;
		
		vector<double> t;
		double t0 = tns;
		int iter = 0;
		while (t0 < tks - 1 / Fs) {
			t.push_back(t0);
			t0 += 1 / Fs;
			iter++;
		}

		int np = floor((tns - Tnp)*Fs);
				
		vector<double> phase(t.size(), 0);
		for (int i = 0; i < t.size(); i++) {
			phase[i] = 2 * M_PI*((f0 - BW / 2)*point.getTau(t[i])
					   +alpha*pow(point.getTau(t[i]),2));
		}

		for (int i = np; i < np + t.size(); i++) {
			signal[i][0] = cos(phase[i - np]);
			signal[i][1] = sin(phase[i - np]);
		}
	}

	return signal;
}

vector<vector<double>> RLS::fft(vector<vector<double>>& signal)
{
	//fft ----------------------------------------------------
	int N = signal.size();

	fftw_complex *in, *out;
	fftw_plan p;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

	for (int i = 0; i < N; i++) {
		in[i][0] = signal[i][0];
		in[i][1] = signal[i][1];
	}

	p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_execute(p);

	vector<vector<double>> spectrum(N, vector<double>(2, 0));
	for (int i = 0; i < N; i++) {
		spectrum[i][0] = out[i][0];
		spectrum[i][1] = out[i][1];
	}

	fftw_destroy_plan(p);
	fftw_free(in); fftw_free(out);
	return spectrum;
}

vector<vector<double>> RLS::ifft(vector<vector<double>>& spectrum)
{
	//fft ----------------------------------------------------
	int N = spectrum.size();

	fftw_complex *in, *out;
	fftw_plan p;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

	for (int i = 0; i < N; i++) {
		in[i][0] = spectrum[i][0];
		in[i][1] = spectrum[i][1];
	}

	p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

	fftw_execute(p);

	vector<vector<double>> signal(N, vector<double>(2, 0));
	for (int i = 0; i < N; i++) {
		signal[i][0] = out[i][0];
		signal[i][1] = out[i][1];
	}

	fftw_destroy_plan(p);
	fftw_free(in); fftw_free(out);
	return signal;
}

RLS::~RLS()
{
}
