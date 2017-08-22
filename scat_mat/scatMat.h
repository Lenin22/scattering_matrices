#pragma once
#include "volumeobject.h"
#include "SignalAndObject.h"

class ScatMat {
public:
	void fieldInBand(vector<double>& freqs, vector<double>& angles, VolumeObject& obj, char polar);
	void readFromFiles(double lfreq, double rfreq, char polar);
};