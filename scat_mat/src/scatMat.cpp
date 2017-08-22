#include "stdafx.h"
#include "scatMat.h"
#include <fstream>
#include <string>
#include <iostream>
using namespace std;

void ScatMat::fieldInBand(vector<double>& freqs, vector<double>& angles, VolumeObject& obj, char polar)
{
	ofstream file;
	string str;
	file.open(str);
	if (file.is_open()) {
		for (int i = 0; i < freqs.size(); i++) {
			vector<vector<vector<double>>> rcs = obj.getRCS(LIGHTSPEED / freqs[i], angles, polar);
			for (int j = 0; j < rcs.size(); j++) {
				for (int k = 0; k < 3; k++) {
					file << rcs[j][k][0] << " " << rcs[j][k][1] << "\t";
				}
			}
			file << "\n";
		}
	}
}

void ScatMat::readFromFiles(double lfreq, double rfreq, char polar)
{

}
