#ifndef VOLUMEOBJECT_H
#define VOLUMEOBJECT_H
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iterator>
using namespace std;

#define _USE_MATH_DEFINES
#ifndef M_PI
    #define M_PI 3.1415f
#endif

class VolumeObject{
public:

    ~VolumeObject();
    virtual void getVisiblePts(double angle) = 0;

    vector<double> complexMult(vector<double>& A, vector<double>& B);

    vector<double> cross(vector<double>& A, vector<double>& B);

    double dot(vector<double>& A, vector<double>& B);

    void matMult(vector<vector<double>>& A,
                 vector<vector<double>>& B,
                 vector<vector<double>>& C);

    void matMult(vector<vector<double>>& A,
                 vector<double>& B,
                 vector<double>& C);

    void rotBody(char axe, vector<double>& pts, vector<double>& out, double angle);
    void rotBody(char axe, vector<vector<double>>& pts, vector<vector<double>>& out, double angle);
	vector<vector<vector<double>>> getRCS(double lambda, vector<double>& angles, char mode);

    void dumpBody(string str);
    void dumpNorm(string str);
    void dumpVectors(string str, vector<vector<double>>& vec);
    void dumpVector(string str, vector<double>& vec);

    vector<vector<double>> body;
    vector<vector<double>> norm;
    vector<vector<vector<double>>> body3d;
    vector<vector<vector<double>>> norm3d;
    vector<bool> visible;

protected:

    double step;
};

class Cone : public VolumeObject{
public:
    Cone(double h, double r, double dh);
    ~Cone(){};
    void getVisiblePts(double angle);//lrb

private:
    double height;
    double radius;
};

class Disk : public VolumeObject{
public:
    Disk(double r, double dh);
    ~Disk(){};
    void getVisiblePts(double angle);//bt
private:
    double radius;
};

class Cylinder : public VolumeObject{
public:
    Cylinder(double h, double radius, double dh);
    ~Cylinder(){};
    void getVisiblePts(double angle);//rbt

private:
    double height;
    double radius;
};

class Rocket : public VolumeObject{
public:
    Rocket(double h1, double h2, double radius, double dh);
    ~Rocket(){};
    void getVisiblePts(double angle); //lco, rco, rcyl, tcyl

private:
    double height1;
    double height2;
    double radius;
};

#endif // VOLUMEOBJECT_H
