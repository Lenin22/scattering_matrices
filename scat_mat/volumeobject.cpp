#include "stdafx.h"
#include "volumeobject.h"

VolumeObject::~VolumeObject(){
    body.clear();
    norm.clear();
}

vector<double> VolumeObject::complexMult(vector<double>& A, vector<double>& B){
    vector<double> res(2,0);
    res[0] = A[0]*B[0]-A[1]*B[1];
    res[1] = A[0]*B[1]+A[1]*B[0];
    return res;
}

vector<double> VolumeObject::cross(vector<double>& A, vector<double>& B){
    vector<double> res(3,0);
    res[0] = A[1]*B[2]-A[2]*B[1];
    res[1] = A[2]*B[0]-A[0]*B[2];
    res[2] = A[0]*B[1]-A[1]*B[0];
    return res;
}

double VolumeObject::dot(vector<double>& A, vector<double>& B){
    double res = 0.0;
    for(int i = 0; i < A.size(); i++){
        res+=A[i]*B[i];
    }
    return res;
}

void VolumeObject::dumpBody(string str){

    ofstream file;
    file.open(str);
    if(file.is_open()){
        for(int i = 0; i < body3d.size(); i++){
            for(int j = 0; j < body3d[i].size(); j++){
                for(int k = 0; k < 3; k++){
                    file << body3d[i][j][k] << "\t";
                }
                file << "\n";
            }
        }
    }else{
        std::cout << str.c_str() << " not opened";
    }
    file.close();
}

void VolumeObject::dumpNorm(string str){

    ofstream file;
    file.open(str);
    if(file.is_open()){
        for(int i = 0; i < norm3d.size(); i++){
            for(int j = 0; j < norm3d[i].size(); j++){
                for(int k = 0; k < 3; k++){
                    file << norm3d[i][j][k] << "\t";
                }
                file << "\n";
            }
        }
    }else{
        cout << str.c_str() << " not opened";
    }
    file.close();
}

void VolumeObject::dumpVectors(string str, vector<vector<double>>& vec){

    ofstream file;
    file.open(str);
    if(file.is_open()){
        for(int i = 0; i < vec.size(); i++){
            for(int j = 0; j < vec[i].size(); j++){

                file << vec[i][j] << "\t";
            }
            file << "\n";
        }
    }else{
        std::cout << str.c_str() << " not opened";
    }
    file.close();
}

void VolumeObject::dumpVector(string str, vector<double>& vec){

    ofstream file;
    file.open(str);
    if(file.is_open()){
        for(int i = 0; i < vec.size(); i++){
            file << vec[i] << "\n";
        }
    }else{
        std::cout << str.c_str() << " not opened";
    }
    file.close();
}

void VolumeObject::matMult(vector<vector<double>>& A,
                           vector<vector<double>>& B,
                           vector<vector<double>>& C){

    for(int i = 0; i < A.size(); i++){
        for(int j = 0; j < B.at(0).size(); j++){
            C[i][j] = 0;
            for(int k = 0; k < A.at(0).size(); ++k){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void VolumeObject::matMult(vector<vector<double>>& A,
                           vector<double>& B,
                           vector<double>& C){

    for(int i = 0; i < A.size(); i++){
        C[i] = 0;
        for(int j = 0; j < B.size(); j++){
            C[i] += A[i][j] * B[j];
        }
    }
}

void VolumeObject::rotBody(char axe, vector<double>& pts,
                           vector<double>& out, double angle){

    vector<vector<double>> rot;
    switch (axe) {
    case 'x':
        rot = {{1,    0     ,    0      },
               {0,cos(angle),-sin(angle)},
               {0,sin(angle),cos(angle)}};
        break;
    case 'z':
        rot = {{cos(angle),  -sin(angle), 0},
               {sin(angle), cos(angle),  0},
               {    0,         0,        1}};
        break;
    case 'y':
        rot = {{cos(angle), 0,sin(angle)},
               {    0,      1,    0     },
               {-sin(angle),0,cos(angle)}};
        break;
    default:
        break;
    }

    matMult(rot, pts, out);
}

void VolumeObject::rotBody(char axe, vector<vector<double>>& pts,
                           vector<vector<double>>& out, double angle){

    vector<vector<double>> rot;
    switch (axe) {
    case 'x':
        rot = {{1,    0     ,    0      },
               {0,cos(angle),sin(angle)},
               {0,-sin(angle),cos(angle)}};
        break;
    case 'z':
        rot = {{cos(angle),  sin(angle), 0},
               {-sin(angle), cos(angle),  0},
               {    0,         0,        1}};
        break;
    case 'y':
        rot = {{cos(angle), 0,-sin(angle)},
               {    0,      1,    0     },
               {sin(angle),0,cos(angle)}};
        break;
    default:
        break;
    }


    matMult(pts, rot, out);
}

Disk::Disk(double r, double dh){

    step = dh;
    radius = r;
    int N = ceil(2*radius/step);
    double X, Y;
    //init
    body3d.push_back(vector<vector<double>>());
    body3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());

    for(int k = 0; k < N; k++){
        for(int l = 0; l < N; l++){
            X = -radius+(double)k*step;
            Y = -radius+(double)l*step;
            if(X*X+Y*Y<radius*radius){
                body3d[0].push_back({X,Y,0.0});
                norm3d[0].push_back({0.0,0.0,-1.0});
            }
        }
    }
    copy(body3d[0].begin(), body3d[0].end(), back_inserter(body3d[1]));
    copy(norm3d[0].begin(), norm3d[0].end(), back_inserter(norm3d[1]));

    for(int i = 0; i < norm3d[0].size(); i++){
        norm3d[0][i][2] = 1.0;
    }

    double nPts = body3d[0].size();
    double square = M_PI*radius*radius/nPts;

    for(int i = 0; i < 2; i++){
        for(int j = 0; j < body3d[0].size(); j++){
            for(int k = 0; k < 3; k++){
                norm3d[i][j][k] *= square;
            }
        }
    }
}

Cone::Cone(double h, double r, double dheight){

    height = h;
    radius = r;
    step = dheight;
    double alpha = atan(radius/height);
    double dh = step*cos(alpha);
    int nh = ceil(height/dh);

    //init
    body3d.push_back(vector<vector<double>>());
    body3d.push_back(vector<vector<double>>());
    body3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());

    for(int k = 0; k < nh; k++){
        double rad = k*dh*tan(alpha);
        if(rad!=0.0){
            double dfi = step/rad;
            int nfi = floor(2*M_PI/dfi);
            for(int l = 0; l < nfi; l++){

                vector<double> normal = {0,0,0};
                vector<double> initNorm = {cos(alpha), 0, -sin(alpha)};
                rotBody('z', initNorm, normal, l*dfi);

                if(cos(l*dfi)>0){
                    body3d[0].push_back({rad*cos(l*dfi), rad*sin(l*dfi), (k+1)*dh});
                    norm3d[0].push_back(normal);
                }
                else{
                    body3d[1].push_back({rad*cos(l*dfi), rad*sin(l*dfi), (k+1)*dh});
                    norm3d[1].push_back(normal);
                }
            }
        }
    }

    double nPts = body3d[0].size()+body3d[1].size();
    double square = M_PI*radius*sqrt(radius*radius+height*height)/nPts;

    for(int i = 0; i < body3d.size(); i++){
        for(int j = 0; j < body3d[i].size(); j++){
            for(int k = 0; k < 3; k++){
                norm3d[i][j][k] *= square;
            }
        }
    }

    Disk disk(radius, step);
    copy(disk.body3d[1].begin(), disk.body3d[1].end(), back_inserter(body3d[2]));
    copy(disk.norm3d[1].begin(), disk.norm3d[1].end(), back_inserter(norm3d[2]));
    for(int i = 0; i < body3d[2].size(); i++){
        body3d[2][i][2] = height;
    }
}

Cylinder::Cylinder(double h, double r, double dheight){

    height = h;
    radius = r;
    step = dheight;

    double dfi = step/radius;
    int nh = floor(height/step-1);
    int nfi = 2*M_PI/dfi;

    //init
    body3d.push_back(vector<vector<double>>());
    body3d.push_back(vector<vector<double>>());
    body3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());

    //edge points
    for(int k = 0; k < nfi; k++){
        for(int l = 0; l < nh; l++){
            if(cos(k*dfi)<0){
                body3d[0].push_back({radius*cos(k*dfi), radius*sin(k*dfi),-height/2+(l+1)*step});
                norm3d[0].push_back({cos(k*dfi), sin(k*dfi), 0});
            }
        }
    }

    //normalize
    double nPts = body3d[0].size();
    double square = M_PI*radius*height/nPts;
    for(int i = 0; i < body3d[0].size(); i++){
        for(int j = 0; j < 3; j++){
            norm3d[0][i][j] *= square;
        }
    }

    //top and bottom
    Disk disk(radius, step);
    copy(disk.body3d[0].begin(), disk.body3d[0].end(), back_inserter(body3d[1]));
    copy(disk.norm3d[0].begin(), disk.norm3d[0].end(), back_inserter(norm3d[1]));
    copy(disk.body3d[1].begin(), disk.body3d[1].end(), back_inserter(body3d[2]));
    copy(disk.norm3d[1].begin(), disk.norm3d[1].end(), back_inserter(norm3d[2]));

    for(int i = 0; i < body3d[1].size(); i++){
        body3d[1][i][2] = -height/2;
    }

    for(int i = 0; i < body3d[2].size(); i++){
        body3d[2][i][2] = height/2;
    }

}

Rocket::Rocket(double h1, double h2, double r, double dheight){

    height1 = h1;
    height2 = h2;
    radius = r;
    step = dheight;

    //init
    body3d.push_back(vector<vector<double>>());
    body3d.push_back(vector<vector<double>>());
    body3d.push_back(vector<vector<double>>());
    body3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());
    norm3d.push_back(vector<vector<double>>());


    Cone cone(height1,radius,step);
    Cylinder cyl(height2, radius, step);

    copy(cone.body3d[0].begin(), cone.body3d[0].end(), back_inserter(body3d[0]));
    copy(cone.norm3d[0].begin(), cone.norm3d[0].end(), back_inserter(norm3d[0]));
    copy(cone.body3d[1].begin(), cone.body3d[1].end(), back_inserter(body3d[1]));
    copy(cone.norm3d[1].begin(), cone.norm3d[1].end(), back_inserter(norm3d[1]));

    copy(cyl.body3d[0].begin(), cyl.body3d[0].end(), back_inserter(body3d[2]));
    copy(cyl.norm3d[0].begin(), cyl.norm3d[0].end(), back_inserter(norm3d[2]));
    copy(cyl.body3d[2].begin(), cyl.body3d[2].end(), back_inserter(body3d[3]));
    copy(cyl.norm3d[2].begin(), cyl.norm3d[2].end(), back_inserter(norm3d[3]));

    for(int i = 2; i < 4; i++){
        for(int j = 0; j < body3d[i].size(); j++){
            body3d[i][j][2] += (height2/2+height1);
        }
    }

}

void Cone::getVisiblePts(double angle){

    double alpha = atan(radius/height);

    if (angle < alpha){

        visible = {1, 1, 0};

    }else if (angle <= M_PI/2){

        visible = {0, 1, 0};

    }else if (angle < M_PI){

        visible = {0, 1, 1};

    }else if (angle >= M_PI){

        visible = {0, 0, 1};
    }

}

void Disk::getVisiblePts(double angle){
    if(angle < M_PI/2){
        visible = {1,0};
	}
	else{
        visible = {0,1};
    }
}

void Cylinder::getVisiblePts(double angle){

    if (angle == 0){
        visible = {0, 1, 0};
    }
    else if (angle == M_PI/2){
        visible = {1, 0, 0};
    }
    else if (angle >= M_PI){
        visible = {0, 0, 1};
    }
    else if (angle < M_PI/2){
        visible = {1, 1, 0};
    }
    else if (angle < M_PI){
        visible = {1, 0, 1};
    }
}

void Rocket::getVisiblePts(double angle){

    double alpha = atan(radius/height1);

    if (angle == 0.0){
        visible = {1, 1, 0, 0};
    }
    else if(angle < alpha){
        visible = {1, 1, 1, 0};
    }
    else if(angle < M_PI/2){
        visible = {0, 1, 1, 0};
    }
    else if(angle < M_PI-alpha){
        visible = {0, 1, 1, 1};
    }
    else if(angle < M_PI){
        visible = {0, 0, 1, 1};
    }
    else if(angle==M_PI){
        visible = {0, 0, 0, 1};
    }

}

vector<vector<vector<double>>> VolumeObject::getRCS(double lambda, vector<double>& angles, char mode){

	double lightspeed = 299792458;

	double k_abs = 2 * M_PI / lambda;
	vector<double> k_vec = { 0, 0, k_abs };

	double B_abs = 1;
	vector<double> B_vec;
	if (mode == 'H'){
		B_vec = { 0, B_abs, 0 };
	}
	else if (mode == 'V'){
		B_vec = { B_abs, 0, 0 };
	}

	vector<double> E_vec = cross(B_vec, k_vec);
	for (int i = 0; i < E_vec.size(); i++){
		E_vec[i] *= lightspeed / k_abs;
	}

	double E_abs = sqrt(dot(E_vec, E_vec));

	double R_abs = 1e3;
	vector<double> R_vec = { 0, 0, -R_abs };

	int N = angles.size();

	vector<vector<vector<double>>> E_field(N,
		vector<vector<double>>(3, vector<double>(2,0)));

	vector<double> coefE = { -lightspeed / (lambda*R_abs)*sin(k_abs*R_abs),
		lightspeed / (lambda*R_abs)*cos(k_abs*R_abs) };

	for (int i = 0; i < N; i++){

		vector<double> B(3, 0);
		vector<double> R_(3, 0);
		rotBody('y', B_vec, B, angles[i]);
		rotBody('y', R_vec, R_, angles[i]);
		for (int j = 0; j < 3; j++){
			R_[j] /= R_abs;
		}

		getVisiblePts(angles[i]);

		vector<double> crossE;
		vector<double> expon;
		double phase;
		vector<vector<double>> Etemp(3, vector<double>(2, 0));

		for (int iVis = 0; iVis < visible.size(); iVis++){
			if (visible[iVis] == 1){
				for (int j = 0; j < norm3d[iVis].size(); j++){

					crossE = cross(norm3d[iVis][j], B);

					phase = 2 * k_abs*dot(body3d[iVis][j], R_);
					expon = { cos(phase), -sin(phase) };

					for (int k = 0; k < 3; k++){
						Etemp[k][0] += crossE[k] * expon[0];
						Etemp[k][1] += crossE[k] * expon[1];
					}

				}
			}
		}

		for (int j = 0; j < 3; j++){
			Etemp[j] = complexMult(Etemp[j], coefE);
			E_field[i][j][0] = Etemp[j][0];
			E_field[i][j][1] = Etemp[j][1];
		}
	}

	return E_field;
}