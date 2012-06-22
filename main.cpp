#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <random>
#include <time.h>
#include "engine.h"
#include "molecularAbsorption.h"

#define DATA_LENGTH 32

using namespace std;

class nanoDevice 
{
public:
	int ID;
	double x,y,z;
	int data[DATA_LENGTH];
	nanoDevice(int, std::tr1::mt19937&);
	int printInfo();
	double distance();
};

nanoDevice::nanoDevice(int devID, std::tr1::mt19937 &eng) {

	ID=devID;
	std::tr1::poisson_distribution<double> poisson(5);
	x=poisson(eng)*1E-1;
	y=poisson(eng)*1E-1;
	z=poisson(eng)*1E-1;

	//Initialize Data Array
	for (int i=0;i<(sizeof(data)/sizeof(int));i++) {
		data[i] = rand() % 2;
	}

}

int nanoDevice::printInfo() {

	cout<<"Device ID: "<<ID<<"   ";
	cout<<"X-Pos: "<<x<<"cm   ";
	cout<<"Y-Pos: "<<y<<"cm   ";
	cout<<"Z-Pos: "<<z<<"cm   "<<endl;
	cout<<"Data: ";
	for (int i=0;i<(sizeof(data)/sizeof(int));i++) {
		cout<<data[i];
	}
	cout<<endl<<endl;

	return ID;
}

double nanoDevice::distance() {

	return (sqrt(pow(x,2)+pow(y,2)+pow(z,2)));

}


vector<nanoDevice*> initializeDevices (int num);
void sendDevMATLAB (std::vector<nanoDevice*> devices);
void sendCoeffMATLAB (vector<double> result, int mode);
vector<double> modulate (int data);
double max(vector<double> vec);
vector<double> sumCoeffs (vector<vector<double>> vec, int mode);
double calcNoise (vector<double> coeffs, double lowerB, double upperB);
vector<double> pathLoss (vector<double> coeffs);
void noisePowerSim (vector<double> coeffs);


void main() {

	//SET SIMULATION PARAMETERS
	int modeMolec = 9;
	int modePlot = 1;
	double Temp = 296; //In degrees Kelvin
	double Press = 1; //In atm
	double conc = 0; //As a percentage
	double lowerB = 3E+12;
	double upperB = 4E+12;

	vector<vector<double>> absCoeffs;
	vector<double> coeffs;
	vector<nanoDevice*> res;
	res=initializeDevices(100);
	sendDevMATLAB(res);

	cout<<endl<<endl<<"Total Number of Devices Initialized: "<<res.size()<<endl<<endl;

	absCoeffs=molecularAbsorption(Press, Temp, conc);

	coeffs=sumCoeffs(absCoeffs, modeMolec);

	if (modePlot==1) {
		sendCoeffMATLAB(coeffs, modePlot);
	}
	if (modePlot==2) {
		vector<double> atten;
		for (int i=0;i<10000;i++) {
			atten.push_back(coeffs[i]*10*log10((double)e));
		}
		sendCoeffMATLAB(atten,modePlot);
	}
	if (modePlot==3) {
		vector<double> atten;
		for (int i=0;i<10000;i++) {
			atten.push_back(296-296*(pow(e,(-1*coeffs[i]))));
		}
		modePlot=2;
		sendCoeffMATLAB(atten,modePlot);
	}
	if (modePlot==4) {
		vector<double> loss;
		loss=pathLoss(coeffs);
		modePlot=2;
		sendCoeffMATLAB(loss,modePlot);
	}

	//cout<<"The Noise Power at the Reciever is: "<<calcNoise(coeffs, lowerB, upperB)<<" dB"<<endl<<endl;
	noisePowerSim(coeffs);

	system("PAUSE");

}

vector<nanoDevice*> initializeDevices (int num) {

	vector<nanoDevice*> devices;
	std::tr1::mt19937 eng;
	eng.seed((unsigned int)time(NULL));

	for (int i=0;i<num;i++) {
		devices.push_back(new nanoDevice(i,eng));
		if (num<=20) {
			devices[i]->printInfo();
		}
	}
	return devices;
}

void sendDevMATLAB (std::vector<nanoDevice*> devices) {
	
	Engine *ep;
	mxArray *X = NULL;
	mxArray *Y = NULL;
	mxArray *Z = NULL;
	double* xMat = new double[devices.size()];
	double* yMat = new double[devices.size()];
	double* zMat = new double[devices.size()];

	for (int i=0;i<devices.size();i++) {
		xMat[i]=devices[i]->x;
		yMat[i]=devices[i]->y;
		zMat[i]=devices[i]->z;
	}

	//Start the MATLAB engine 
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
	}

	X = mxCreateDoubleMatrix(1, devices.size(), mxREAL);
	Y = mxCreateDoubleMatrix(1, devices.size(), mxREAL);
	Z = mxCreateDoubleMatrix(1, devices.size(), mxREAL);
	memcpy((void *)mxGetPr(X), (void *)xMat, sizeof(double)*devices.size());
	memcpy((void *)mxGetPr(Y), (void *)yMat, sizeof(double)*devices.size());
	memcpy((void *)mxGetPr(Z), (void *)zMat, sizeof(double)*devices.size());

	engPutVariable(ep, "X", X);
	engPutVariable(ep, "Y", Y);
	engPutVariable(ep, "Z", Z);

	engEvalString(ep, "figure");
	engEvalString(ep, "scatter3(X(:),Y(:),Z(:),'filled')");

}

void sendCoeffMATLAB (vector<double> result, int mode) {
	
	Engine *ep;
	mxArray *mat = NULL;
	double cMat[10000];

	for (int i=0;i<10000;i++) {
		cMat[i] = result[i];
	}

	//Start the MATLAB engine 
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
	}

	mat = mxCreateDoubleMatrix(1, 10000, mxREAL);
	memcpy((void *)mxGetPr(mat), (void *)cMat, sizeof(cMat));

	engPutVariable(ep, "mat", mat);

	engEvalString(ep, "figure");
	engEvalString(ep, "c = 2.997925e8");
	engEvalString(ep, "v = linspace(3,333,10000)");
	if (mode == 1 || 3) {
		engEvalString(ep, "plot(v*c*1e2,sum(mat,1))");
		engEvalString(ep, "ylabel('Absorption Coefficients [/m]')");
	}
	if (mode == 2) {
		engEvalString(ep, "semilogy(v*c*1e2,sum(mat,1))");
		engEvalString(ep, "ylabel('Signal Attenuation [dB/m]')");
	}
	engEvalString(ep, "xlabel('Frequency [Hz]')");
	
}

vector<double> modulate (int data[]) {

	vector<double> t, y, p, pulse;
	double temp;
	
	//CONSTANTS
	double F=30E+12;
	double T=5E-12;
	double E=10;

	return pulse;
}

double max(vector<double> vec) {

	double res=0;

	for (int i=0; i<vec.size();i++) {
		if (abs(vec[i])>res)
			res=abs(vec[i]);
	}

	return res;
}

vector<double> sumCoeffs (vector<vector<double>> vec, int mode) {
	
	vector<double> sum;
	//INITIALIZE SUM VECTOR
	for (int i=0;i<10000;i++) sum.push_back(0);
	

	switch (mode) {

	case 1 :
	for (int i=0;i<10000;i++) {
		sum[i]=(sum[i]+vec[0][i]);
	}
	return sum;

	case 2 :
	for (int i=0;i<10000;i++) {
		sum[i]=(sum[i]+vec[1][i]);
	}
	return sum;

	case 3 :
	for (int i=0;i<10000;i++) {
		sum[i]=(sum[i]+vec[2][i]);
	}
	return sum;

	case 4 :
	for (int i=0;i<10000;i++) {
		sum[i]=(sum[i]+vec[3][i]);
	}
	return sum;

	case 5 :
	for (int i=0;i<10000;i++) {
		sum[i]=(sum[i]+vec[4][i]);
	}
	return sum;

	case 6 :
	for (int i=0;i<10000;i++) {
		sum[i]=(sum[i]+vec[5][i]);
	}
	return sum;

	case 7 :
	for (int i=0;i<10000;i++) {
		sum[i]=(sum[i]+vec[6][i]);
	}
	return sum;

	case 8 :
	for (int i=0;i<10000;i++) {
		sum[i]=(sum[i]+vec[7][i]);
	}
	return sum;

	case 9 :
	for (int i=0;i<10000;i++) {
		for (int j=0;j<8;j++)
			sum[i]=(sum[i]+vec[j][i]);
	}
	return sum;

	default : 
	cout<<"Invalid Selection..."<<endl;
	break;
	}


}

double calcNoise (vector<double> coeffs, double lowerB, double upperB) {

	if (lowerB<100E+9 || lowerB>10E+12) {
		cout<<"Lower Bandwidth limit must be between 100 GHz and 10 THz... Aborting Noise Calculation"<<endl;
		return 0;
	}
	if (upperB<100E+9 || upperB>10E+12) {
		cout<<"Upper Bandwidth limit must be between 100 GHz and 10 THz... Aborting Noise Calculation"<<endl;
		return 0;
	}
	if (lowerB>upperB) {
		cout<<"Lower Bandwidth limit is greater than Upper Bandwidth limit... Aborting Noise Calculation"<<endl;
		return 0;
	}
	
	double kB = 1.3806E-23;
	int indexLowerB, indexUpperB;
	double noise;
	double offset = ((10E+12)-(100E+9))/9999;
	indexLowerB = (lowerB/offset)-101;
	indexUpperB = (upperB/offset)-101;	

	for (int i=indexLowerB;i<=indexUpperB;i++) {

		noise+=kB*(296-296*(pow(e,(-coeffs[i]))));

	}

	return (10*log10(noise));

}

vector<double> pathLoss (vector<double> coeffs) {

	//PREPARE FREQUENCY LINEAR SPACE
	std::vector<double> v;
	v.push_back((double) 3);
	double offset=(double)330/(double)9999;
	for (int i=1;i<10000;i++) {
		v.push_back((v[i-1])+offset);
	}

	double spread;
	double atten;
	vector<double> pathLoss;

	for (int i=0;i<10000;i++) {
		spread=20*log((4*PI*v[i]*100));
		atten=coeffs[i]*10*log10((double)e);
		pathLoss.push_back(spread+atten);
	}

	return pathLoss;

}

void noisePowerSim (vector<double> coeffs) {

	double bandArray[12] = { 100E+9, 300E+9, 1E+12, 2E+12, 3E+12, 4E+12, 5E+12, 6E+12, 7E+12, 8E+12, 9E+12, 10E+12 };
	double lowerB, upperB, noise;

	std::ofstream myfile;
	myfile.open ("noisePower.txt", std::ios::app);

	for (int i=0;i<11;i++) {
		lowerB=bandArray[i];
		upperB=bandArray[i+1];
		noise=calcNoise(coeffs, lowerB, upperB);
		cout<<"The Noise Power at the Reciever is: "<<noise<<" dB For the Range: "<<lowerB<<" - "<<upperB<<endl;
		myfile<<noise<<"\n";
	}

	myfile.close();
	cout<<endl;

}