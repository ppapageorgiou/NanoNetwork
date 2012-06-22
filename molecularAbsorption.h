#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>
#include "engine.h"

#define PI 3.1415926536
#define e 2.71828183


//MOLECULE STRUCTURE
struct molecule {
	std::vector<double> igas;
	std::vector<double> iso;
	std::vector<double> wnum;
	std::vector<double> inti;
	std::vector<double> Acoeff;
	std::vector<double> abroad;
	std::vector<double> sbroad;
	std::vector<double> els;
	std::vector<double> abcoef;
	std::vector<double> tsp;
	std::vector<double> gn;
};

	//FUNCTION PROTOTYPES
	std::vector<std::vector<double>> molecularAbsorption(double P, double T, double conc);
	void importData(std::string name, int size, int index, molecule &O2, molecule &CO2, molecule &CH4,
					molecule &N2O, molecule &O3, molecule &NO2, molecule &CO, molecule &H2O);
	int countLines(std::string name);
	std::vector<double> parseData(std::string line);
	void printMolecule(molecule name);
	molecule chooseMolecule(int index, molecule &O2, molecule &CO2, molecule &CH4,
							molecule &N2O, molecule &O3, molecule &NO2, molecule &CO, molecule &H2O);
	void vectorToFile(std::vector<std::vector<double>> vector);
	std::vector<double> chooseKvVector (int index, std::vector<double> &O2kv, std::vector<double> &CO2kv, std::vector<double> &CH4kv,
									std::vector<double> &N2Okv, std::vector<double> &O3kv, std::vector<double> &NO2kv,
									std::vector<double> &COkv, std::vector<double> &H2Okv);
	//std::vector<double> sumCoeffs (std::vector<std::vector<double>> vec);
	//void sendMATLAB (std::vector<double> result);