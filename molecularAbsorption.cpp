// TERAHERTZ BAND CHANNEL MODEL (100GHz - 10THz)
// Author: Panayiotis Papageorgiou 2012
// Imperial College London
// Supervisor: Dr. Bruno Clerckx
// 
// DESCRIPTION
// ===========
// This is an attempt to construct a terahertz channel model from HITRAN
// data as described in the paper Channel Modeling and Capacity Analysis 
// for Electromagnetic Wireless Nanonetworks in the Terahertz Band by
// Josep Miquel Jornet, Student Member, IEEE, and Ian F. Akyildiz, Fellow, IEEE
//
// Special Thanks: Special thanks to Josep Miquel Jornet for his valuable help
// in understanding the implementation of the equations in the above paper.

//INCLUDED LIBRARIES

#include "molecularAbsorption.h"


std::vector<std::vector<double>> molecularAbsorption(double P, double T, double conc)
{
	//INITIALIZE CONSTANTS
	double NA = 6.022141E+23;
	double Tstp = 273.15; 
	double T0 = 296;
	double P0 = 1;
	double k = 1.380658E-23;
	double h = 6.62618E-34;
	double c = 2.997925E+8;
	double R = 82.0575;

	//SYSTEM PARAMETERS
	//double P = 1;
	//double T = 296;
	//INITIALIZE SPECTRA
	double bands[] = {2011,331,9720,446,102440,32064,500,4765}; //Number of lines per molecule
	double con[] = {209460,383,1.745,0.3,0.05,0.02,0.1,conc*10000}; //Concentration of each molecule in the atmosphere
	double mass[] = {32/NA,44.01/NA,16.04/NA,44.01/NA,48/NA,46.01/NA,28.01/NA,18.02/NA}; //Molar mass

	//INITIALIZE MOLECULE STRUCRURES
	molecule O2, CO2, CH4, N2O, O3, NO2, CO, H2O;

	//INITIALIZE PARAMETER ARRAY
	std::vector<std::vector<double> > kvt ( 8, std::vector<double> ( 10000 ) );
	for ( int i = 0; i < 8; i++ ) {
		for ( int j = 0; j < 10000; j++ )
			kvt[i][j] = 0;
	}

	std::vector<double> v;
	v.push_back((double) 3);
	double offset=(double)330/(double)9999;
	for (int i=1;i<10000;i++) {
		v.push_back((v[i-1])+offset);
	}
	
	//IMPORT HITRAN PARAMETERS
	for (int i=1;i<9;i++) {
		std::string ext=".out";
		std::string fileName = static_cast<std::ostringstream*>( &(std::ostringstream() << i) )->str();
		fileName.append(ext);
		int size=countLines(fileName);
		importData(fileName, size, i, O2, CO2, CH4, N2O, O3, NO2, CO, H2O);		
	}

	//INITIALZIE ABSORBTION COEFFICIENT VECTOR
	std::vector<double> coeffSum;

	//PARAMETER PROCESSING
	double q;
	double qq;
	double m;
	double vc;
	double S0;
	double gamma_air;
	double gamma_self;
	double n;
	double alpha_l;
	double fv;
	double sigmav;
	double kv;


	for (int i=0;i<8;i++) {
		q=(con[i]*1E-6);
		qq=q*P/R/T*NA;
		m=mass[i];
		molecule current=chooseMolecule(i, O2, CO2, CH4, N2O, O3, NO2, CO, H2O);
		std::cout<<current.igas.size()<<" Lines found for current molecule"<<std::endl<<std::endl;

		for (int j=0;j<current.igas.size();j++) {
			//Set required parameters from HITRAN Data
			vc = current.wnum[j];
			S0 = current.inti[j];
			gamma_air = current.abroad[j];
			gamma_self = current.sbroad[j];
			n = current.abcoef[j];
			//Calculate the Lorentz half-width
			alpha_l = ((1-q)*gamma_air+q*gamma_self)*(P/P0)*pow((T/T0),n);
			for (int k=0;k<10000;k++) {
				//Vleck-Weisskopf assymetric line shape
				fv = alpha_l/PI*pow((v[k]/vc),2)*((1/(pow((v[k]-vc),2)+pow(alpha_l,2)))+(1/(pow((v[k]+vc),2)+pow(alpha_l,2))));
				//Adjusting the far ends of the line shape
				sigmav = (v[k]/vc*tanh(h*c*v[k]/(2*k*T))/tanh(h*c*vc/(2*k*T))*fv)*S0;
				kv = P/P0*T0/T*qq*sigmav;
				//Output to Absorption Coefficient Array
				kvt[i][k] = kvt[i][k] + kv;
			}
		}
	}
	return kvt;
}

void importData(std::string name, int size, int index, molecule &O2, molecule &CO2, molecule &CH4,
				molecule &N2O, molecule &O3, molecule &NO2, molecule &CO, molecule &H2O) {

  std::string* line=new std::string[size];
  std::string buffer;
  std::vector<double> lineResults;
  std::ifstream myfile(name);
  if (myfile.is_open())
  {
    while ( myfile.good() )
    {
		getline(myfile,buffer);
		lineResults=parseData(buffer);
		if (index==1) {
				O2.igas.push_back(lineResults[0]);
				O2.iso.push_back(lineResults[1]);
				O2.wnum.push_back(lineResults[2]);
				O2.inti.push_back(lineResults[3]);
				O2.Acoeff.push_back(lineResults[4]);
				O2.abroad.push_back(lineResults[5]);
				O2.sbroad.push_back(lineResults[6]);
				O2.els.push_back(lineResults[7]);
				O2.abcoef.push_back(lineResults[8]);
				O2.tsp.push_back(lineResults[9]);
				O2.gn.push_back(lineResults[10]);
		}
		if (index==2) {
				CO2.igas.push_back(lineResults[0]);
				CO2.iso.push_back(lineResults[1]);
				CO2.wnum.push_back(lineResults[2]);
				CO2.inti.push_back(lineResults[3]);
				CO2.Acoeff.push_back(lineResults[4]);
				CO2.abroad.push_back(lineResults[5]);
				CO2.sbroad.push_back(lineResults[6]);
				CO2.els.push_back(lineResults[7]);
				CO2.abcoef.push_back(lineResults[8]);
				CO2.tsp.push_back(lineResults[9]);
				CO2.gn.push_back(lineResults[10]);
		}
		if (index==3) {
				CH4.igas.push_back(lineResults[0]);
				CH4.iso.push_back(lineResults[1]);
				CH4.wnum.push_back(lineResults[2]);
				CH4.inti.push_back(lineResults[3]);
				CH4.Acoeff.push_back(lineResults[4]);
				CH4.abroad.push_back(lineResults[5]);
				CH4.sbroad.push_back(lineResults[6]);
				CH4.els.push_back(lineResults[7]);
				CH4.abcoef.push_back(lineResults[8]);
				CH4.tsp.push_back(lineResults[9]);
				CH4.gn.push_back(lineResults[10]);
		}
		if (index==4) {
				N2O.igas.push_back(lineResults[0]);
				N2O.iso.push_back(lineResults[1]);
				N2O.wnum.push_back(lineResults[2]);
				N2O.inti.push_back(lineResults[3]);
				N2O.Acoeff.push_back(lineResults[4]);
				N2O.abroad.push_back(lineResults[5]);
				N2O.sbroad.push_back(lineResults[6]);
				N2O.els.push_back(lineResults[7]);
				N2O.abcoef.push_back(lineResults[8]);
				N2O.tsp.push_back(lineResults[9]);
				N2O.gn.push_back(lineResults[10]);
		}
		if (index==5) {
				O3.igas.push_back(lineResults[0]);
				O3.iso.push_back(lineResults[1]);
				O3.wnum.push_back(lineResults[2]);
				O3.inti.push_back(lineResults[3]);
				O3.Acoeff.push_back(lineResults[4]);
				O3.abroad.push_back(lineResults[5]);
				O3.sbroad.push_back(lineResults[6]);
				O3.els.push_back(lineResults[7]);
				O3.abcoef.push_back(lineResults[8]);
				O3.tsp.push_back(lineResults[9]);
				O3.gn.push_back(lineResults[10]);
		}
		if (index==6) {
				NO2.igas.push_back(lineResults[0]);
				NO2.iso.push_back(lineResults[1]);
				NO2.wnum.push_back(lineResults[2]);
				NO2.inti.push_back(lineResults[3]);
				NO2.Acoeff.push_back(lineResults[4]);
				NO2.abroad.push_back(lineResults[5]);
				NO2.sbroad.push_back(lineResults[6]);
				NO2.els.push_back(lineResults[7]);
				NO2.abcoef.push_back(lineResults[8]);
				NO2.tsp.push_back(lineResults[9]);
				NO2.gn.push_back(lineResults[10]);
		}
		if (index==7) {
				CO.igas.push_back(lineResults[0]);
				CO.iso.push_back(lineResults[1]);
				CO.wnum.push_back(lineResults[2]);
				CO.inti.push_back(lineResults[3]);
				CO.Acoeff.push_back(lineResults[4]);
				CO.abroad.push_back(lineResults[5]);
				CO.sbroad.push_back(lineResults[6]);
				CO.els.push_back(lineResults[7]);
				CO.abcoef.push_back(lineResults[8]);
				CO.tsp.push_back(lineResults[9]);
				CO.gn.push_back(lineResults[10]);
		}
		if (index==8) {
				H2O.igas.push_back(lineResults[0]);
				H2O.iso.push_back(lineResults[1]);
				H2O.wnum.push_back(lineResults[2]);
				H2O.inti.push_back(lineResults[3]);
				H2O.Acoeff.push_back(lineResults[4]);
				H2O.abroad.push_back(lineResults[5]);
				H2O.sbroad.push_back(lineResults[6]);
				H2O.els.push_back(lineResults[7]);
				H2O.abcoef.push_back(lineResults[8]);
				H2O.tsp.push_back(lineResults[9]);
				H2O.gn.push_back(lineResults[10]);
		}
    }
    myfile.close();
  }

  else std::cout << "Unable to open file" << std::endl; 

}

int countLines(std::string name) {

	std::ifstream myfile(name);
	if (myfile.is_open()) {
		std::string temp;
		int count=0;
		while ( !myfile.eof() )
		{
			getline(myfile,temp);
			count++;
		}
		myfile.close();
		return count;
	}
	std::cout<<"Unable to open file" << std::endl;
	return 0;
}

std::vector<double> parseData(std::string line) {

	double temp;
	std::vector<double> data;
	std::string buffer;
	buffer=line.substr(1,1);
	if ( ! (std::istringstream(buffer) >> temp) ) temp = 0;
	data.push_back(temp);

	buffer=line.substr(2,1);
	if ( ! (std::istringstream(buffer) >> temp) ) temp = 0;
	data.push_back(temp);

	buffer=line.substr(3,12);
	if ( ! (std::istringstream(buffer) >> temp) ) temp = 0;
	data.push_back(temp);

	buffer=line.substr(15,10);
	if ( ! (std::istringstream(buffer) >> temp) ) temp = 0;
	data.push_back(temp);

	buffer=line.substr(25,10);
	if ( ! (std::istringstream(buffer) >> temp) ) temp = 0;
	data.push_back(temp);

	buffer=line.substr(35,5);
	if ( ! (std::istringstream(buffer) >> temp) ) temp = 0;
	data.push_back(temp);

	buffer=line.substr(40,5);
	if ( ! (std::istringstream(buffer) >> temp) ) temp = 0;
	data.push_back(temp);

	buffer=line.substr(45,10);
	if ( ! (std::istringstream(buffer) >> temp) ) temp = 0;
	data.push_back(temp);

	buffer=line.substr(55,4);
	if ( ! (std::istringstream(buffer) >> temp) ) temp = 0;
	data.push_back(temp);

	buffer=line.substr(59,8);
	if ( ! (std::istringstream(buffer) >> temp) ) temp = 0;
	data.push_back(temp);

	buffer=line.substr(156,4);
	if ( ! (std::istringstream(buffer) >> temp) ) temp = 0;
	data.push_back(temp);

	return data;
}

void printMolecule(molecule name) {

	std::cout<<" >> "<<name.igas[0]<<" | ";
	std::cout<<name.iso[0]<<" | ";
	std::cout<<name.wnum[0]<<" | ";
	std::cout<<name.inti[0]<<" | ";
	std::cout<<name.Acoeff[0]<<" | ";
	std::cout<<name.abroad[0]<<" | ";
	std::cout<<name.sbroad[0]<<" | ";
	std::cout<<name.els[0]<<" | ";
	std::cout<<name.abcoef[0]<<" | ";
	std::cout<<name.tsp[0]<<" | ";
	std::cout<<name.gn[0]<<" << ";

}

molecule chooseMolecule(int index, molecule &O2, molecule &CO2, molecule &CH4,
						molecule &N2O, molecule &O3, molecule &NO2, molecule &CO, molecule &H2O) {

	switch (index) {

	case 0 : 
		std::cout<<"Processing Molecule: O2..."<<std::endl;
		return O2;

	case 1 :
		std::cout<<"Processing Molecule: CO2..."<<std::endl;
		return CO2;

	case 2 :
		std::cout<<"Processing Molecule: CH4..."<<std::endl;
		return CH4;

	case 3 :
		std::cout<<"Processing Molecule: N2O..."<<std::endl;
		return N2O;

	case 4 :
		std::cout<<"Processing Molecule: O3..."<<std::endl;
		return O3;

	case 5 :
		std::cout<<"Processing Molecule: NO2..."<<std::endl;
		return NO2;

	case 6 :
		std::cout<<"Processing Molecule: CO..."<<std::endl;
		return CO;

	case 7 :
		std::cout<<"Processing Molecule: H2O..."<<std::endl;
		return H2O;

	  default : 
		std::cout<<"Invalid Molecule Selected"<<std::endl;

	}

}

void vectorToFile(std::vector<std::vector<double>> vector) {

	std::ofstream myfile;
	myfile.open ("results.txt", std::ios::app);
	for ( int i = 0; i < 8; i++ ) {
		for ( int j = 0; j < 10000; j++ ) {
			myfile << vector[i][j] << "\n";
		}
		myfile<<"============"<<"\n";
	}
	myfile.close();
	std::cout<<("Succesfully Printed to file!\n");
}