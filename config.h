#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <cstdlib>

class Config {
public:
	int getnx() { return nx; }
	int getny() { return ny; }
	int getnz() { return nz; }

	void setnx(int _nx) { this->nx = _nx; }
	void setny(int _ny) { this->ny = _ny; }
	void setnz(int _nz) { this->nz = _nz; }


	//double get() { return; }
	double getWallXPosition() { return  wallXPosition; }
	double getdzdy() { return dzdy; }
	double getWallLength() { return wallLength; }
	double getLambda() { return lambda; }
	double getA() { return A; }
	double getk() { return k; }
	double getOmega() { return omega; }
	double getTheta() { return theta; }
	double getPhi() { return phi; }
	double getPsi() { return psi; }

	std::string getSaveString() { return saveString; }

	//void set(double _) { this-> = _; }
	void setWallXPosition(double _wallXPosition) { this->wallXPosition = _wallXPosition; }
	void setdzdy(double _dzdy) { this->dzdy = _dzdy; }
	void setWallLength(double _wallLength) { this->wallLength = _wallLength; }
	void setLambda(double _lambda) { this->lambda = _lambda; }
	void setA(double _A) { this->A = _A; }
	void setk(double _k) { this->k = _k; }
	void setOmega(double _omega) { this->omega = _omega; }
	void setTheta(double _theta) { this->theta = _theta; }
	void setPhi(double _phi) { this->phi = _phi; }
	void setPsi(double _psi) { this->psi = _psi; }

	void setCellName(std::string _setcellName) { this->cellName = _setcellName; }
	void setCellType(std::string _setcellType) { this->cellType = _setcellType; }
	void setCellCentering(std::string _cellCentering) { this->cellCentering = _cellCentering; }
	
	void setManualConfig();
	void setFullConfig(double _wallXPosition, double _dzdy, double _wallLength, double _lambda, double _A, double _k, double _omega, double _theta, double _phi, double _psi, std::string _cellCentering, std::string _cellName, std::string _cellType);

	std::string getCellName() { return cellName; }
	std::string getCellType() { return cellType; }
	std::string getCellCentering() { return cellCentering; }

	void setSaveString(std::string _saveString) { this->saveString = _saveString; }

private:
	int nx = 1;
	int ny = 1;
	int nz = 1;

	double wallXPosition = 0.01;// wall at x = wallXPos
	double dzdy = 0.0001;// dzdy->step size for wall spacing
	double wallLength = .1; // wall is wallLen x wallLen meters
	double lambda = 1e-9; // wavelength (meters)
	double A = 1;  // Amplitude (arbitrary units)
	double k = 1 ; // wavenumber (1/m)
	double omega = 0; //
	double theta = 0;
	double phi = 0;
	double psi = 0;


	// Defaults 
	std::string cellName = "salt";
	std::string cellType = "";
	std::string cellCentering = "face-centered";
	std::string saveString = "C:\\Users\\Michael\\Documents\\Programming\\laueDiffractionResults\\EigenResultsFaceHighRes.csv";

};

#endif 