#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
#include <Eigen/Dense>

class Config {
public:
	double getnx() { return nx; }
	double getny() { return ny; }
	double getnz() { return nz; }

	void setnx(double _nx) { this->nx = _nx; }
	void setny(double _ny) { this->ny = _ny; }
	void setnz(double _nz) { this->nz = _nz; }


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

	std::string getCellName() { return cellName; }
	std::string getCellType() { return cellType; }
	std::string getCellCentering() { return cellCentering; }

	void setCellName(double _setcellName) { this->cellName = _setcellName; }
	void setCellType(double _setcellType) { this->cellType = _setcellType; }
	void setCellCentering(double _cellCentering) { this->cellCentering = _cellCentering; }

private:
	int nx = 1;
	int ny = 1;
	int nz = 1;

	double wallXPosition = 0.01;// wall at x = wallXPos
	double dzdy = 0.0005;// dzdy->step size for wall spacing
	double wallLength = .2; // wall is wallLen x wallLen meters
	double lambda = 1e-9; // wavelength (meters)
	double A;  // Amplitude (arbitrary units)
	double k; // wavenumber (1/m)
	double omega = 0; //
	double theta = 0;
	double phi = 0;
	double psi = 0;

	std::string cellName = "salt";
	std::string cellType = "";
	std::string cellCentering = "face-centered";
};

#endif 