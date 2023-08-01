#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <cstdlib>

class Config {
public:
	const int getnx() const { return nx; }
	const int getny() const { return ny; }
	const int getnz() const { return nz; }

	void setnx(int _nx) { this->nx = _nx; }
	void setny(int _ny) { this->ny = _ny; }
	void setnz(int _nz) { this->nz = _nz; }


	const double getWallXPosition() const { return  wallXPosition; }
	const double getdzdy() const { return dzdy; }
	const double getWallLength() const { return wallLength; }
	const double getLambda() const { return lambda; }
	const double getk() const { return k; }
	const double getTheta() const { return theta; }
	const double getPhi() const { return phi; }
	const double getPsi() const { return psi; }
	const int getl() const { return l; }
	const double geta() const { return a; } 
	const int getWallDivisions() const { return wallDivisions; }

	std::string getSaveString() { return saveString; }

	void setWallXPosition(double _wallXPosition) { this->wallXPosition = _wallXPosition; }
	void setdzdy(double _dzdy) { this->dzdy = _dzdy; }
	void setWallLength(double _wallLength) { this->wallLength = _wallLength; }
	void setLambda(double _lambda) { this->lambda = _lambda; }
	void setk(double _k) { this->k = _k; }
	void setTheta(double _theta) { this->theta = _theta; }
	void setPhi(double _phi) { this->phi = _phi; }
	void setPsi(double _psi) { this->psi = _psi; }
	void setl(int _l) { this->l = _l; }
	void seta(double _a) { this->a = _a; }


	void setWallDivisions(double _dzdy, double _wallLength) {this->wallDivisions =  static_cast<int>(std::floor(_wallLength / _dzdy));}

	void setCellName(std::string _setcellName) { this->cellName = _setcellName; }
	void setCellType(std::string _setcellType) { this->cellType = _setcellType; }
	void setCellCentering(std::string _cellCentering) { this->cellCentering = _cellCentering; }
	void setConfigTag(std::string _configTag) { this->configTag = _configTag; }
	
	void setSaveString(Config config);
	void setManualConfig();
	void setManualConfigTag(Config &config);
	void setFullConfig(double _wallXPosition, double _dzdy, double _wallLength, double _lambda, double _theta, double _phi, double _psi, int _l, double _a, int _nx, int _ny, int _nz, std::string _cellCentering, std::string _cellName, std::string _configTag);
	
	const std::string getCellName() const { return cellName; }
	const std::string getCellType() const { return cellType; }
	const std::string getCellCentering() const { return cellCentering; }
	const std::string getConfigTag() const { return configTag; }

private:
	int nx = 1;
	int ny = 1;
	int nz = 1;

	double wallXPosition = 0.01;// wall at x = wallXPos meters
	double dzdy = 0.0002;// dzdy->step size for wall spacing on a flat surface
	double wallLength = .1; // wall is wallLen x wallLen meters in area
	double lambda = 1e-9; // wavelength (meters)
	double k = 6.2832 / 3.e-9; // wavenumber (2 * pi / lambda)
	double theta = 0; // Rotation around __ axis
	double phi = 0; // Rotation around __ axis
	double psi = 0; // Rotation around __ axis
	int wallDivisions; // number of bins in the square (flat) wall
	int l = 3; // number of spherical harmonics to include (eventually will investigate this, but higher than l=1 or l=2 may )
	double a = 1.;// a refers to the atomic size

	// Defaults 
	std::string cellName = "salt";
	std::string cellType = "";
	std::string cellCentering = "face-centered";
	std::string saveString = "C:\\Users\\Michael\\Documents\\Programming\\laueDiffractionResults\\EigenResultsFaceHighRes.csv";
	std::string configTag = "_DefaultParameters_";
};

#endif 