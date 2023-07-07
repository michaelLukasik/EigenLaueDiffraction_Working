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


	//double get() { return; }
	const double getWallXPosition() const { return  wallXPosition; }
	const double getdzdy() const { return dzdy; }
	const double getWallLength() const { return wallLength; }
	const double getLambda() const { return lambda; }
	const double getA() const { return A; }
	const double getk() const { return k; }
	const double getOmega() const { return omega; }
	const double getTheta() const { return theta; }
	const double getPhi() const { return phi; }
	const double getPsi() const { return psi; }
	const int getWallDivisions() const { return wallDivisions; }

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
	void setWallDivisions(double _dzdy, double _wallLength) {this->wallDivisions =  static_cast<int>(std::floor(_wallLength / _dzdy));}

	void setCellName(std::string _setcellName) { this->cellName = _setcellName; }
	void setCellType(std::string _setcellType) { this->cellType = _setcellType; }
	void setCellCentering(std::string _cellCentering) { this->cellCentering = _cellCentering; }
	void setConfigTag(std::string _configTag) { this->configTag = _configTag; }
	
	void setSaveString(Config config);
	void setManualConfig();
	void setFullConfig(double _wallXPosition, double _dzdy, double _wallLength, double _lambda, double _A, double _k, double _omega, double _theta, double _phi, double _psi, std::string _cellCentering, std::string _cellName, std::string _configTag);
	
	const std::string getCellName() const { return cellName; }
	const std::string getCellType() const { return cellType; }
	const std::string getCellCentering() const { return cellCentering; }
	const std::string getConfigTag() const { return configTag; }

private:
	int nx = 1;
	int ny = 1;
	int nz = 1;

	double wallXPosition = 0.01;// wall at x = wallXPos
	double dzdy = 0.0002;// dzdy->step size for wall spacing
	double wallLength = .1; // wall is wallLen x wallLen meters
	double lambda = 1e-9; // wavelength (meters)
	double A = 1;  // Amplitude (arbitrary units)
	double k = 1 ; // wavenumber (1/m)
	double omega = 0; //
	double theta = 0;
	double phi = 0;
	double psi = 0;
	int wallDivisions; //= static_cast<int>(std::floor(wallLength / dzdy));


	// Defaults 
	std::string cellName = "salt";
	std::string cellType = "";
	std::string cellCentering = "face-centered";
	std::string saveString = "C:\\Users\\Michael\\Documents\\Programming\\laueDiffractionResults\\EigenResultsFaceHighRes.csv";
	std::string configTag = "_DefaultParameters_";
};

#endif 