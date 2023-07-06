#include <Diffraction.h>
#include <Bravais.h>
#include <Crystal.h>
#include <config.h>
#include <math.h>
#include <Eigen\Dense>

//setManualConfig is the easiest way to produce one image for testing purposes
void Config::setManualConfig() {
	Config::setWallXPosition(0.01);
	Config::setdzdy(0.0001);
	Config::setWallLength(.1);
	Config::setLambda(1.e-9);
	Config::setA(1.);
	Config::setk(1.);
	Config::setOmega(0.);
	Config::setTheta(0.);
	Config::setPhi(0.);
	Config::setPsi(0.);

	Config::setCellCentering("face-centered");
	Config::setCellName("salt");
	Config::setCellType("");
	
}

//setFullConfig is meant to be used to automate large numbers of pictures at different distances, spacings, wavelengths etc... 
void Config::setFullConfig(double wallXPosition, double dzdy, double wallLength, double lambda, double A, double k, double omega, double theta, double phi, double psi , std::string cellCentering , std::string cellName , std::string cellType) {
	Config::setWallXPosition(wallXPosition);
	Config::setdzdy(dzdy);
	Config::setWallLength(wallLength);
	Config::setLambda(lambda);
	Config::setA(A);
	Config::setk(k);
	Config::setOmega(omega);
	Config::setTheta(theta);
	Config::setPhi(phi);
	Config::setPsi(psi);

	Config::setCellCentering(cellCentering);
	Config::setCellName(cellName);
	Config::setCellType(cellType);
}


