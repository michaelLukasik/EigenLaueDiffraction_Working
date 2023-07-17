#include <Diffraction.h>
#include <Bravais.h>
#include <Crystal.h>
#include <config.h>
#include <fstream>
#include <iostream> 
#include <math.h>
#include <Eigen\Dense>

//setManualConfig is the easiest way to produce one image for testing purposes
void Config::setManualConfigTag(Config &config) {
	std::stringstream configTag;
	double xpos_str = config.getWallXPosition();
	double dzdy_str = config.getdzdy();
	double len_str = config.getWallLength();
	double lambda_str = config.getLambda();
	double theta_str = config.getTheta();
	double psi_str = config.getPsi();
	double phi_str = config.getPhi();
	double nx_str = config.getnx();
	double ny_str = config.getny();
	double nz_str = config.getnz();

	std::string name_str = config.getCellName();
	std::string centering_str = config.getCellCentering();
	 
	// Combine them all in a similar way s.t. we can use the plotting script for the rotating case
	
	std::fixed;
	configTag << "CONFIG_" << xpos_str << "XPOS_" << dzdy_str << "DZDY_" << len_str << "LEN_" << theta_str << "THETA_" << phi_str << "PHI_" << psi_str << "PSI_" << nx_str << "NX_" << ny_str << "NY_" << nz_str << "NZ_" << lambda_str << "LAMBDA_";
	std::string fullConfigTag = configTag.str() + centering_str + "CENTERING_" + name_str + "NAME_CONFIG_";
	config.setConfigTag(fullConfigTag);
}

void Config::setManualConfig() {

	Config::setWallXPosition(1.e-2);
	Config::setdzdy(1.e-3);
	Config::setWallLength(1.e-1);
	Config::setLambda(1.e-10);
	Config::setk(8 * std::atan(1.) / Config::getLambda()); //  (2*pi) /  lambda
	Config::setTheta(0.);
	Config::setPhi(0.);
	Config::setPsi(0.);
	Config::setl(3.); // Number of l-terms in ylm to include (drops off rather hard after ~2)
	Config::seta(1.e-10); //1.e-10 // "Atomic size" used for scattering calculation, set to about 1 angstrom  ** ENFORCE A LIMIT that this must be small, or generalize with another flag if we want to use this for macroscpoic objects
	Config::setnx(1);  // number of copies of cells in x,y,z 
	Config::setny(1);
	Config::setnz(1);
	Config::setWallDivisions(Config::getdzdy(), Config::getWallLength());

	//Config::setSaveString("C:\\Users\\Michael\\Documents\\Programming\\laueDiffractionResults\\rotationTests\\EigenResultsFaceHighRes.csv");
	Config::setCellCentering("primative");
	Config::setCellName("salt");
	Config::setCellType("");
	Config::setConfigTag("_ShouldNotAppear_");
}

//setFullConfig is meant to be used to automate large numbers of pictures at different distances, spacings, wavelengths etc... 
void Config::setFullConfig(double wallXPosition, double dzdy, double wallLength, double lambda, double theta, double phi, double psi, int l, double a, int nx, int ny, int nz, std::string cellCentering, std::string cellName, std::string configTag) { //std::string cellType, std::string configTag) {
	Config::setWallXPosition(wallXPosition);
	Config::setdzdy(dzdy);
	Config::setWallLength(wallLength);
	Config::setLambda(lambda);
	Config::setk(8 * std::atan(1.) / Config::getLambda());
	Config::setTheta(theta);
	Config::setPhi(phi);
	Config::setPsi(psi);
	Config::setl(l);
	Config::seta(a);
	Config::setWallDivisions(dzdy, wallLength);
	Config::setnx(static_cast<int>(nx));
	Config::setny(static_cast<int>(ny));
	Config::setnz(static_cast<int>(nz));

	Config::setCellCentering(cellCentering);
	Config::setCellName(cellName);
	Config::setConfigTag(configTag);
}



