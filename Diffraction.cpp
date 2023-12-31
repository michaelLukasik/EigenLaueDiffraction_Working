#include <Diffraction.h>
#include <Bravais.h>
#include <Crystal.h>
#include <config.h>
#include <math.h>
#include <windows.h>
#include <ppl.h>
#include <cmath>
#include <ctime>
#include <complex>
#include <fstream>
#include <vector>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <complex>
#include <numeric>
#include <Eigen\Dense>
#include <chrono>
#include <execution>
#include <algorithm>
#include <Eigen/../unsupported/Eigen/FFT>


std::string returnTime()
{
	auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	return std::ctime(&t);

}
// Create a Screen to project the transmition onto
Eigen::MatrixXd getScreen(const Config& config) {
	Eigen::MatrixXd screen(static_cast<int>(std::floor(std::pow(config.getWallDivisions(), 2))), 3);
	for (int i = 0; i < std::pow(config.getWallDivisions(), 2); ++i) {
		screen(i, 0) = config.getWallXPosition();
		screen(i, 1) = (std::floor(i / config.getWallDivisions())) * config.getdzdy() - (config.getWallLength() / 2.0); // Second term to center the screen on 0y
		screen(i, 2) = (i % config.getWallDivisions()) * config.getdzdy() - (config.getWallLength() / 2.0);				// Second term to center the screen on 0z
	}
	return screen;
}

Eigen::MatrixXd buildWave(const Eigen::MatrixXd& screen, const Config& config){
	Eigen::MatrixXd waveFunction(static_cast<int>(std::pow(config.getWallDivisions(), 2)), 5);
	for (int i = 0; i < screen.rows(); i++) {
		waveFunction(i, 0) = 0; // Real part of wavefunction
		waveFunction(i, 1) = 0; // Imaginary part of wavefunction
		waveFunction(i, 2) = screen(i, 0); // X
		waveFunction(i, 3) = screen(i, 1); // Y
		waveFunction(i, 4) = screen(i, 2); // Z
	}
	return waveFunction;
}
bool isAlmostEqual(const std::vector<double>& v1, const std::vector<double>& v2, const double atomicDistance) {
	std::cout << ((v2[0] - v1[0]) < 0.1 * atomicDistance and (v2[1] - v1[1]) < 0.1 * atomicDistance and (v2[2] - v1[2]) < 0.1 * atomicDistance) << std::endl;
	return ((v2[0] - v1[0]) < 0.1 * atomicDistance and (v2[1] - v1[1]) < 0.1 * atomicDistance and (v2[2] - v1[2]) < 0.1 * atomicDistance);
}
Eigen::MatrixXd buildLattice(const Crystal& crystal, const Config& config, const Eigen::MatrixXd& cellStructure){
	double atomicDistance = config.geta();
	std::vector<std::vector<double>> fullLatticeVectors = {}; // We use a standard container here so we can use pushback (and later unique). Eigen::Matrix doesnt have a pushback feature due to memory concerns
	for (int a = 0; a < cellStructure.rows(); a++) {
		for (int nxc = 0; nxc <= config.getnx(); nxc++) {
			for (int nyc = 0; nyc <= config.getny(); nyc++) {
				for (int nzc = 0; nzc <= config.getnz(); nzc++) {
					Eigen::Vector3d latticePointPosition(cellStructure(a, 0), cellStructure(a, 1), cellStructure(a, 2));
					Eigen::Vector3d offsetVector(crystal.getAxialDistanceA() * nxc, crystal.getAxialDistanceB() * nyc, crystal.getAxialDistanceC() * nzc);
					Eigen::Vector3d centeringVector(( (config.getnx() + 1.) / 2.) * crystal.getAxialDistanceA(), ( (config.getny() + 1.) / 2.) * crystal.getAxialDistanceB(), ( (config.getnz() + 1. ) / 2.) * crystal.getAxialDistanceC());
					Eigen::Vector3d finalLatticePointPosition = latticePointPosition + offsetVector - centeringVector;
					std::vector<double> finalPositionVector(finalLatticePointPosition.data(), finalLatticePointPosition.data() + finalLatticePointPosition.size());
					fullLatticeVectors.push_back(finalPositionVector);

				}
			} 
		}
	}
	std::sort(fullLatticeVectors.begin(), fullLatticeVectors.end());
	fullLatticeVectors.erase(std::unique(fullLatticeVectors.begin(), fullLatticeVectors.end(), [&atomicDistance](const std::vector<double>& v1, const std::vector<double>& v2)
		{return (std::abs(v2[0] - v1[0]) < atomicDistance and std::abs(v2[1] - v1[1]) <  atomicDistance and std::abs(v2[2] - v1[2]) < atomicDistance);}), fullLatticeVectors.end()); // Remove duplicates in vector form

	std::cout << std::endl << "There are " << fullLatticeVectors.size() << " sites in the final lattice after culling duplicate points \n" << std::endl;
	
	Eigen::MatrixXd fullLattice(fullLatticeVectors.size(), 3);

	for (int i = 0; i < fullLatticeVectors.size(); i++) {
		fullLattice.row(i) << fullLatticeVectors[i][0], fullLatticeVectors[i][1], fullLatticeVectors[i][2]; //put the vectors back into Eigen form for easy use later
	}
	return fullLattice;
}
Eigen::MatrixXd rotateLattice(const Eigen::MatrixXd& fullLattice, const Config& config) {
	Eigen::Matrix3d rotmatPhi; 
	Eigen::Matrix3d rotmatTheta;
	Eigen::Matrix3d rotmatPsi;

	double phi = atan(1.) * 4. / 180. * config.getPhi();
	double psi = atan(1.) * 4. / 180. * config.getPsi();
	double theta = atan(1.) * 4. / 180. * config.getTheta();


	rotmatPhi << std::cos(phi), 0., std::sin(phi), 0., 1., 0., -std::sin(phi), 0, std::cos(phi);
	rotmatTheta << std::cos(theta), -std::sin(theta), 0., std::sin(theta), std::cos(theta),0., 0.,0.,1. ;
	rotmatPsi << 1., 0., 0., 0., std::cos(psi), -std::sin(psi), 0., std::sin(psi), std::cos(psi);
	return fullLattice * rotmatPhi * rotmatTheta * rotmatPsi;
}

Eigen::VectorXd getNormToScreenPosition(const Eigen::MatrixXd& screenPositions, const Eigen::MatrixXd& rotatedLattice, int atomIndex) {
	Eigen::VectorXd normToScreenPosition(screenPositions.rows());	
	normToScreenPosition =	 pow(screenPositions.col(0).array() - rotatedLattice(atomIndex, 0),2)
						 +   pow(screenPositions.col(1).array() - rotatedLattice(atomIndex, 1), 2)
						 +   pow(screenPositions.col(2).array() - rotatedLattice(atomIndex, 2), 2);

	return normToScreenPosition.array().sqrt();
	
}

Eigen::VectorXcd sphBesselFirstKind(int l, Eigen::VectorXcd x) {
	//if (l == 0) { return ;}
	if (l == 0) { return x.array().sin() / x.array(); }
	if (l == 1) { return (x.array().sin() / pow(x.array(),2)) - (x.array().cos()) / x.array(); }
	if (l == 2) { return ((3. / x.array().square()) - 1.) * (x.array().sin() / x.array()) - (3. * x.array().cos() / x.array().square()); }
	if (l == 3) { return ((15. / x.array().cube()) - (6. / x.array())) * (x.array().sin() / x.array()) - ((15. / x.array().square() -1)*(x.array().cos() / x.array())); }
	if (l >= 4) { throw "Currently not taking l > 3"; } // Computation/results payoff not worth it. Could be interesting for l-spectrum study, but currently not needed. 

}

double sphBesselFirstKind(int l, double  x) {
	//if (l == 0) { return ;}
	if (l == 0) { return std::sin(x) / x; }
	else if (l == 1) { return (std::sin(x) / pow(x, 2)) - (std::cos(x) / x); }
	else if (l == 2) { return ((3. / pow(x, 2)) - 1.) * (std::sin(x) / x) - (3. * std::cos(x) / pow(x, 2)); }
	else if (l == 3) { return ((15. / pow(x, 3)) - (6. / x)) * (std::sin(x) / x) - ((15. / pow(x, 2)) - 1.) * (std::cos(x) / x); }
	else if (l >= 4) { throw "Currently not taking l > 3"; }

}

Eigen::VectorXcd sphBesselSecondKind(int l, Eigen::VectorXcd x) {
	//if (l == 0) { return ;}
	if (l == 0) { return -1.*x.array().cos() / x.array(); }
	else if (l == 1) { return (x.array().sin() / x.array().square()) - (x.array().sin() / x.array()); }
	else if (l == 2) { return ((-3./ x.array().square()) + 1.) * (x.array().cos() / x.array()) - (3. * x.array().sin() / pow(x.array(), 2)); }
	else if (l == 3) { return ((-15./ x.array().cube()) + (6. / x.array())) * (x.array().cos() / x.array()) - ((15. / x.array().square() - 1) * (x.array().sin() / x.array())); }
	else if (l >= 4) { throw "Currently not taking l > 3"; }
}

double  sphBesselSecondKind(int l, double x) {
	//if (l == 0) { return ;}
	if (l == 0) { return -std::cos(x) / x; }
	else if (l == 1) { return (std::sin(x) / pow(x, 2)) - (std::sin(x) / x); }
	else if (l == 2) { return ((-3. / pow(x, 2)) + 1.) * (std::cos(x) / x) - (3. * std::sin(x) / pow(x, 2)); }
	else if (l == 3) { return ((-15. / pow(x, 3)) + (6. / x)) * (std::cos(x) / x) - ((15. / pow(x, 2)) - 1.) * (std::sin(x) / x); }
	else if (l >= 4) { throw "Currently not taking l > 3"; }

}


Eigen::VectorXcd legandrePoly(int l, Eigen::VectorXcd x) {

	if (l == 0) { return pow(x.array(), 0); }
	else if (l == 1) { return x; }
	else if (l == 2) { return (1. / 2.) * (3. * x.array().square() - 1.); }
	else if (l == 3) { return (1. / 2.) * (5. * pow(x.array(), 3.) - 3. * x.array()); }
	//if (l == 4) { return (1. / 8.) * (35. * pow(x.array(), 4) - 30. * pow(x.array(),2) + 3); }
	else if (l >= 4) { throw "Currently not taking l > 3"; }

}

Eigen::VectorXcd hankelFirstKind(int l , Eigen::VectorXcd x) {
	std::complex<double> i(0.0, 1.0);
	return sphBesselFirstKind(l, x) + i * sphBesselSecondKind(l, x);
}

std::complex<double> hankelFirstKind(int l, double x) {
	std::complex<double> i(0.0, 1.0);
	return sphBesselFirstKind(l, x) + i * sphBesselSecondKind(l, x);
}

std::complex<double> scatteringAmplitude(const Config& config, const int rows, int l) {

	std::complex<double> kc(config.getk(), 0.0);
	long double k = config.getk();
	long double a = config.geta();
	long double ka = k * a;
	Eigen::VectorXcd iVector =Eigen::VectorXcd::Constant(rows,(1.0,0.0));
	std::complex<double> i(0.0, 1.0);
	double besselka = sphBesselFirstKind(l, ka);
	std::complex<double> numerator = i * besselka;
	std::complex<double> denominator = hankelFirstKind(l, ka) * kc;
	return numerator / denominator; 

}
void propogateWave(const Config& config, const Eigen::VectorXd& normsToScreen, const Eigen::MatrixXd& screenPositions, Eigen::MatrixXd& waveFunction, const Eigen::MatrixXd& rotatedLattice, int atomIndex) { // https://physics.stackexchange.com/questions/633609/understanding-quantum-hard-sphere-scattering

	const std::complex<double> i(0.0, 1.0);
	const std::complex<double> r(1.0, 0.0);

	Eigen::VectorXd adjacentDistanceForAngleFromAtom(screenPositions.rows());
	Eigen::VectorXd hypotenuseDistanceForAngleFromAtom(screenPositions.rows());

	adjacentDistanceForAngleFromAtom = screenPositions.col(0).array() - rotatedLattice(atomIndex, 0);
	hypotenuseDistanceForAngleFromAtom = pow(screenPositions.col(1).array() - rotatedLattice(atomIndex, 1), 2)
										+ pow(screenPositions.col(2).array() - rotatedLattice(atomIndex, 2), 2);
	hypotenuseDistanceForAngleFromAtom = hypotenuseDistanceForAngleFromAtom.array().sqrt();
	Eigen::VectorXd angleFromAtom = adjacentDistanceForAngleFromAtom.array() / hypotenuseDistanceForAngleFromAtom.array();

	Eigen::VectorXcd angleFromAtomCos = angleFromAtom.array().cos()*r;


	for (int l = 0; l <= config.getl(); l++) {
		const int nRows = screenPositions.rows();
		Eigen::VectorXcd kr = config.getk() * normsToScreen;
		Eigen::VectorXcd innerTerms = sphBesselFirstKind(l,kr)  + i * scatteringAmplitude(config, screenPositions.rows(),l) * hankelFirstKind(l, kr);
		Eigen::VectorXcd lPoly = legandrePoly(l, angleFromAtomCos);
		Eigen::VectorXcd partialWaveFunction = pow(i, l) * ((2. * l) + 1.) * innerTerms.array() * lPoly.array();
		waveFunction.col(0) += partialWaveFunction.real();
		waveFunction.col(1) += partialWaveFunction.imag();	
	}
}
//void propogateWaveFFT
void exportData(Eigen::MatrixXd waveFunction, const Config& config) {
	std::vector<std::vector<double>> waveFunctionVectors = {};
	for (int i = 0; i < waveFunction.rows(); i++) {
		Eigen::VectorXd waveFunctionAtScreenPosition(5);
		waveFunctionAtScreenPosition << waveFunction(i, 0), waveFunction(i, 1), waveFunction(i, 2), waveFunction(i, 3), waveFunction(i, 4);
		std::vector<double> finalWaveFunctionVector(waveFunctionAtScreenPosition.data(), waveFunctionAtScreenPosition.data() + waveFunctionAtScreenPosition.size());
		waveFunctionVectors.push_back(finalWaveFunctionVector);
	}

	std::string saveTag = config.getConfigTag();
	std::string savePath = "C:\\Users\\Michael\\Documents\\Programming\\laueDiffractionResults\\csvFiles\\EigenResults_" +saveTag + ".csv";
	
	std::ofstream out(savePath);
	for (auto &row : waveFunctionVectors) {
		for (auto col : row)
			out << col << ',';
		out << '\n';
	}
}


int main(const int argc, char* argv[]) {

	auto timeStart = std::chrono::high_resolution_clock::now();
	auto timeToPropogation = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - timeStart);
	const double approxTimePerCalc = 6.275e-6;
	bool waitForResult = 0; //  Set to true to see timing data for the propogation loop, will include in csv later
	
	// Configs are now defined within config.cpp
	Config config;
	if (argc == 16) {
		const char* arglist[16] = {"executable", "double wallXPosition", "double dzdy", "double wallLength,", "double lambda", "double theta",
											"double phi", "double psi", "int l (ylm's for harmonics)", "double a (atomic size)",
											"int nx", "int ny" , "int nz" ,"string cellCentering", "string cellName", "string configTag" };
		std::cout << "Paramets inputted: \n";
		for (int i = 0; i < argc; i++) {
			std::cout << i << ": " << argv[i] << " -- " << arglist[i] << std::endl;
		}
		config.setFullConfig(std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3]), 
			std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6]),
			std::atof(argv[7]), std::atof(argv[8]), std::atof(argv[9]), 
			std::atoi(argv[10]), std::atoi(argv[11]), std::atoi(argv[12]),
			argv[13], argv[14], argv[15] );
		waitForResult = 0;


	}
	else {
		std::cout << "Number of Arguments inputted: " << argc << std::endl;
		std::cout << "A full list of argument is 15 arugments long and follows the order: \n"
			" [1] double wallXPosition, [2] double dzdy, [3] double wallLength,\n"
			" [4] double lambda, [5] double theta, [6] double phi,\n"
			" [7] double psi, [8] int l (ylm), [9] double a (atomic length),\n"
			" [10] int nx, [11] int ny, [12] int nz, \n"
			" [13] string cellCentering, [14] string cellName, [15] string configTag," << std::endl;
		std::cout << "Using the parameters as defined in setManualConfig in config.cpp" << std::endl;
		config.setManualConfig();
		config.setManualConfigTag(config);

	}
	// Get properties of individual crystal with the proper ortientation
	Crystal crystal;
	crystal.setCellStrings(config);
	crystal.setCellType(crystal);
	crystal.setCellProperties(crystal);
	Eigen::MatrixXd cellStructure = getCellStructure(crystal);
	
	// Build up the lattice and screen
	Eigen::MatrixXd fullLattice = buildLattice(crystal, config, cellStructure);
	Eigen::MatrixXd rotatedLattice = rotateLattice(fullLattice, config);
	Eigen::MatrixXd wave = buildWave(getScreen(config), config);


	std::cout << "Using cellName = " << config.getCellName() << " || Using cellType = " << crystal.getCellType() << " || Using cellCentering = " << config.getCellCentering() << std::endl;
	std::cout << pow(config.getWallLength() / config.getdzdy(), 2) << " points to calculate for " << rotatedLattice.rows() << " atoms in the lattice \n";
	std::cout << "\n ================================================== \n";
	std::cout << "Start time = " << returnTime();
	std::cout << "CSV will be about  " << std::setprecision(3) << 0.0478 * pow(config.getWallLength() / config.getdzdy(), 2) / 1000. << "MBs\n";
	std::cout << "Come back in *roughly* " << std::setprecision(3) << approxTimePerCalc * pow(config.getWallLength() / config.getdzdy(), 2) * rotatedLattice.rows() / 60. << " minutes" << std::endl;

	
	for(int aIndex = 0; aIndex < rotatedLattice.rows(); aIndex++){
		auto timeFromLoopStart = std::chrono::high_resolution_clock::now();
		Eigen::VectorXd normsToScreen = getNormToScreenPosition(getScreen(config), rotatedLattice, aIndex);
		propogateWave(config, normsToScreen, getScreen(config), wave, rotatedLattice, aIndex);
		};

	
	
	auto timeToPropogationEnd = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - timeStart);
	std::cout << "Total time:  " << timeToPropogationEnd.count() / 1000. << " s" << std::endl;
	std::cout << "Time per atom:  ~" << (timeToPropogationEnd.count() - timeToPropogation.count())/ (1000. * static_cast<float>(rotatedLattice.rows())) << " s" << std::endl;

	exportData(wave, config);
	std::cout << "Time end: " << returnTime();
	if (waitForResult == 1 ) { std::cin.ignore(); }
	return 0;
}	

