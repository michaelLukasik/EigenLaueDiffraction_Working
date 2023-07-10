#include <Diffraction.h>
#include <Bravais.h>
#include <Crystal.h>
#include <config.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <iterator>
#include <iostream>
#include <complex>
#include <numeric>
#include <Eigen\Dense>
#include <chrono>


// Create a Screen to project the transmition onto


Eigen::MatrixXd getScreen(const Config& config) {
	Eigen::MatrixXd screen(static_cast<int>(std::floor(std::pow(config.getWallDivisions(), 2))), 3);
	for (int i = 0; i < std::pow(config.getWallDivisions(), 2); ++i) {
		screen(i, 0) = config.getWallXPosition();
		screen(i, 1) = (std::floor(i / config.getWallDivisions())) * config.getdzdy() - (config.getWallLength() / 2.0); // Second term to center the screen on 0y
		screen(i, 2) = (i % config.getWallDivisions()) * config.getdzdy() - (config.getWallLength() / 2.0);				// Second term to center the screen on 0z
	}
	return screen	;
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
bool isAlmostEqual(const std::vector<double>& v1, const std::vector<double>& v2) {
	return std::sqrt(pow(v2[0] - v1[0], 2) + pow(v2[1] - v1[1], 2) + pow(v2[2] - v1[2], 2)) < 0.10*std::sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
}
Eigen::MatrixXd buildLattice(const Crystal& crystal, const Config& config, const Eigen::MatrixXd& cellStructure){
//Eigen::MatrixXd buildLattice(int xCells, int yCells, int zCells, double axialDistanceA, double axialDistanceB, double axialDistanceC, const Eigen::MatrixXd& cellStructure) {
	std::vector<std::vector<double>> fullLatticeVectors = {}; // We use a standard container here so we can use pushback (and later unique). Eigen::Matrix doesnt have a pushback feature due to memory concerns
	for (int a = 0; a < cellStructure.rows(); a++) {
		for (int nxc = 0; nxc <= config.getnx(); nxc++) {
			for (int nyc = 0; nyc <= config.getny(); nyc++) {
				for (int nzc = 0; nzc <= config.getnz(); nzc++) {
					Eigen::Vector3d latticePointPosition(cellStructure(a, 0), cellStructure(a, 1), cellStructure(a, 2));
					Eigen::Vector3d offsetVector(crystal.getAxialDistanceA() * nxc, crystal.getAxialDistanceB() * nyc, crystal.getAxialDistanceC() * nzc);
					Eigen::Vector3d centeringVector((config.getnx() + 1. / 2.) * crystal.getAxialDistanceA(), (config.getny() + 1. / 2.) * crystal.getAxialDistanceB(), (config.getnz() + 1. / 2.) * crystal.getAxialDistanceC());
					Eigen::Vector3d finalLatticePointPosition = latticePointPosition + offsetVector - centeringVector;
					std::vector<double> finalPositionVector(finalLatticePointPosition.data(), finalLatticePointPosition.data() + finalLatticePointPosition.size());
					fullLatticeVectors.push_back(finalPositionVector);

				}
			} 
		}
	}
	std::sort(fullLatticeVectors.begin(), fullLatticeVectors.end());
	fullLatticeVectors.erase(std::unique(fullLatticeVectors.begin(), fullLatticeVectors.end(), isAlmostEqual), fullLatticeVectors.end()); // Remove duplicates in vector form
	
	std::cout << std::endl << "There are " << fullLatticeVectors.size() << " sites in the final lattice after culling duplicate points" << std::endl;
	
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
void propogateWave(double lambda, const Eigen::VectorXd& normsToScreen, const Eigen::MatrixXd& screenPositions, Eigen::MatrixXd& waveFunction) {
	Eigen::VectorXd numerator = (std::atan(1) * 2. / lambda) * normsToScreen;
	Eigen::VectorXd cosResult = numerator.array().cos() / normsToScreen.array();
	Eigen::VectorXd sinResult = numerator.array().sin() / normsToScreen.array();
	waveFunction.col(0) += cosResult;
	waveFunction.col(1) += sinResult;
}

void exportData(Eigen::MatrixXd waveFunction, const Config& config) {
	std::vector<std::vector<double>> waveFunctionVectors = {};
	for (int i = 0; i < waveFunction.rows(); i++) {
		Eigen::VectorXd waveFunctionAtScreenPosition(5);
		waveFunctionAtScreenPosition << waveFunction(i, 0), waveFunction(i, 1), waveFunction(i, 2), waveFunction(i, 3), waveFunction(i, 4);
		std::vector<double> finalWaveFunctionVector(waveFunctionAtScreenPosition.data(), waveFunctionAtScreenPosition.data() + waveFunctionAtScreenPosition.size());
		waveFunctionVectors.push_back(finalWaveFunctionVector);
	}
	std::string saveTag = config.getCellCentering() + "_" + config.getCellName() + "_" + config.getCellType() + "_" + std::to_string(config.getWallLength())
													+ "Length_" + std::to_string(config.getWallDivisions()) + "Divisions_" + config.getConfigTag();
	std::string savePath = "C:\\Users\\Michael\\Documents\\Programming\\laueDiffractionResults\\EigenResults_" + saveTag + ".csv";
	std::ofstream out(savePath);
	for (auto& row : waveFunctionVectors) {
		for (auto col : row)
			out << col << ',';
		out << '\n';
	}
}

//int main() {
  int main(int argc, char* argv[]) {

	auto timeStart = std::chrono::high_resolution_clock::now();

	// Configs are now defined within config.cpp
	Config config;
	if (argc == 14) {
		config.setFullConfig(std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]),
			std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]),
			std::atof(argv[9]), std::atof(argv[10]), argv[11], argv[12], argv[13]);
	}
	else {
		std::cout << argc << std::endl;
		std::cout << "A full list of argument is 13 arugments long and follows the order: \n"
			" [1] double wallXPosition, [2] double dzdy, [3] double wallLength,\n"
			" [4] double lambda, [5] double A, [6] double k,\n"
			" [7] double omega, [8] double theta, [9] double phi,\n" 
			" [10] double psi, [11] string cellName, [12] string cellType "
			" [12] string cellCentering, [13] string configTag"<< std::endl;
		config.setManualConfig();
	}
	// Get properties of individual crystal with the proper ortientation
	Crystal crystal; 
	crystal.setCellStrings(config);
	crystal.setCellType(crystal);
	crystal.setCellProperties(crystal);

	std::cout << "Using cellName = " << config.getCellName() << " || Using cellType = " << crystal.getCellType() << " || Using cellCentering = " << config.getCellCentering() << std::endl;
	
	Eigen::MatrixXd cellStructure = getCellStructure(crystal);
	
	// Build up the lattice and screen
	Eigen::MatrixXd fullLattice = buildLattice(crystal, config, cellStructure);
	Eigen::MatrixXd rotatedLattice = rotateLattice(fullLattice, config);
	Eigen::MatrixXd wave = buildWave(getScreen(config), config);
	
	auto timeToPropogation = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - timeStart);
	std::cout << "Time to propogation loop: " << timeToPropogation.count() << " ms" << std::endl;

	// For each atom, get the contribution from scattering to every screen position. This is the bulk of the calculation 
	for (int atomIndex = 0 ; atomIndex < rotatedLattice.rows(); atomIndex++) {
		auto timeFromLoopStart = std::chrono::high_resolution_clock::now();
		std::cout << "Working on atom number " << atomIndex;
		Eigen::VectorXd normsToScreen = getNormToScreenPosition(getScreen(config), rotatedLattice, atomIndex);
		std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - timeFromLoopStart).count() << " ms to build norm Vector) ";
		propogateWave(config.getLambda(), normsToScreen, getScreen(config), wave);
		std::cout << " (" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - timeFromLoopStart).count() << " ms Total)" << std::endl;
	}
	
	
	auto timeToPropogationEnd = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - timeStart);
	std::cout << "Total time:  " << timeToPropogationEnd.count() / 1000. << " s" << std::endl;
	std::cout << "Time per atom:  ~" << (timeToPropogationEnd.count() - timeToPropogation.count())/ (1000. * static_cast<float>(rotatedLattice.rows())) << " s" << std::endl;

	exportData(wave, config);


	return 0;
}	

