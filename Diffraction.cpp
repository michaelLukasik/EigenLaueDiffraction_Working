#include <Diffraction.h>
#include <Bravais.h>
#include <Crystal.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <math.h>
#include <iterator>
#include <iostream>
#include <complex>
#include <numeric>
#include <Eigen\Dense>

// Create a Screen to project the transmition onto

Eigen::MatrixXd getScreen(){
	Eigen::MatrixXd screen(static_cast<int>(std::floor(std::pow(wallDivisions, 2))), 3);
	for (int i = 0; i < std::pow(wallDivisions,2); ++i) {
		screen(i, 0) = wallXPosition;
		screen(i, 1) = (std::floor(i / wallDivisions)) * dzdy - (wallLength / 2.0); // Second term to center the screen on 0y
		screen(i, 2) = (i % wallDivisions) * dzdy - (wallLength / 2.0);				// Second term to center the screen on 0z
		if (debug == true) {
			std::cout << screen(i, 0) << "  " << screen(i, 1) << "  " << screen(i, 2) << " \n";
		}
	}
	return screen	;
}

Eigen::MatrixXd buildWave(Eigen::MatrixXd screen){
	Eigen::MatrixXd waveFunction(static_cast<int>(std::pow(wallDivisions, 2)), 5);
	for (int i = 0; i < screen.rows(); i++) {
		waveFunction(i, 0) = 0; // Real part of wavefunction
		waveFunction(i, 1) = 0; // Imaginary part of wavefunction
		waveFunction(i, 2) = screen(i, 0); // X
		waveFunction(i, 3) = screen(i, 1); // Y
		waveFunction(i, 4) = screen(i, 2); // Z
		if (debug) {
			std::cout << "Vector " << i << " looks like : \n";
			std::cout << waveFunction(i, 0) << " " << waveFunction(i, 1) << " " << waveFunction(i, 2) << " " << waveFunction(i, 3) << " " << waveFunction(i, 4) << "\n";
		}
	}
	return waveFunction;
}
bool isAlmostEqual(std::vector<double> v1, std::vector<double> v2) {
	return std::sqrt(pow(v2[0] - v1[0], 2) + pow(v2[1] - v1[1], 2) + pow(v2[2] - v1[2], 2)) < 0.10*std::sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
}
Eigen::MatrixXd buildLattice(int xCells, int yCells, int zCells, double axialDistanceA, double axialDistanceB, double axialDistanceC, Eigen::MatrixXd cellStructure) {
	std::vector<std::vector<double>> fullLatticeVectors = {}; // We use a standard container here so we can use pushback (and later unique). Eigen::Matrix doesnt have a pushback feature due to memory concerns
	for (int a = 0; a < cellStructure.rows(); a++) {
		for (int nxc = 0; nxc <= xCells; nxc++) {
			for (int nyc = 0; nyc <= yCells; nyc++) {
				for (int nzc = 0; nzc <= zCells; nzc++) {
					Eigen::Vector3d latticePointPosition(cellStructure(a, 0), cellStructure(a, 1), cellStructure(a, 2));
					Eigen::Vector3d offsetVector(axialDistanceA * nxc, axialDistanceB * nyc, axialDistanceC * nzc);
					Eigen::Vector3d centeringVector((xCells + 1. / 2.) * axialDistanceA, (yCells + 1. / 2.) * axialDistanceB, (zCells + 1. / 2.) * axialDistanceC);
					Eigen::Vector3d finalLatticePointPosition = latticePointPosition + offsetVector - centeringVector;
					std::vector<double> finalPositionVector(finalLatticePointPosition.data(), finalLatticePointPosition.data() + finalLatticePointPosition.size());
					fullLatticeVectors.push_back(finalPositionVector);

				}
			} //[](double l, double r) { return std::abs(l - r) < 0.01; }
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

Eigen::MatrixXd rotateLattice(Eigen::MatrixXd fullLattice,double phideg,double thetadeg, double psideg) {
	Eigen::Matrix3d rotmatPhi; 
	Eigen::Matrix3d rotmatTheta;
	Eigen::Matrix3d rotmatPsi;

	rotmatPhi << std::cos(phi), 0., std::sin(phi), 0., 1., 0., -std::sin(phi), 0, std::cos(phi);
	rotmatTheta << std::cos(theta), -std::sin(theta), 0., std::sin(theta), std::cos(theta),0., 0.,0.,1. ;
	rotmatPsi << 1., 0., 0., 0., std::cos(psi), -std::sin(psi), 0., std::sin(psi), std::cos(psi);
	return fullLattice * rotmatPhi * rotmatTheta * rotmatPsi;
}

//int getNormToScreenPosition(Eigen::MatrixXd screenPositions, Eigen::MatrixXd rotatedLattice, int atomIndex) {
Eigen::VectorXd getNormToScreenPosition(Eigen::MatrixXd screenPositions, Eigen::MatrixXd rotatedLattice, int atomIndex) {
	
	Eigen::VectorXd normToScreenPosition(screenPositions.rows());
	for (int i = 0; i < screenPositions.rows(); i++) {
		normToScreenPosition(i) = std::sqrt(  pow(screenPositions(i, 0) - rotatedLattice(atomIndex, 0), 2)
											+ pow(screenPositions(i, 1) - rotatedLattice(atomIndex, 1), 2)
											+ pow(screenPositions(i, 2) - rotatedLattice(atomIndex, 2), 2));
	}
	return normToScreenPosition;
}
Eigen::MatrixXd propogateWave(double lambda, Eigen::VectorXd normsToScreen, Eigen::MatrixXd screenPositions, Eigen::MatrixXd waveFunction) {
	for (int i = 0; i < normsToScreen.rows(); i++) {
		waveFunction(i,0) += std::cos(std::atan(1) * 2. * normsToScreen(i) / lambda) / normsToScreen(i);
		waveFunction(i,1) += std::sin(std::atan(1) * 2. * normsToScreen(i) / lambda) / normsToScreen(i);
	}
	return waveFunction;
}

void exportData(Eigen::MatrixXd waveFunction) {
	std::vector<std::vector<double>> waveFunctionVectors = {};
	for (int i = 0; i < waveFunction.rows(); i++) {
		Eigen::VectorXd waveFunctionAtScreenPosition(5);
		waveFunctionAtScreenPosition << waveFunction(i, 0), waveFunction(i, 1), waveFunction(i, 2), waveFunction(i, 3), waveFunction(i, 4);
		std::vector<double> finalWaveFunctionVector(waveFunctionAtScreenPosition.data(), waveFunctionAtScreenPosition.data() + waveFunctionAtScreenPosition.size());
		waveFunctionVectors.push_back(finalWaveFunctionVector);
	}
	
	std::ofstream out("C:\\Users\\Michael\\Documents\\Programming\\laueDiffractionResults\\EigenResults.csv");
	for (auto& row : waveFunctionVectors) {
		for (auto col : row)
			out << col << ',';
		out << '\n';
	}
}


int main() {
	
	// only change the values between here and the next comment line 
	
	nx = 1; // number of COPYS of unit cells in x/y/z directions (set to zero, you still maintain the base unit cell)
	ny = 1;
	nz = 1;
	lambda = 1e-9;
	debug = 0;
	wallXPosition = 0.01;
	wallLength = .2;
	dzdy = 0.0005;
	cellName = "salt"; //Current Options (salt, graphite) [this nwo sets the cellType ]
	cellCentering = "face-centered"; // Current Options:  {primative, body-centered, face-centered, base-centered} (if centering doesnt exist for type, it assumes primative)
	
	phi = 0.;
	theta = 0.;
	psi = 0.;
	
	// Get properties of individual crystal with the proper ortientation
	Crystal crystal; 
	crystal.setCellName(cellName);
	crystal.setCellType(cellName, cellCentering);
	crystal.setCellProperties(cellName, cellCentering, crystal.getCellType());
	std::cout << "Using cellName = " << cellName << " || Using cellType = " << crystal.getCellType() << " || Using cellCentering = " << cellCentering << std::endl;
	Eigen::MatrixXd cellStructure = getCellStructure(crystal.getCellType(), cellCentering, crystal.getAxialDistanceA(), crystal.getAxialDistanceB(), crystal.getAxialDistanceC(),
													crystal.getAxialAngleAlpha(), crystal.getAxialAngleBeta(), crystal.getAxialAngleGamma());
	// Build up the lattice and screen
	Eigen::MatrixXd fullLattice = buildLattice(nx, ny, nz, crystal.getAxialDistanceA(), crystal.getAxialDistanceB(), crystal.getAxialDistanceC(), cellStructure);
	Eigen::MatrixXd rotatedLattice = rotateLattice(fullLattice, phi,  theta,  psi);
	wallDivisions = static_cast<int>(std::floor(wallLength / dzdy));
	Eigen::MatrixXd wave = buildWave(getScreen());
	
	// For each atom, get the contribution from scattering to every screen position. This is the bulk of the calculation 
	for (int a = 0 ; a < rotatedLattice.rows(); a++) {
		std::cout << "Working on atom number " << a << std::endl;
		Eigen::VectorXd normsToScreen = getNormToScreenPosition(getScreen(), rotatedLattice, a);
		wave = propogateWave(lambda, normsToScreen, getScreen(), wave);

	}
	
	exportData(wave);


	return 0;
}	

