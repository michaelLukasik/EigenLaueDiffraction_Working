#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <config.h>
#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <string>
#include <vector>
#include <Eigen/Dense>

class Crystal {
public:
	//Crystal();
	//virtual ~Crystal();

	const double getAxialDistanceA() const { return axialDistanceA; }
	const double getAxialDistanceB() const { return axialDistanceB; }
	const double getAxialDistanceC() const { return axialDistanceC; }

	const double getAxialAngleAlpha() const { return axialAngleAlpha; }
	const double getAxialAngleBeta() const { return axialAngleBeta; }
	const double getAxialAngleGamma() const { return axialAngleGamma; }

	void setAxialDistanceA(double _axialDistanceA) { this->axialDistanceA = _axialDistanceA;}
	void setAxialDistanceB(double _axialDistanceB) { this->axialDistanceB = _axialDistanceB; }
	void setAxialDistanceC(double _axialDistanceC) { this->axialDistanceC = _axialDistanceC; }

	void setAxialAngleAlpha(double _axialAngleAlpha) { this->axialAngleAlpha = _axialAngleAlpha; }
	void setAxialAngleBeta(double _axialAngleBeta) { this->axialAngleBeta = _axialAngleBeta; }
	void setAxialAngleGamma(double _axialAngleGamma) { this->axialAngleGamma = _axialAngleGamma; }

	void setCellCentering(std::string _cellCentering) { this->cellCentering = _cellCentering; }
	void setCellName(std::string _cellName) { this->cellName = _cellName; }


	const Eigen::MatrixXd getCellStructure() const { return cellStructure; }
	const std::string getCellCentering() const { return cellCentering; }
	const std::string getCellName() const { return cellName; }
	const std::string getCellType() const{ return cellType; }

	void setCellStructure(Eigen::MatrixXd _cellStructure) { this->cellStructure = _cellStructure; }

	// The following are defined in crystal.cpp, requires look up based on inputs
	void setCellStrings(const Config& config);
	void setCellType(Crystal& crystal);
	void setCellProperties(Crystal& crystal);






private:
	std::string cellName;
	std::string cellCentering;
	std::string cellType;
	Eigen::MatrixXd cellStructure;

	double axialDistanceA;
	double axialDistanceB;
	double axialDistanceC;

	double axialAngleAlpha;
	double axialAngleBeta;
	double axialAngleGamma;


};

#endif