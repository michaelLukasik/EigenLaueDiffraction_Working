#ifndef CRYSTAL_H
#define CRYSTAL_H
#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <string>
#include <vector>
#include <Eigen/Dense>

class Crystal {
public:
	//Crystal();
	//virtual ~Crystal();

	double getAxialDistanceA() { return axialDistanceA; }
	double getAxialDistanceB() { return axialDistanceB; }
	double getAxialDistanceC() { return axialDistanceC; }

	double getAxialAngleAlpha() { return axialAngleAlpha; }
	double getAxialAngleBeta() { return axialAngleBeta; }
	double getAxialAngleGamma() { return axialAngleGamma; }

	void setAxialDistanceA(double _axialDistanceA) { this->axialDistanceA = _axialDistanceA;}
	void setAxialDistanceB(double _axialDistanceB) { this->axialDistanceB = _axialDistanceB; }
	void setAxialDistanceC(double _axialDistanceC) { this->axialDistanceC = _axialDistanceC; }

	void setAxialAngleAlpha(double _axialAngleAlpha) { this->axialAngleAlpha = _axialAngleAlpha; }
	void setAxialAngleBeta(double _axialAngleBeta) { this->axialAngleBeta = _axialAngleBeta; }
	void setAxialAngleGamma(double _axialAngleGamma) { this->axialAngleGamma = _axialAngleGamma; }



	Eigen::MatrixXd getCellStructure() { return cellStructure; }
	std::string getCallCentering() { return cellCentering; }
	std::string getCellName() { return cellName; }
	std::string getCellType() { return cellType; }

	void setCellStructure(Eigen::MatrixXd _cellStructure) { this->cellStructure = _cellStructure; }
	void setCellCentering(std::string _cellCentering) { this->cellCentering = _cellCentering; }
	void setCellName(std::string _cellName) { this->cellName = _cellName; }
	
	// The following are defined in crystal.cpp, requires look up based on inputs
	void setCellType(std::string _cellName, std::string _cellCentering);
	void setCellProperties(std::string _cellName, std::string _cellCentering, std::string _cellType);






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