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

std::string validCrystals[2] = { "salt", "graphite" };
std::string validCenterings[2] = { "primative", "body-centered" };

void Crystal::setCellType(std::string cellName, std::string cellCentering) {
	this->cellCentering = cellCentering;
	std::string monoclinicArray[1] = {"jadeite"};
	std::string orthorhombicArray[1] = {"olivine"};
	std::string tetragonalArray[1] = {"pinnoite"};
	std::string rhombohedralArray[1] = {"dolomite"};
	std::string hexagonalArray[1] = { "graphite" };
	std::string cubicArray[1] = { "salt" };

	if (std::find(std::begin(monoclinicArray), std::end(monoclinicArray), cellName) != std::end(monoclinicArray) && (cellCentering == "primative" || "base-centered")) {
		
		this->cellType = "monoclinic";
	}
	else if (std::find(std::begin(orthorhombicArray), std::end(orthorhombicArray), cellName) != std::end(orthorhombicArray) && (cellCentering == "primative" || "base-centered" || "body-centered" || "face-centered")) {

		this->cellType = "orthorhombic";
	}
	else if (std::find(std::begin(tetragonalArray), std::end(tetragonalArray), cellName) != std::end(tetragonalArray) && (cellCentering == "primative" || "body-centered")) {

		this->cellType = "tetragonal";
	}
	else if (std::find(std::begin(rhombohedralArray), std::end(rhombohedralArray), cellName) != std::end(rhombohedralArray) && (cellCentering == "primative")) {

		this->cellType = "rhombohedral";
	}
	else if (std::find(std::begin(hexagonalArray), std::end(hexagonalArray), cellName) != std::end(hexagonalArray) && (cellCentering == "primative"))
	{
		this->cellType = "hexagonal";
	}
	else if (std::find(std::begin(cubicArray), std::end(cubicArray), cellName) != std::end(cubicArray) && (cellCentering  == "primative" || "body-centered" || "face-centered"))
	{
		this->cellType = "cubic";
	}
	else
	{
		std::cout << "Cell name not recognized, Name: " << cellName << " and Centering: " << cellCentering << " is not a valid configuration, check again.";
		throw;
	}
}

void Crystal::setCellProperties(std::string cellName, std::string cellCentering, std::string cellType) { //https://next-gen.materialsproject.org/materials/mp-22862
	if (cellName == "graphite") { // https://som.web.cmu.edu/structures/S022-C-graphite.html

		Crystal::setAxialDistanceA(2.456e-10);
		Crystal::setAxialDistanceB(2.456e-10);
		Crystal::setAxialDistanceC(6.696e-10);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(120);
	}
	else if (cellName == "salt") { // https://next-gen.materialsproject.org/materials/mp-22862

		Crystal::setAxialDistanceA(pow(20, -5));//5.59e-10);
		Crystal::setAxialDistanceB(5.59e-10);
		Crystal::setAxialDistanceC(5.59e-10);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(90);
	}
	else if (cellName == "jadeite") { // https://www.mindat.org/min-2062.html

		Crystal::setAxialDistanceA(9.418e-10);
		Crystal::setAxialDistanceB(8562e-10);
		Crystal::setAxialDistanceC(5.219e-10);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(107.58);
		Crystal::setAxialAngleGamma(90);
	}
	else if (cellName == "olivine") { // https://www.webmineral.com/data/Olivine.shtml

		Crystal::setAxialDistanceA(4.78e-10);
		Crystal::setAxialDistanceB(10.25e-10);
		Crystal::setAxialDistanceC(6.3e-10);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(90);
	}
	else if (cellName == "pinnoite") { // https://www.mindat.org/min-3217.html

		Crystal::setAxialDistanceA(7.617e-10);
		Crystal::setAxialDistanceB(7.617e-10);
		Crystal::setAxialDistanceC(8.19e-10);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(90);
	}
	else if (cellName == "dolomite") { // https://www.mindat.org/min-1304.html

		Crystal::setAxialDistanceA(4.815e-10);
		Crystal::setAxialDistanceB(4.815e-10);
		Crystal::setAxialDistanceC(16.119e-10);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(120);
	}
	else {
		throw "No match to cellName, check cellName and rerun the program";
	}

}