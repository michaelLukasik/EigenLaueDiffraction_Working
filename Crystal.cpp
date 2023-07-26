#include <Crystal.h>
#include <config.h>
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

void Crystal::setCellStrings(const Config& config) {
	Crystal::setCellCentering(config.getCellCentering());
	Crystal::setCellName(config.getCellName());
}


void Crystal::setCellType(Crystal& crystal) {
	
	const std::string cellCentering = crystal.getCellCentering();
	const std::string cellName = crystal.getCellName();

	std::string monoclinicArray[1] = {"jadeite"};
	std::string orthorhombicArray[1] = {"olivine"};
	std::string tetragonalArray[1] = {"pinnoite"};
	std::string rhombohedralArray[1] = {"dolomite"};
	std::string hexagonalArray[1] = { "graphite" };
	std::string cubicArray[1] = { "salt" };
	std::string uniqueArray[5] = { "twoCrystal", "loss", "DNA", "BDNA", "B10DNA" };


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
	else if (std::find(std::begin(uniqueArray), std::end(uniqueArray), cellName) != std::end(uniqueArray) && (cellCentering == "primative" || "body-centered" || "face-centered"))
	{
		this->cellType = "unique";
	}
	else
	{
		std::cout << "Cell name not recognized, Name: " << cellName << " and Centering: " << cellCentering << " is not a valid configuration, check again.";
		throw;
	}
}
void Crystal::setCellProperties(Crystal& crystal) {
//void Crystal::setCellProperties(std::string cellName, std::string cellCentering, std::string cellType) { //https://next-gen.materialsproject.org/materials/mp-22862
	
	const std::string cellName = crystal.getCellName();

	const std::string cellType = crystal.getCellType();
	const std::string cellCentering = crystal.getCellCentering();

	if (cellName == "graphite") { // https://som.web.cmu.edu/structures/S022-C-graphite.html
		Crystal::setAxialDistanceA(2.456e-10);
		Crystal::setAxialDistanceB(2.456e-10);
		Crystal::setAxialDistanceC(6.696e-10);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(120);
	}
	else if (cellName == "salt") { // https://next-gen.materialsproject.org/materials/mp-22862

		Crystal::setAxialDistanceA(5.59e-10);
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
	else if (cellName == "twoCrystal") { // https://www.mindat.org/min-1304.html

		Crystal::setAxialDistanceA(5.59e-10);
		Crystal::setAxialDistanceB(5.59e-10);
		Crystal::setAxialDistanceC(5.59e-10);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(90);
	}

	else if (cellName == "loss") { // :(
		Crystal::setAxialDistanceA(5.59e-10);
		Crystal::setAxialDistanceB(5.59e-10);
		Crystal::setAxialDistanceC(5.59e-10);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(90);
	}

	else if (cellName == "graphite2D") { // 
		Crystal::setAxialDistanceA(2.456e-10);
		Crystal::setAxialDistanceB(2.456e-10);
		Crystal::setAxialDistanceC(0);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(120);
	}
	else if (cellName == "DNA") { // https://pubs.aip.org/aapt/ajp/article/86/2/95/1057814/Rosalind-Franklin-s-X-ray-photo-of-DNA-as-an
		Crystal::setAxialDistanceA(1.e-9);
		Crystal::setAxialDistanceB(3.4e-9);
		Crystal::setAxialDistanceC(1.e-9);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(90);
	}
	else if (cellName == "BDNA") { // https://pubs.aip.org/aapt/ajp/article/86/2/95/1057814/Rosalind-Franklin-s-X-ray-photo-of-DNA-as-an
		Crystal::setAxialDistanceA(1.e-9);
		Crystal::setAxialDistanceB(3.4e-9);
		Crystal::setAxialDistanceC(1.e-9);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(90);
	}
	else if (cellName == "B10DNA") { // https://pubs.aip.org/aapt/ajp/article/86/2/95/1057814/Rosalind-Franklin-s-X-ray-photo-of-DNA-as-an
		Crystal::setAxialDistanceA(1.e-9);
		Crystal::setAxialDistanceB(3.4e-9);
		Crystal::setAxialDistanceC(1.e-9);

		Crystal::setAxialAngleAlpha(90);
		Crystal::setAxialAngleBeta(90);
		Crystal::setAxialAngleGamma(90);
	}
	else {
		throw "No match to cellName, check cellName and rerun the program";
	}

}
