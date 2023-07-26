#include <Bravais.h>
#include <math.h>
#include <config.h>
#include <Crystal.h>
#include <iostream>
#include <complex>
#include <numeric>
#include <Eigen\Dense>
Eigen::MatrixXd getCellStructure(const Crystal& crystal){
	
	std::string cellType = crystal.getCellType();
	std::string cellCentering = crystal.getCellCentering();
	std::string cellName = crystal.getCellName();
	
	double a = crystal.getAxialDistanceA();
	double b = crystal.getAxialDistanceB();
	double c = crystal.getAxialDistanceC();


	double pi = 3.14159265358979323846;
	double alpha = crystal.getAxialAngleAlpha()* (pi/180.);
	double beta = crystal.getAxialAngleBeta() * (pi / 180.);
	double gamma = crystal.getAxialAngleGamma() * (pi / 180.);

	//Eigen::MatrixXd cellStructure;

	//Fix the general Case!!
	/*if (cellType == "triclinic") { // 3 different atomic spacing, 3 different axial angles
		Eigen::MatrixXd cellStructure(8, 3);
		cellStructure << 0, 0, 0 ;
		cellStructure << a*std::sin(gamma),a*std::cos(gamma),

	}
	*/

	if (cellType == "monoclinic") {
		if (cellCentering != "primative" or cellCentering != "base-centered") {
			cellCentering = "primative";
		}
		if (cellCentering == "primative") {
			Eigen::MatrixXd cellStructure(8, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << 0, b* std::sin(beta), c* std::cos(beta);
			cellStructure.row(3) << a, b* std::sin(beta), c* std::cos(beta);
			cellStructure.row(4) << 0, 0, c* std::sin(beta);
			cellStructure.row(5) << a, 0, c* std::sin(beta) + b * std::cos(beta);
			cellStructure.row(6) << 0, b* std::sin(beta), (c * std::sin(beta)) + b * std::cos(beta);
			cellStructure.row(7) << a, b* std::sin(beta), (c * std::sin(beta)) + b * std::cos(beta);
			return cellStructure;
		}
		else if (cellCentering == "base-centered") {
			Eigen::MatrixXd cellStructure(14, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << 0, b* std::sin(beta), c* std::cos(beta);
			cellStructure.row(3) << a, b* std::sin(beta), c* std::cos(beta);
			cellStructure.row(4) << 0, 0, c* std::sin(beta);
			cellStructure.row(5) << a, 0, c* std::sin(beta) + b * std::cos(beta);
			cellStructure.row(6) << 0, b* std::sin(beta), c* std::sin(beta) + b * std::cos(beta);
			cellStructure.row(7) << a, b* std::sin(beta), c* std::sin(beta) + b * std::cos(beta);
			cellStructure.row(8) << a / 2, b* std::sin(beta) / 2, (c * std::cos(beta) / 2);
			cellStructure.row(9) << a / 2, b* std::sin(beta) / 2, (c * std::cos(beta) / 2) + (b * std::cos(beta) / 2);
			return cellStructure;
		}
	}
	else if (cellType == "orthorhombic") { //
		if (cellCentering == "primative") {
			Eigen::MatrixXd cellStructure(8, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << 0, b, 0;
			cellStructure.row(3) << 0, 0, c;
			cellStructure.row(4) << a, b, 0;
			cellStructure.row(5) << a, 0, c;
			cellStructure.row(6) << 0, b, c;
			cellStructure.row(7) << a, b, c;
			return cellStructure;
		}
		if (cellCentering == "base-centered") {
			Eigen::MatrixXd cellStructure(10, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << 0, b, 0;
			cellStructure.row(3) << 0, 0, c;
			cellStructure.row(4) << a, b, 0;
			cellStructure.row(5) << a, 0, c;
			cellStructure.row(6) << 0, b, c;
			cellStructure.row(7) << a, b, c;
			cellStructure.row(8) << a / 2., b / 2., 0;
			cellStructure.row(9) << a / 2., b / 2., c;
			return cellStructure;
		}
		if (cellCentering == "body-centered") {
			Eigen::MatrixXd cellStructure(9, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << 0, b, 0;
			cellStructure.row(3) << 0, 0, c;
			cellStructure.row(4) << a, b, 0;
			cellStructure.row(5) << a, 0, c;
			cellStructure.row(6) << 0, b, c;
			cellStructure.row(7) << a, b, c;
			cellStructure.row(8) << a / 2., b / 2., c / 2.;
			return cellStructure;
		}
		if (cellCentering == "face-centered") {
			Eigen::MatrixXd cellStructure(14, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << 0, b, 0;
			cellStructure.row(3) << 0, 0, c;
			cellStructure.row(4) << a, b, 0;
			cellStructure.row(5) << a, 0, c;
			cellStructure.row(6) << 0, b, c;
			cellStructure.row(7) << a, b, c;
			cellStructure.row(8) << a / 2., b / 2., 0;
			cellStructure.row(9) << a / 2., b / 2., c / 2;
			cellStructure.row(10) << 0., b / 2., c / 2;
			cellStructure.row(11) << a, b / 2., c / 2;
			cellStructure.row(12) << a / 2., 0., c / 2.;
			cellStructure.row(13) << a / 2., b, c / 2;
			return cellStructure;
		}
		else {
			std::cout << "WARNING: centering not recognized for this cellType \n";
		}
	}
	else if (cellType == "tetragonal") {
		if (cellCentering != "primative" or cellCentering != "body-centered") {
				std::cout << "WARNING: centering not recognized for this cellType, using primative case instead \n";
				cellCentering = "primative";
			}
		if (cellCentering == "primative") {
			Eigen::MatrixXd cellStructure(8, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << 0, a, 0;
			cellStructure.row(3) << 0, 0, c;
			cellStructure.row(4) << a, a, 0;
			cellStructure.row(5) << a, 0, c;
			cellStructure.row(6) << 0, a, c;
			cellStructure.row(7) << a, a, c;
			return cellStructure;
		}
		if (cellCentering == "body-centered") {
			Eigen::MatrixXd cellStructure(9, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << 0, a, 0;
			cellStructure.row(3) << 0, 0, c;
			cellStructure.row(4) << a, a, 0;
			cellStructure.row(5) << a, 0, c;
			cellStructure.row(6) << 0, a, c;
			cellStructure.row(7) << a, a, c;
			cellStructure.row(8) << a / 2., a / 2., c / 2.;
			return cellStructure;
		}

		else {
				std::cout << "WARNING: centering not recognized for this cellType \n";
		}
	}
	else if (cellType == "rhombohedral") {
		if (cellCentering != "primative") {
				std::cout << "WARNING: centering not recognized for this cellType, using primative case instead [only primative centering exists for this cell] \n";
				cellCentering = "primative";
		}
		Eigen::MatrixXd cellStructure(8, 3);
		cellStructure.row(0) << 0, 0, 0;
		cellStructure.row(1) << a * std::sin(alpha), a* std::cos(alpha), a* std::cos(alpha);
		cellStructure.row(2) << a * std::cos(alpha), a* std::sin(alpha), a* std::cos(alpha);
		cellStructure.row(3) << a * std::sin(alpha), a* std::sin(alpha), a* std::cos(alpha);
		cellStructure.row(4) << a * std::cos(alpha), a* std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha);
		cellStructure.row(5) << a * std::sin(alpha) + a * std::cos(alpha), a* std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha);
		cellStructure.row(6) << a * std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha);
		cellStructure.row(7) << a * std::sin(alpha) + a * std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha);
		return cellStructure;

	}
	else if (cellType == "hexagonal") {
		if (cellCentering != "primative") {
				std::cout << "WARNING: centering not recognized for this cellType, using primative case instead [only primative centering exists for this cell] \n";
				cellCentering = "primative";
		}
		Eigen::MatrixXd cellStructure(8, 3);
		cellStructure.row(0) << 0, 0, 0;
		cellStructure.row(1) << a, 0, 0;
		cellStructure.row(2) << -a * std::sin(pi / 6.), a* std::cos(pi / 6.), 0;
		cellStructure.row(3) << a - a * std::sin(pi / 6.), a* std::cos(pi / 6.), 0;
		cellStructure.row(4) << 0, 0, c;
		cellStructure.row(5) << a, 0, c;
		cellStructure.row(6) << -a * std::sin(pi / 6.), a* std::cos(pi / 6.), c;
		cellStructure.row(7) << a - a * std::sin(pi / 6.), a* std::cos(pi / 6.), c;
		return cellStructure;

	}
	else if (cellType == "cubic") {
		if (cellCentering == "primative") {
			Eigen::MatrixXd cellStructure(8, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << 0, a, 0;
			cellStructure.row(3) << 0, 0, a;
			cellStructure.row(4) << a, a, 0;
			cellStructure.row(5) << a, 0, a;
			cellStructure.row(6) << 0, a, a;
			cellStructure.row(7) << a, a, a;
			return cellStructure;
		}

		else if (cellCentering == "body-centered") {
			Eigen::MatrixXd cellStructure(9, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << 0, a, 0;
			cellStructure.row(3) << 0, 0, a;
			cellStructure.row(4) << a, a, 0;
			cellStructure.row(5) << a, 0, a;
			cellStructure.row(6) << 0, a, a;
			cellStructure.row(7) << a, a, a;
			cellStructure.row(8) << a / 2., a / 2., a / 2.;
			return cellStructure;

		}
		else if (cellCentering == "face-centered") {
			Eigen::MatrixXd cellStructure(14, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << 0, a, 0;
			cellStructure.row(3) << 0, 0, a;
			cellStructure.row(4) << a, a, 0;
			cellStructure.row(5) << a, 0, a;
			cellStructure.row(6) << 0, a, a;
			cellStructure.row(7) << a, a, a;
			cellStructure.row(8) << 0, a / 2, a / 2;
			cellStructure.row(9) << a, a / 2, a / 2;
			cellStructure.row(10) << a / 2, 0, a / 2;
			cellStructure.row(11) << a / 2, a, a / 2;
			cellStructure.row(12) << a / 2, a / 2, 0;
			cellStructure.row(13) << a / 2, a / 2, a;
			return cellStructure;

		}
		else {
				std::cout << "WARNING: centering not recognized for this cellType \n";
		}
	}

	else if (cellType == "unique") {
		if (cellName == "twoCrystal") {
			if (cellCentering == "primative") {
				Eigen::MatrixXd cellStructure(1, 3);
				cellStructure.row(0) << 0, 0, 0;
				return cellStructure;
			}
			else {
				std::cout << "The two crystal testing objet can only be set to centering = primative, any other centering is nonsense. Fix and rerun. ";
				throw;
			}
		}
		else if (cellName == "loss") {
			Eigen::MatrixXd cellStructure(20, 3);
			//  |
			cellStructure.row(0) << 0, -2. * a, 3. * a;
			cellStructure.row(1) << 0, -2. * a, 2. * a;
			cellStructure.row(2) << 0, -2. * a, 1. * a;
			//  ||
			cellStructure.row(3) << 0, a, 3. * a;
			cellStructure.row(4) << 0, a, 2. * a;
			cellStructure.row(5) << 0, a, 1. * a;
			cellStructure.row(6) << 0, 2. * a, 2. * a;
			cellStructure.row(7) << 0, 2. * a, a;
			//  ||
			cellStructure.row(8) << 0, -2. * a, -3. * a;
			cellStructure.row(9) << 0, -2. * a, -2. * a;
			cellStructure.row(10) << 0, -2. * a, -1. * a;
			cellStructure.row(11) << 0, -1. * a, -3. * a;
			cellStructure.row(12) << 0, -1. * a, -2. * a;
			cellStructure.row(13) << 0, -1. * a, -1. * a;
			// |_
			cellStructure.row(14) << 0, a, -1. * a;
			cellStructure.row(15) << 0, a, -2. * a;
			cellStructure.row(16) << 0, a, -3. * a;
			cellStructure.row(17) << 0, 1.5 * a, -2.5 * a;
			cellStructure.row(18) << 0, 2.5 * a, -2.5 * a;
			cellStructure.row(19) << 0, 3.5 * a, -2.5 * a;
			return cellStructure;
		}
		else if (cellName == "graphite2D") {
			Eigen::MatrixXd cellStructure(4, 3);
			cellStructure.row(0) << 0, 0, 0;
			cellStructure.row(1) << a, 0, 0;
			cellStructure.row(2) << -a * std::sin(pi / 6.), a* std::cos(pi / 6.), 0;
			cellStructure.row(3) << a - a * std::sin(pi / 6.), a* std::cos(pi / 6.), 0;
			return cellStructure;
		}
		else if (cellName == "DNA") {
			const float pitch = 3.4e-9;
			const float r = 1.0e-9;
			const float inner_r = r / 2.;
			const int nPoints = 10;
			const float nPointsf = float(nPoints);
			const float phase = (180.) * pi / 180.;

			Eigen::MatrixXd cellStructure(4*nPoints, 3); // Pointing Up Towards +y
			for (int i = 0; i < nPoints; ++i) {
				cellStructure.row(i) << r * std::sin( (float(i) * 2. * pi) / nPointsf), i* (pitch / nPointsf), r * std::cos((float(i) * 2. * pi) / nPointsf);
				cellStructure.row(i + nPoints) << r * std::sin(-1.*((float(i) * 2. * pi + phase ) / nPointsf) + phase), i* (pitch / nPointsf), r* std::cos(-1.*((float(i) * 2. * pi) / nPointsf) + phase);
				cellStructure.row(i + (2*nPoints)) << inner_r * std::sin((float(i) * 2. * pi) / nPointsf), i* (pitch / nPointsf), inner_r* std::cos((float(i) * 2. * pi) / nPointsf);
				cellStructure.row(i + (3*nPoints)) << inner_r * std::sin(((float(i) * 2. * pi) / nPointsf) + phase), i* (pitch / nPointsf), inner_r* std::cos(((float(i) * 2. * pi) / nPointsf) + phase);
			}
			for (int i = 0; i < cellStructure.rows(); ++i) {
				std::cout << "[" << cellStructure(i, 0) << "," << cellStructure(i, 1) << "," << cellStructure(i, 2) << "]," << std::endl;
			}
			return cellStructure;
		}

		else if (cellName == "BDNA") {
			const float pitch = 3.4e-9;
			const float r = 1.0e-9;
			const float inner_r = r / 2.;
			const float inner_inner_r = r / 4.;
			const int nPoints = 20;
			const float nPointsf = float(nPoints);
			const float phase = (90.) * pi / 180.;

			Eigen::MatrixXd cellStructure(6 * nPoints, 3); // Pointing Up Towards +y
			for (int i = 0; i < nPoints; ++i) {
				cellStructure.row(i) << r * std::sin((float(i) * 2. * pi) / nPointsf), i* (pitch / nPointsf), r* std::cos((float(i) * 2. * pi) / nPointsf);
				cellStructure.row(i + nPoints) << r * std::sin(-1. * ((float(i) * 2. * pi + phase) / nPointsf) + phase), i* (pitch / nPointsf), r* std::cos(-1. * ((float(i) * 2. * pi) / nPointsf) + phase);
				cellStructure.row(i + (2 * nPoints)) << inner_r * std::sin((float(i) * 2. * pi) / nPointsf), i* (pitch / nPointsf), inner_r* std::cos((float(i) * 2. * pi) / nPointsf);
				cellStructure.row(i + (3 * nPoints)) << inner_r * std::sin(((float(i) * 2. * pi) / nPointsf) + phase), i* (pitch / nPointsf), inner_r* std::cos(((float(i) * 2. * pi) / nPointsf) + phase);
				cellStructure.row(i + (4 * nPoints)) << inner_inner_r * std::sin((float(i) * 2. * pi) / nPointsf), i* (pitch / nPointsf), inner_inner_r* std::cos((float(i) * 2. * pi) / nPointsf);
				cellStructure.row(i + (5 * nPoints)) << inner_inner_r * std::sin(((float(i) * 2. * pi) / nPointsf) + phase), i* (pitch / nPointsf), inner_inner_r* std::cos(((float(i) * 2. * pi) / nPointsf) + phase);
			}
			//for (int i = 0; i < cellStructure.rows(); ++i) {
			//	std::cout << "[" << cellStructure(i, 0) << "," << cellStructure(i, 1) << "," << cellStructure(i, 2) << "]," << std::endl;
			//}
			//std::cin.ignore();
			return cellStructure;
		}

		else if (cellName == "B10DNA") {
			const float pitch = 3.4e-9;
			const float r = 1.0e-9;
			const float inner_r = r / 2.;
			const float inner_inner_r = r / 4.;
			const int nPoints = 10;
			const float nPointsf = float(nPoints);
			const float phase = (180.) * pi / 180.;

			Eigen::MatrixXd cellStructure(6 * nPoints, 3); // Pointing Up Towards +y
			for (int i = 0; i < nPoints; ++i) {
				cellStructure.row(i) << r * std::sin((float(i) * 2. * pi) / nPointsf), i* (pitch / nPointsf), r* std::cos((float(i) * 2. * pi) / nPointsf);
				cellStructure.row(i + nPoints) << r * std::sin(-1. * ((float(i) * 2. * pi + phase) / nPointsf) + phase), i* (pitch / nPointsf), r* std::cos(-1. * ((float(i) * 2. * pi) / nPointsf) + phase);
				cellStructure.row(i + (2 * nPoints)) << inner_r * std::sin((float(i) * 2. * pi) / nPointsf), i* (pitch / nPointsf), inner_r* std::cos((float(i) * 2. * pi) / nPointsf);
				cellStructure.row(i + (3 * nPoints)) << inner_r * std::sin(((float(i) * 2. * pi) / nPointsf) + phase), i* (pitch / nPointsf), inner_r* std::cos(((float(i) * 2. * pi) / nPointsf) + phase);
				cellStructure.row(i + (4 * nPoints)) << inner_inner_r * std::sin((float(i) * 2. * pi) / nPointsf), i* (pitch / nPointsf), inner_inner_r* std::cos((float(i) * 2. * pi) / nPointsf);
				cellStructure.row(i + (5 * nPoints)) << inner_inner_r * std::sin(((float(i) * 2. * pi) / nPointsf) + phase), i* (pitch / nPointsf), inner_inner_r* std::cos(((float(i) * 2. * pi) / nPointsf) + phase);
			}
			//for (int i = 0; i < cellStructure.rows(); ++i) {
			//	std::cout << "[" << cellStructure(i, 0) << "," << cellStructure(i, 1) << "," << cellStructure(i, 2) << "]," << std::endl;
			//}
			//std::cin.ignore();
			return cellStructure;
		}
		else {
			std::cout << "WARNING: cellType not recognized, returing empty cell structure \n";
		}

	}

}