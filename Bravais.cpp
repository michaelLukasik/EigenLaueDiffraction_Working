#include <Bravais.h>
#include <math.h>
#include <iostream>
#include <complex>
#include <numeric>
#include <Eigen\Dense>

Eigen::MatrixXd getCellStructure(std::string cellType, std::string cellCentering, double a, double b, double c, double alphadeg, double betadeg, double gammadeg) { // https://en.wikipedia.org/wiki/Crystal_system
	double pi = 3.14159265358979323846;
	double alpha = alphadeg* (pi/180.);
	double beta = betadeg * (pi / 180.);
	double gamma = gammadeg * (pi / 180.);
	Eigen::MatrixXd cellStructure;

	//Fix the general Case!!
	/*if (cellType == "triclinic") { // 3 different atomic spacing, 3 different axial angles
		Eigen::MatrixXd cellStructure(8, 3);
		cellStructure << 0, 0, 0 ;
		cellStructure << a*std::sin(gamma),a*std::cos(gamma),

	}
	*/

	if (cellType == "monoclinic") {
		if (cellCentering != "primative" or cellCentering != "base-centered") {
			std::cout << "WARNING: centering not recognized for this cellType, using primative case instead \n";
			cellCentering = "primative";
		}
		if (cellCentering == "primative") {
			Eigen::MatrixXd cellStructure(8, 3);
			cellStructure << 0, 0, 0;
			cellStructure << a, 0, 0;
			cellStructure << 0, b* std::sin(beta), c* std::cos(beta);
			cellStructure << a, b* std::sin(beta), c* std::cos(beta);
			cellStructure << 0, 0, c* std::sin(beta);
			cellStructure << a, 0, c* std::sin(beta) + b * std::cos(beta);
			cellStructure << 0, b* std::sin(beta), (c * std::sin(beta)) + b * std::cos(beta);
			cellStructure << a, b* std::sin(beta), (c * std::sin(beta)) + b * std::cos(beta);
		}
		else if (cellCentering == "base-centered") {
			Eigen::MatrixXd cellStructure(14, 3);
			cellStructure << 0, 0, 0;
			cellStructure << a, 0, 0;
			cellStructure << 0, b* std::sin(beta), c* std::cos(beta);
			cellStructure << a, b* std::sin(beta), c* std::cos(beta);
			cellStructure << 0, 0, c* std::sin(beta);
			cellStructure << a, 0, c* std::sin(beta) + b * std::cos(beta);
			cellStructure << 0, b* std::sin(beta), c* std::sin(beta) + b * std::cos(beta);
			cellStructure << a, b* std::sin(beta), c* std::sin(beta) + b * std::cos(beta);
			cellStructure << a / 2, b* std::sin(beta) / 2, (c * std::cos(beta) / 2);
			cellStructure << a / 2, b* std::sin(beta) / 2, (c * std::cos(beta) / 2) + (b * std::cos(beta) / 2);
		}
	}
	else if (cellType == "orthorhombic") { //
		if (cellCentering == "primative") {
			Eigen::MatrixXd cellStructure(8, 3);
			cellStructure << 0, 0, 0;
			cellStructure << a, 0, 0;
			cellStructure << 0, b, 0;
			cellStructure << 0, 0, c;
			cellStructure << a, b, 0;
			cellStructure << a, 0, c;
			cellStructure << 0, b, c;
			cellStructure << a, b, c;
			}
		if (cellCentering == "base-centered") {
			Eigen::MatrixXd cellStructure(10, 3);
			cellStructure << 0, 0, 0;
			cellStructure << a, 0, 0;
			cellStructure << 0, b, 0;
			cellStructure << 0, 0, c;
			cellStructure << a, b, 0;
			cellStructure << a, 0, c;
			cellStructure << 0, b, c;
			cellStructure << a, b, c;
			cellStructure << a / 2., b / 2., 0;
			cellStructure << a / 2., b / 2., c;
			}
		if (cellCentering == "body-centered") {
			Eigen::MatrixXd cellStructure(9, 3);
			cellStructure << 0, 0, 0;
			cellStructure << a, 0, 0;
			cellStructure << 0, b, 0;
			cellStructure << 0, 0, c;
			cellStructure << a, b, 0;
			cellStructure << a, 0, c;
			cellStructure << 0, b, c;
			cellStructure << a, b, c;
			cellStructure << a / 2., b / 2., c / 2.;
			}
		if (cellCentering == "face-centered") {
			Eigen::MatrixXd cellStructure(14, 3);
			cellStructure << 0, 0, 0;
			cellStructure << a, 0, 0;
			cellStructure << 0, b, 0;
			cellStructure << 0, 0, c;
			cellStructure << a, b, 0;
			cellStructure << a, 0, c;
			cellStructure << 0, b, c;
			cellStructure << a, b, c;
			cellStructure << a / 2., b / 2., 0;
			cellStructure << a / 2., b / 2., c / 2;
			cellStructure << 0., b / 2., c / 2;
			cellStructure << a, b / 2., c / 2;
			cellStructure << a / 2., 0., c / 2.;
			cellStructure << a / 2., b, c / 2;
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
			cellStructure << 0, 0, 0;
			cellStructure << a, 0, 0;
			cellStructure << 0, a, 0;
			cellStructure << 0, 0, c;
			cellStructure << a, a, 0;
			cellStructure << a, 0, c;
			cellStructure << 0, a, c;
			cellStructure << a, a, c;
			}
		if (cellCentering == "body-centered") {
			Eigen::MatrixXd cellStructure(9, 3);
			cellStructure << 0, 0, 0;
			cellStructure << a, 0, 0;
			cellStructure << 0, a, 0;
			cellStructure << 0, 0, c;
			cellStructure << a, a, 0;
			cellStructure << a, 0, c;
			cellStructure << 0, a, c;
			cellStructure << a, a, c;
			cellStructure << a / 2., a / 2., c / 2.;
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
		cellStructure << 0, 0, 0;
		cellStructure << a * std::sin(alpha), a* std::cos(alpha), a* std::cos(alpha);
		cellStructure << a * std::cos(alpha), a* std::sin(alpha), a* std::cos(alpha);
		cellStructure << a * std::sin(alpha), a* std::sin(alpha), a* std::cos(alpha);
		cellStructure << a * std::cos(alpha), a* std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha);
		cellStructure << a * std::sin(alpha) + a * std::cos(alpha), a* std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha);
		cellStructure << a * std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha);
		cellStructure << a * std::sin(alpha) + a * std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha), a* std::sin(alpha) + a * std::cos(alpha);


	}
	else if (cellType == "hexagonal") {
		if (cellCentering != "primative") {
				std::cout << "WARNING: centering not recognized for this cellType, using primative case instead [only primative centering exists for this cell] \n";
				cellCentering = "primative";
		}
		Eigen::MatrixXd cellStructure(8, 3);
		cellStructure << 0, 0, 0;
		cellStructure << a, 0, 0;
		cellStructure << -a * std::sin(pi / 6.), a* std::cos(pi / 6.), 0;
		cellStructure << a - a * std::sin(pi / 6.), a* std::cos(pi / 6.), 0;
		cellStructure << 0, 0, c;
		cellStructure << a, 0, c;
		cellStructure << -a * std::sin(pi / 6.), a* std::cos(pi / 6.), c;
		cellStructure << a - a * std::sin(pi / 6.), a* std::cos(pi / 6.), c;

	}
	else if (cellType == "cubic") {
		//std::cout << " Inside salt, primative " << std::endl;

		if (cellCentering == "primative") {
			//std::cout << " Inside salt, primative " << std::endl;
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
	else {
		std::cout << "WARNING: cellType not recognized, returing empty cell structure \n";
	}
	
	//return cellStructure;	
}


