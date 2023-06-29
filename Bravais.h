#ifndef BRAVAIS_H
#define BRAVAIS_H

#include <string>
#include <vector>
#include <Eigen/Dense>

Eigen::MatrixXd getCellStructure(std::string cellType, std::string cellCentering, double a, double b, double c, double alpha, double beta, double gamma);  // https://en.wikipedia.org/wiki/Crystal_system

#endif