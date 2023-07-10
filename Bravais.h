#ifndef BRAVAIS_H
#define BRAVAIS_H

#include <Crystal.h>
#include <string>
#include <vector>
#include <Eigen/Dense>

Eigen::MatrixXd getCellStructure(const Crystal& crystal);
//Eigen::MatrixXd getCellStructure(std::string cellType, std::string cellCentering, double a, double b, double c, double alpha, double beta, double gamma);  // https://en.wikipedia.org/wiki/Crystal_system

#endif