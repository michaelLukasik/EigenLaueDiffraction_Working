#ifndef BRAVAIS_H
#define BRAVAIS_H

#include <Crystal.h>
#include <string>
#include <vector>
#include <Eigen/Dense>

Eigen::MatrixXd getCellStructure(const Crystal& crystal);   // https://en.wikipedia.org/wiki/Crystal_system

#endif